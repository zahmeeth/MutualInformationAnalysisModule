import time
import json
import os
import errno
import uuid
import csv
import math
import pandas as pd
import numpy as np
import collections
import natsort
import shutil
import itertools
import operator
from itertools import islice
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import OrderedDict
from itertools import groupby

from Workspace.WorkspaceClient import Workspace as Workspace
from DataFileUtil.DataFileUtilClient import DataFileUtil
from KBaseReport.KBaseReportClient import KBaseReport

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class MutualInfoUtil:
    def __init__(self, config):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.shock_url = config['shock-url']
        self.dfu = DataFileUtil(self.callback_url)
        self.ws = Workspace(self.ws_url, token=self.token)
        self.scratch = config['scratch']

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def _validate_run_flux_mutual_information_analysis_params(self, params):
        """
        _validate_run_flux_mutual_information_analysis_params:
                validates params passed to run_flux_mutual_information_analysis method
        """

        log('start validating validate_run_flux_mutual_information_analysis params')

        # check for required parameters
        for p in ['fbamodel_id', 'compounds', 'media_id', 'workspace_name']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))
    

    def _get_file_from_ws(workspace, obj_name):
        try:
            file_path = ws_client.get_objects(
                [{'name': obj_name,
                  'workspace': workspace}])[0]
        except Exception as e:
            raise ValueError(
                'Unable to get object from workspace: (' +
                workspace + '/' + obj_name + ')' + str(e))
        return file_path       

    def _generate_html_report(self, result_directory, mutual_info_dict):
        
        """
        _generate_html_report: generate html summary report
        """
        log('start generating html report')
        
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'mutual_information_report.html')
        
        shutil.copy(os.path.join(result_directory, 'MI_plot.png'),
                    os.path.join(output_directory, 'MI_plot.png'))
        
        overview_content = ''
        overview_content += '<table><tr><th>Mutual Information for various chemical compound combinations'
        overview_content += ' Object</th></td>'
        overview_content += '<tr><th>Input Chemical Compound Combination</th>'
        overview_content += '<th>Mutual Information (in Bits)</th>'
        overview_content += '</tr>'
        for k, v in mutual_info_dict.items():
            overview_content += '<tr><td>{}</td><td>{}</td></tr>'.format(k, v)
        overview_content += '</table>'

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          overview_content)
                result_file.write(report_template)

        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for Mutual Information App'})
        
        return html_report

    def _generate_output_file_list(self, result_directory, params):
        """
        _generate_output_file_list: zip result files and generate file_links for report
        """

        log('start packing result files')
        output_files = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file = os.path.join(output_directory, 'DESeq2_result.zip')
        plot_file = os.path.join(output_directory, 'DESeq2_plot.zip')

        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:
            for root, dirs, files in os.walk(result_directory):
                for file in files:
                    if not (file.endswith('.zip') or
                            file.endswith('.png') or
                            file.endswith('.DS_Store')):
                        zip_file.write(os.path.join(root, file), file)

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'File(s) generated by DESeq2 App'})

        with zipfile.ZipFile(plot_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:
            for root, dirs, files in os.walk(result_directory):
                for file in files:
                    if file.endswith('.png'):
                        zip_file.write(os.path.join(root, file), file)

        output_files.append({'path': plot_file,
                             'name': os.path.basename(plot_file),
                             'label': os.path.basename(plot_file),
                             'description': 'Visualization plots by DESeq2 App'})

        return output_files

    def _generate_report(self, result_directory, mutual_info_dict, params):
        """
        _generate_report: generate summary report
        """

        log('creating report')
        output_html_files = self._generate_html_report(result_directory,
                                                       mutual_info_dict)
        report_params = {'message': '',
                         'workspace_name': params.get('workspace_name'),
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 333,
                         'report_object_name': 'MutualInfomation_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output       

    def _generate_mutual_info(self, compounds_file, fba_file):

        df1 = pd.read_csv(fba_file)
        df1.as_matrix()

        df2 = pd.read_csv(compounds_file)
        df2.as_matrix()

        #----Input validation of Media/FBAs with Binary Matrix FBAs------
        # 1.0 Number of rows in Media.csv file =  (Number of columns -1)
        #   1.0. If they are different: Through an ERROR saying missed match number of FBAs in media and binary matrix. 
        # 1.1 Check whether the elements in Media.csv file contains only binary values (i.e. 0 and 1)
        #   1.1. If the elements are different: Through an ERROR saying not approapriate input values
        # 1.2 Check whether the compounds in Media.csv file match with number of FBAs
        #   1.2. If the compounds are different from number of FBAs: Through an ERROR saying not appropriate input values

        s_df1 = df1.shape
        s_df2 = df2.shape


        Temp_df2 = np.array(df2.values)
        # Create matrix with only the elements remove first column and all the rows
        Temp_df2 = Temp_df2[0:,1:]

        Bin_val_check =  np.array_equal(Temp_df2, Temp_df2.astype(bool))
        num_compounds = (s_df2[1])-1

        if ((s_df1[1]-1) != s_df2[0]) or (Bin_val_check != True) or (int(math.log(s_df2[0],2)) != num_compounds):
            print ('invalid input values')

        #-----All possible combination of the chemical compounds----------------------
        # 2.0 Sperating m0 from rest of the lables

        Temp1_df2 = df2

        cols = Temp1_df2.columns
        for i in range(1,len(cols)):
            Temp1_df2.loc[Temp1_df2[cols[i]] == 1 , cols[i]] = cols[i]

        print Temp1_df2

        # 2.1 Creating a disctionary for all FBAs except m0
        print len(Temp1_df2)
        mydict = {}
        for x in range(0, len(Temp1_df2)):
            for i in range(1,s_df2[1]):
                currentvalue = Temp1_df2.iloc[x,i]
                currentid = Temp1_df2.iloc[x,0]
                currentvalue = Temp1_df2.iloc[x,i]
                mydict.setdefault(currentid,[])
                if currentvalue > 0:
                    mydict[currentid].append(currentvalue)
                    
        # Add the first key as m0
        media_0_name = 'm0'
        mydict[media_0_name] = "['0']"
        #Sort the keys 
        mydict = collections.OrderedDict(natsort.natsorted(mydict.items())) 
        print mydict

        for k,v in mydict.iteritems():
            print k,v

        # List of Compounds combination in the list
        my_combi_list = []
        Compounds_Combi = list(range(1,num_compounds+1))
        for L in range(0, len(Compounds_Combi)+1):
            for subset in itertools.combinations(Compounds_Combi, L):
                my_combi_list.append(list(subset))
        print my_combi_list


        # Created a dictionary where the keys: 
        # list of compounds combination 
        # values are corresponding FBAs list in df2
        result_dict = {}
        for element in my_combi_list[1:]:
            for k, v in mydict.iteritems():
                if set(v).issubset(set(map(lambda x:str(x), element))):
                    key = ','.join(map(lambda x:str(x), element))
                    if result_dict.get(key):
                        media_list = result_dict[key]
                        media_list.append(k)
                        media_list = list(set(media_list))
                        result_dict.update({key: media_list})
                    else:
                        result_dict.update({key: [media_0_name, k]})            
        print result_dict

        # Created a dictionary where the keys are: 
        # list of compounds combination 
        # values are compounds combination FBAs with df1 vaules 
        All_Comp_Combi_dic = {}
        for column, value in result_dict.items():
            All_Comp_Combi_dic.update({column : df1.get(value)})   


        #To print an item from the All_Comp_Combi_dic
        df = (pd.DataFrame(All_Comp_Combi_dic.items()))

        #print df[0]
        #print df[1][7]

        MI_dict = {}
        for k in range(0, len(df[0])):          
            drop_rows_df = df[1][k].drop_duplicates(keep="first")
            drop_columns_df = drop_rows_df.T.drop_duplicates(keep="first").T
            remove = []
            removed = {}
            cols = df[1][k].columns
            for i in range(len(cols)-1):
                duplicated = []
                v = df[1][k][cols[i]].values
                for j in range(i+1,len(cols)):
                    if np.array_equal(v,df[1][k][cols[j]].values):
                        remove.append(cols[j])
                        duplicated.append(cols[j])
                if duplicated and cols[i] not in remove:
                    removed.update({cols[i]:duplicated})
                count = {}
                for key, value in removed.items():
                    count.update({key: len(value)})

                #print v

                # print drop_columns_df
                values = count.values()
                # print values
                values = map(lambda x: x+1, values)
                # print values
                d = {x:values.count(x) for x in values}     

                #-------Mutual Inforamtion (MI) calculation-------------
                FBAs = len(df[1][k].columns)
                pure_entropy = math.log(FBAs,2)
                #print pure_entropy


                # If No duplicates exist and list "value" is empty
                if not values:
                    #print("List is empty")
                    No_duplicate_FBAs = len(drop_columns_df.columns)
                    conditional_entropy = -1 * (No_duplicate_FBAs*((1/No_duplicate_FBAs)*((1/1)*math.log(1.0/1.0,2))));
                    Mutual_Info = pure_entropy - conditional_entropy
                    #print('Mutaul Info:', Mutual_Info)

                if values:
                # If duplicates exist and list "value" is not empty
                    conditional_entropy = 0
                    for key in d:
                        #print key, d[key]
                        Temp = -1 * d[key] * (key/float(FBAs)) * key * (1.0/key) * math.log(1.0/key,2)
                        conditional_entropy = Temp + conditional_entropy
                    #print "%3f" %Temp
                    Mutual_Info = pure_entropy - conditional_entropy
                    #print('Mutaul Info:', Mutual_Info)

                MI_dict[df[0][k]] = Mutual_Info

        #Sorted MI_dict
        MI_dict = sorted(MI_dict.items(), key=lambda x: (len(x[0]),x[1]), reverse=True)
        MI_dict = OrderedDict(MI_dict)
        
        print("Plot MI_dict")
        plt.bar(range(len(MI_dict)), MI_dict.values(), align='center', alpha=0.5, width=0.7)
        plt.xticks(range(len(MI_dict)), MI_dict.keys(), rotation='vertical')
        plt.xlabel('Compund Combinations')
        plt.ylabel('Mutual Information (in Bits)')
        plt.title("Organism:XYZ")
        fig1 = plt.gcf()
        fig1.savefig(os.path.join(self.scratch, 'MI_plot.png'), dpi=100)  

        return MI_dict
   
