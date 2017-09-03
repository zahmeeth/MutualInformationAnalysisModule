import time
import os
import errno
import uuid
import math
import pandas as pd
import numpy as np
import collections
import natsort
import uuid
import shutil
import itertools
import json
from itertools import combinations
import matplotlib
from fba_tools.fba_toolsClient import fba_tools
import matplotlib.pyplot as plt
from collections import OrderedDict
from copy import deepcopy

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

	def _get_file_from_ws(self, workspace, obj_name):
		try:
			file_path = self.ws.get_objects(
				[{'name': obj_name,
				  'workspace': workspace}])[0]
		except Exception as e:
			raise ValueError(
				'Unable to get object from workspace: (' +
				workspace + '/' + obj_name + ')' + str(e))
		return file_path

	def _make_media_files(self, ws_name, base, compounds):
		"""
		Build and store media objects for each combination of compound added to the base media.
		:param base: The base media file
		:param compounds: the set of compound to test
		:return: A list of media ids and a matrix with each media combination defined
		"""
		base_media = self._get_file_from_ws(ws_name, base)['data']

		media_ids = [base_media['id']]
		new_media_list = []
		media_matrix = [[""]+compounds]
		media_matrix.append([[base_media['id']]+[0]*len(compounds)])
		for n_comp in range(1, len(compounds)+1):
			for combo in combinations(compounds, n_comp):
				new_media_id = base_media['id'] + '_v%s' % len(media_matrix)
				media_ids.append(new_media_id)
				media_matrix.append([new_media_id]+[1 if comp in combo else 0 for comp in compounds])
				new_media = deepcopy(base_media)
				new_media['id'] = new_media_id
				new_media['name'] = new_media_id
				for new_comp in combo:
					new_media['mediacompounds'].append(
						{'compound_ref': '48/1/1/compounds/id/%s' % new_comp.split('_')[0],
						 'concentration': 1.0, 'maxFlux': 1000, 'minFlux': -1000})
				new_media_list.append(new_media)

		print("Made %s Media Files" % (len(media_ids)-1))
		info = self.ws.save_objects(
			{'workspace': ws_name,
			 "objects": [{
				 "type": "KBaseBiochem.Media",
				 "data": media,
				 "name": media['name']
			 } for media in new_media_list]
			 })
		print info
		return media_ids, media_matrix

	def _run_fba(self, workspace_name, media_id_list, fbamodel_id):
		fba_tool_obj = fba_tools(self.callback_url)
		fba_tool_obj.run_flux_balance_analysis({
			"workspace" : workspace_name,
			"fbamodel_id" : fbamodel_id,
			"fba_output_id" : fbamodel_id + ".mifba",
			"fbamodel_workspace" : workspace_name,
			"media_id_list" : media_list,
			"target_reaction" : "bio1",
			"minimize_flux" : 1
			})
		output = self.ws.get_objects2({
			'objects' : [{
				'ref' : workspace_name + "/" + fbamodel_id + '.mifba'
			}]
			})
		
		fba = output['data'][0]['data']
		biomass_data = "FBAs,Biomass\n"
		secretion_file = ","+','.join(media_list)+"\n"
		full_secretion_file = ","+','.join(media_list)+"\n"
		full_flux_file = ","+','.join(media_list)+"\n"
		flux_file = ","+','.join(media_list)+"\n"
		objectives = fba['other_objectives']
		for i in range(0, len(objectives)):
			biomass_data = biomass_data + media_list[i] + "," + objectives[i] + "\n"
		
		flux_vars = fba['FBAReactionVariables']
		for var in flux_vars:
			id = var['modelreaction_ref'].split("/").pop()
			flux_file = flux_file + id
			full_flux_file = full_flux_file + id
			fluxes = var['other_values']
			for i in range(0, len(fluxes)):
				full_flux_file = full_flux_file + "," + fluxes[i]
				if abs(fluxes[i]) < 1e-7:
					flux_file = flux_file + ",0"
				else:
					flux_file = flux_file + ",1"
			flux_file = flux_file + "\n"
			full_flux_file = full_flux_file + "\n"
		
		secretion_vars = fba['FBACompoundVariables']
		for var in secretion_vars:
			id = var['modelcompound_ref'].split("/").pop()
			secretion_file = secretion_file + id
			full_secretion_file = full_secretion_file + id
			fluxes = var['other_values']
			for i in range(0, len(fluxes)):
				full_secretion_file = full_secretion_file + "," + fluxes[i]
				if abs(fluxes[i]) < 1e-7:
					secretion_file = secretion_file + ",0"
				elif fluxes[i] < 0:
					secretion_file = secretion_file + ",-1"
				else:
					secretion_file = secretion_file + ",1"
			secretion_file = secretion_file + "\n"
			full_secretion_file = full_secretion_file + "\n"
			
		output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
		self._mkdir_p(output_directory)
		biomass_path = os.path.join(output_directory, 'biomass.csv')
		secretion_path = os.path.join(output_directory, 'secretion.csv')
		flux_path = os.path.join(output_directory, 'flux.csv')
		full_secretion_path = os.path.join(output_directory, 'full_secretion.csv')
		full_flux_path = os.path.join(output_directory, 'full_flux.csv')
		
		with open(biomass_path, 'w') as biomass_f:
			biomass_f.write(biomass_data)
		
		with open(secretion_path, 'w') as secretion_f:
			secretion_f.write(secretion_file)
		
		with open(flux_path, 'w') as flux_f:
			flux_f.write(flux_file)
			
		with open(full_secretion_path, 'w') as full_secretion_f:
			full_secretion_f.write(full_secretion_file)
		
		with open(full_flux_path, 'w') as full_flux_f:
			full_flux_f.write(full_flux_file) 
		
		return [biomass_path,secretion_path,flux_path,full_secretion_path,full_flux_path]
	
	def _generate_html_report(self, result_directory, mutual_info_dict):
		
		"""
		_generate_html_report: generate html summary report
		"""
		#scratch, uui, datafileutil, file_to_shock, shockId, extended report

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

		report_shock_id = self.dfu.file_to_shock({'file_path': output_directory,
												  'pack': 'zip'})['shock_id']

		html_report.append({'shock_id': report_shock_id,
							'name': os.path.basename(result_file_path),
							'label': os.path.basename(result_file_path),
							'description': 'HTML summary report for Mutual Information App'})
		
		return html_report

	def _generate_report(self, result_directory, mutual_info_dict, params,paths):
		"""
		_generate_report: generate summary report
		"""

		uuidStr = str(uuid.uuid4())
		self._mkdir_p(result_directory + '/' + uuidStr)

		shutil.copy('/kb/module/data/index.html', result_directory + '/' + uuidStr + '/index.html')
		shutil.copy('pdata.json', result_directory + '/' + uuidStr + '/pdata.json')

		# DataFileUtils to shock
		report_shock_id = self.dfu.file_to_shock({'file_path': result_directory + '/' + uuidStr,
												  'make_handler': 0,
												  'pack': 'zip'})['shock_id']

		report_file = {'name': 'index.html',
					   'description': 'the report',
					   'shock_id': report_shock_id}
					   
		biomass_file = {'name': 'biomass_file.csv',
					   'description': 'biomass_file',
					   'path': paths[0]}
		
		flux_file = {'name': 'flux_file.csv',
					   'description': 'flux_file',
					   'path': paths[1]}
					   
		full_flux_file = {'name': 'full_flux_file.csv',
					   'description': 'full_flux_file',
					   'path': paths[2]}
					   
		secretion_file = {'name': 'secretion_file.csv',
					   'description': 'secretion_file',
					   'path': paths[3]}
					   
		full_secretion_file = {'name': 'full_secretion_file.csv',
					   'description': 'full_secretion_file',
					   'path': paths[4]}

		log('creating report')
		#output_html_files = self._generate_html_report(result_directory,
		#											   mutual_info_dict)
		report_params = {'message': '',
						 'workspace_name': params.get('workspace_name'),
						 'html_links': [report_file],
						 'file_links': [biomass_file,flux_file,full_flux_file,secretion_file,full_secretion_file],
						 'direct_html_link_index': 0,
						 'html_window_height': 333,
						 'report_object_name': 'MutualInfomation_report_' + uuidStr}

		kbase_report_client = KBaseReport(self.callback_url)
		output = kbase_report_client.create_extended_report(report_params)

		report_output = {'report_name': output['name'], 'report_ref': output['ref']}

		return report_output	   

	def _generate_mutual_info(self, media_matrix, fba_file, mi_options):

		df1 = pd.read_csv(fba_file)
		df1.as_matrix()


	   #----Input validation of Media/FBAs with Binary Matrix FBAs------
		# 1.0 Number of rows in Media.csv file =  (Number of columns -1)
		#   1.0. If they are different: Through an ERROR saying missed match number of FBAs in media and binary matrix. 
		# 1.1 Check whether the elements in Media.csv file contains only binary values (i.e. 0 and 1)
		#   1.1. If the elements are different: Through an ERROR saying not approapriate input values
		# 1.2 Check whether the compounds in Media.csv file match with number of FBAs
		#   1.2. If the compounds are different from number of FBAs: Through an ERROR saying not appropriate input values

		print media_matrix

		s_df1 = df1.shape
		s_df2 = media_matrix.shape


		Temp_df2 = np.array(media_matrix.values)
		# Create matrix with only the elements remove first column and all the rows
		Temp_df2 = Temp_df2[0:,1:]

		Bin_val_check =  np.array_equal(Temp_df2, Temp_df2.astype(bool))
		num_compounds = (s_df2[1])-1

		if ((s_df1[1]-1) != s_df2[0]) or (Bin_val_check != True) or (int(math.log(s_df2[0],2)) != num_compounds):
			print ('invalid input values')

		#-----All possible combination of the chemical compounds----------------------
		# 2.0 Sperating m0 from rest of the lables

		Temp1_df2 = media_matrix

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

				MI_dict[df[0][k]] = Mutual_Info

		#Sorted MI_dict
		MI_dict = sorted(MI_dict.items(), key=lambda x: (-len(x[0]), x[0]))
		MI_dict = OrderedDict(MI_dict)

		x_groups = [[] for x in range(num_compounds)]
		y_groups = [[] for x in range(num_compounds)]
		names = [[] for x in range(num_compounds)]
		Comp_Mapping = [[] for x in range(num_compounds)]

		for key, val in MI_dict.iteritems():
			del_count = key.count(',')
			x_groups[del_count].append(key)
			y_groups[del_count].append(val)

			# for x, y in zip(x_groups, y_groups):
			# data.append(go.Bar(x=x, y=y, name='test'))

		compound_IDs = ['H2', 'Vitamin K', 'Hematine', 'Glucose', 'Acetate', 'Formate', 'B12']

		pdata = []
		for i in range(0, len(x_groups)):
			names[i] = str(i + 1) + ' Compound Combination'
			Comp_Mapping = str(i + 1) + '-' + compound_IDs[i]

			record = {}
			record["x"] = []
			for e in x_groups[i]:
				record["x"].append("c" + e)
			record["y"] = y_groups[i]
			record["names"] = names[i]
			record["Comp_Mapping"] = Comp_Mapping
			pdata.append(record)

		print pdata

		with open('pdata.json', 'w') as outfile:
			json.dump(pdata, outfile)
		return MI_dict
   
