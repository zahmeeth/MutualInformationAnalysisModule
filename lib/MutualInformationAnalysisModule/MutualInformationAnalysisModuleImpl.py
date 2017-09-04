# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import json

from MutualInformationAnalysisModule.Utils.MutualInfoUtil import MutualInfoUtil


#END_HEADER


class MutualInformationAnalysisModule:
    '''
    Module Name:
    MutualInformationAnalysisModule

    Module Description:
    A KBase module: MutualInformationAnalysisModule
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.0.0"
    GIT_URL = "https://github.com/zahmeeth/MutualInformationAnalysisModule.git"
    GIT_COMMIT_HASH = "5a6e49e71ea48c5286b2cd52749c412a3a7669b3"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        self.scratch = config['scratch']
        #END_CONSTRUCTOR
        pass

    def run_flux_mutual_information_analysis(self, ctx, params):
        """
         :param params: instance of type
            "RunFluxMutualInformationAnalysisParams" -> structure: parameter
            "fbamodel_id" of type "ws_fbamodel_id" (The workspace ID for a
            FBAModel data object. @id ws KBaseFBA.FBAModel), parameter
            "compounds" of list of type "compound_id" (A string representing a
            compound id.), parameter "workspace_name" of String, parameter
            "media_id" of String, parameter "mi_options" of String
         :returns: instance of type "RunFluxMutualInformationAnalysisResults"
            -> structure: parameter "report_name" of String, parameter
            "report_ref" of type "ws_report_id" (The workspace ID for a Report
            object @id ws KBaseReport.Report)
         """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_flux_mutual_information_analysis
        print('Starting flux mutual information analysis method.')

        MI_runner = MutualInfoUtil(self.config)
        MI_runner._validate_run_flux_mutual_information_analysis_params(params)
        fbamodel_id = params.get('fbamodel_id')
        # compounds is maybe a string delimited by ','
        compounds = params.get('compounds').split(',')
        #if isinstance(params.get('compounds'), str):
        #    compounds = params.get('compounds').split(',')
        media_id = params.get('media_id')
        workspace_name = params.get('workspace_name')

        print(params)

        print('Making Media Objects')
        media_id_list, media_matrix, myuuid = MI_runner._make_media_files(workspace_name, media_id, compounds)
        [biomass_path,secretion_path,flux_path,full_secretion_path,full_flux_path] = MI_runner._run_fba(workspace_name, media_id_list, fbamodel_id, myuuid, media_id)
        # Loading media matrix - which is shared by all three modes of the function
        import pandas as pd
        media_matrix = pd.read_csv('/kb/module/data/AllFBAs_7.csv')
        media_matrix.as_matrix()
        # Loading fluxes when running in flux mode
        data_file = "";
        if params['mi_options'] == "flux":
            data_file = flux_path
        elif params['mi_options'] == "biomass":
            data_file = biomass_path
        elif params['mi_options'] == "secretion":
            data_file = secretion_path
        # Running core mutual information function
        #Temporary measure for debugging
        mutual_info = {}
        for media in media_id_list:
        	mutual_info[media] = 1
        #mutual_info = MI_runner._generate_mutual_info(media_matrix, data_file, params['mi_options'])
        # Writing output report
        output = MI_runner._generate_report(self.scratch, mutual_info, params,[biomass_path,secretion_path,flux_path,full_secretion_path,full_flux_path])

        #END run_flux_mutual_information_analysis

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_flux_mutual_information_analysis return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]