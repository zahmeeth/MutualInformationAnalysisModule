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
    GIT_COMMIT_HASH = "90a8e7e2d30b52d9b4c6b8c079770be791c41eb3"

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
           compound id.), parameter "workspace" of String, parameter
           "media_id" of String
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
        compounds = params.get('compounds')
        media_id = params.get('media_id')
        workspace_name = params.get('workspace_name')

        #compounds_file = MutualInfoUtil._get_file_from_ws(compounds)
        #fba_file = MutualInfoUtil._get_file_from_ws(fba_object_ref)
        fba_file = '/kb/module/data/BT_7bits.csv'
        compounds_file = '/kb/module/data/AllFBAs_7.csv'
        mutual_info = MI_runner._generate_mutual_info(compounds_file, fba_file)
        print('YAY!')
        output = MI_runner._generate_report(self.scratch, mutual_info, params)    
        
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
