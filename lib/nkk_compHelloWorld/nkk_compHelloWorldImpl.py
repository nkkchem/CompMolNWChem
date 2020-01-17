# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
import subprocess as _subprocess
import csv
from pybel import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from csv import DictReader
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil 

#END_HEADER


class nkk_compHelloWorld:
    '''
    Module Name:
    nkk_compHelloWorld

    Module Description:
    A KBase module: nkk_compHelloWorld
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_nkk_compHelloWorld(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_nkk_compHelloWorld

        sim_dir = '~/../simulation'
        os.system('ls')
        import pandas as pd
        
        # Read inputs from .tsv file

        df = pd.read_csv(params['Input_File'], sep ='\t')
        ids = df['id']
        InChIes = df['structure']

        import inchi_to_submission as its
        import extract_properties_mulliken_charges_mol2 as mul
        import compound_parsing as com
        
        its.inchi_to_dft(ids,InChIes)

        length = len(ids)
        for i in range(length):
            os.chdir('./'+ids[i]+'/dft')
            file1 = open('nwchem.out', 'r')
            nAtoms = mul.getNumberOfAtoms(file1)
            energy = mul.getInternalEnergy0K(file1)
            charge =mul.getMullikenCharge(file1,nAtoms)
            file1.close()
           
            mul.nAtoms = nAtoms
            mul.E0K = energy

            mul.calculate(ids[i])

        for j in range(length):
            os.chdir('./'+ids[j]+'/dft')
            os.system('ls')

            with open(ids[j]+'_Mulliken.mol2') as IN:
                
            
                AllChem.MolToSmiles(IN,True)
            #xx = com._make_compound_info(open(ids[j]+'_Mulliken.mol2'))
            #print(xx)
                
            #os.chdir('../..')
            
        report = KBaseReport(self.callback_url)
        report_info = report.create({'report':{'objects_created':[],
                                               'text_message': params['Input_File'],
                                               'text_message':params['calculation_type']},
                                               'workspace_name': params['workspace_name']})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }

        
        return [output]
        
        #END run_nkk_compHelloWorld

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_nkk_compHelloWorld return value ' +
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
