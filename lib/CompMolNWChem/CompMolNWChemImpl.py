# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
import subprocess as _subprocess
import csv
import inchi_to_submission as its
import extract_properties_mulliken_charges_mol2 as mul
import compound_parsing as com
import pandas as pd
import compound_parsing as parse
import export as ex

from pybel import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from csv import DictReader
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.CompoundSetUtilsClient import CompoundSetUtils

#from CompMolNWChem import CompoundSetUtil

#END_HEADER


class CompMolNWChem:
    '''
    Module Name:
    CompMolNWChem

    Module Description:
    A KBase module: CompMolNWChem
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/nkkchem/CompMolNWChem.git"
    GIT_COMMIT_HASH = "0e157ea3f395c544a04c1542473be59ec39129ef"

    #BEGIN_CLASS_HEADER
    # staging file prefix
    #STAGING_GLOBAL_FILE_PREFIX = '/data/bulk/'
    #STAGING_USER_FILE_PREFIX = '/staging/'
    
   # def _get_staging_file_path(self, token_user, staging_file_subdir_path):
   #     """
   #     _get_staging_file_path: return staging area file path
   #     directory pattern:
   #         perfered to return user specific path: /staging/sub_dir/file_name
   #         if this path is not visible to user, use global bulk path: /data/bulk/user_name/sub_dir/file_name
   #     """
#
#        user_path = os.path.join(self.STAGING_USER_FILE_PREFIX, staging_file_subdir_path.strip('/'))
#
#        if os.path.exists(user_path):
#            return user_path
#        else:
#            return os.path.join(self.STAGING_GLOBAL_FILE_PREFIX, token_user,
#                                staging_file_subdir_path.strip('/'))
    
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        self.comp = CompoundSetUtils(self.callback_url)
        self.scratch = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #self.scratch = config['scratch']

        #END_CONSTRUCTOR
        pass


    def run_CompMolNWChem(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_CompMolNWChem

        # Initial Tests to Check for Proper Inputs
        for name in ['Input_File','calculation_type','workspace_name']:
            if name not in params:
                raise ValueError('Parameter "' + name + '"is required but missing')
        if not isinstance(params['Input_File'], str):
            raise ValueError('Input_File must be a string')

        
        ### Test Space
        
        scratch_file_path = self.dfu.download_staging_file({'staging_file_subdir_path':params['Input_File']}
                                       ).get('copy_file_path')

        print('Scratch File Path: ',scratch_file_path)

        mol2_file_dir = None        
        ext = os.path.splitext(scratch_file_path)[1]
        file_name = os.path.basename(scratch_file_path)
        if ext == '.sdf':
            compounds = parse.read_sdf(scratch_file_path,
                                       mol2_file_dir=mol2_file_dir,
                                       callback_url=self.callback_url)
        elif ext == '.tsv':
            compounds = parse.read_tsv(scratch_file_path,
                                       mol2_file_dir=mol2_file_dir,
                                       callback_url=self.callback_url)
        else:
            raise ValueError('Invalid input file type. Expects .tsv or .sdf')


        #print('Compounds:',compounds)

        compoundset = {
            'id': params['Input_File'],
            'name': params['Input_File'],
            'description': 'Compound Set produced from %s' % file_name,
            'compounds': compounds,
        }

        print(compounds)
        ###
        
        # Read ids and smiles from compound set for nwchem input
        
        ids = []
        smiles = []

        for d in compounds:
           ids.append(d['id'])
           smiles.append(d['smiles'])
        print(ids)
        print(smiles)

        # TODO: Initialize CompoundSetUtils and pass inputs 
        
        # Read the ids and structures of the compounds
                       
        its.inchi_to_dft(ids,smiles)
        
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
        
        #END run_CompMolNWChem

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_CompMolNWChem return value ' +
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
