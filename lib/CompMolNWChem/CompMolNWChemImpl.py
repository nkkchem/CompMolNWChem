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
import re
import zipfile
import uuid
import copy
import shutil
import argparse

os.system('pip install snakemake')

from pybel import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from csv import DictReader
#from mulitprocessing import cpu_count
from os.path import *
from pkg_resources import resource_filename
from cme import *
from snakemake import snakemake
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
    GIT_COMMIT_HASH = "496b69853840cfaacbd146fce50466e73701bcc7"

    #BEGIN_CLASS_HEADER
    def _generate_output_file_list(self, result_directory):
        """
        _generate_output_file_list: zip result files and generate file_links for report
        """

        #log('start packing result files')
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
                        zip_file.write(os.path.join(root, file), 
                                       os.path.join(os.path.basename(root), file))

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
                        zip_file.write(os.path.join(root, file), 
                                       os.path.join(os.path.basename(root), file))

        output_files.append({'path': plot_file,
                             'name': os.path.basename(plot_file),
                             'label': os.path.basename(plot_file),
                             'description': 'Visualization plots by DESeq2 App'})

        return output_files

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

    def _save_to_ws_and_report(self, ws_id, source, compoundset, message=None):
        """Save compound set to the workspace and make report"""
        info = self.dfu.save_objects(
            {'id': ws_id,
             "objects": [{
                 "type": "KBaseBiochem.CompoundSet",
                 "data": compoundset,
                 "name": compoundset['name']
             }]})[0]
        compoundset_ref = "%s/%s/%s" % (info[6], info[0], info[4])
        if not message:
            message = 'Imported %s as %s' % (source, info[1])
        report_params = {
            'objects_created': [{'ref': compoundset_ref,
                                 'description': 'Compound Set'}],
            'message': message,
            'workspace_name': info[7],
            'report_object_name': 'compound_set_creation_report'
        }

        # Construct the output to send back
        report_client = KBaseReport(self.callback_url)
        report_info = report_client.create_extended_report(report_params)
        output = {'report_name': report_info['name'],
                  'report_ref': report_info['ref'],
                  'compoundset_ref': compoundset_ref}
        return output

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

        
        # Load the tsv file into a compound set using DataFileUtil methods
        
        scratch_file_path = self.dfu.download_staging_file({'staging_file_subdir_path':params['Input_File']}
                                       ).get('copy_file_path')

        #print('Scratch File Path: ',scratch_file_path)

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
        #elif ext == '.csv':
        #    compounds = parse.read_csv(scratch_file_path,
        #                               mol2_file_dir=mol2_file_dir,
        #                               callback_url=self.callback_url)
        #else:
        #    raise ValueError('Invalid input file type. Expects .tsv or .sdf')

        #DEBUG::
        #print('Compounds:',compounds)

#        compoundset = {
#            'id': params['Input_File'],
#            'name': params['Input_File'],
#            'description': 'Compound Set produced from %s' % file_name,
#            'compounds': compounds,
#        }

        # Finish Reading in Compound Set

        # Read ids and smiles from compound set for nwchem input
        
#        ids = []
#        smiles = []

#        for d in compounds:
#           ids.append(d['id'])
#           smiles.append(d['smiles'])
        #print(ids)
        #print(smiles)
        

        
        # Read the ids and structures of the compounds
        
#        its.inchi_to_dft(ids,smiles)

        #DEBUG::
        #os.system('pwd')
        #os.system('ls')
        
#        length = len(ids)
#        for i in range(length):
#            os.chdir('./'+ids[i]+'/dft')
#            x = ids[i] + '_nwchem.out'
            #print('x:',x)
#            file1 = open(x, 'r')
#            nAtoms = mul.getNumberOfAtoms(file1)
#            energy = mul.getInternalEnergy0K(file1)
#            charge =mul.getMullikenCharge(file1,nAtoms)
#            file1.close()
           
#            mul.nAtoms = nAtoms
#            mul.E0K = energy

#            mul.calculate(ids[i])

       
        from snakemake import snakemake

        inchilist = scratch_file_path
        
        with open(inchilist,'r') as f:
#            if 'inchi.csv' in inchilist:
            reader = csv.reader(f,delimiter=(','))

            lineCount = 0
            for row in reader:
                if lineCount == 0:
                    InChIKey = row.index("InChI-Key")
                    InChI = row.index("InChI")
                    lineCount+=1
                else:
                    print(row)
                    inchikey_str=row[InChIKey].split('=')[1]
                    print('inchikey string ', inchikey_str)

                    moldir = inchikey_str
                    if not os.path.exists(moldir):
                        os.mkdir(moldir)
    
                    initial_structure_dir = moldir + '/initial_structure'
                    if not os.path.exists(initial_structure_dir):
                        os.mkdir(initial_structure_dir)

                    md_structure_dir = moldir + '/md'
                    if not os.path.exists(md_structure_dir):
                        os.mkdir(md_structure_dir)

                    dft_structure_dir = moldir + '/dft'
                    if not os.path.exists(dft_structure_dir):
                        os.mkdir(dft_structure_dir)

                    inchifile_str = initial_structure_dir + '/' + inchikey_str + '.inchi'
                    with open(inchifile_str,'w+') as f:
                        f.write(row[InChI])

        os.system('snakemake -p --cores 2 --snakefile snakemake/final_pipeline.snakemake -w 300')

        # Build KBase Output. Should output entire /simulation directory and build a CompoundSet with Mol2 Files

        result_directory = '/simulation/'

        ## Build CompoundSet with Mol2 Files... similarly to fetch_mol2_files_from_zinc (CompoundSetUtils)....

#        compoundset_copy = copy.deepcopy(compoundset)

#        count = 0

#        for compound in compoundset_copy.get('compounds'):
#            if not compound.get('mol2_handle_ref'):
#                mol2_file_path = result_directory+compound.get('id')
#                SMILES = compound.get('smiles')

#                shutil.move(mol2_file_path,self.scratch)

#                os.chdir(self.scratch)
               
#                mol2_file_path = self.scratch + '/'+ compound.get('id')+'/dft/' + compound.get('id')+'_Mulliken.mol2'              
#                handle_id = self.dfu.file_to_shock({'file_path': mol2_file_path,
#                                                    'make_handle': True})['handle']['hid']
#                print('Handle ID:',handle_id)
#                compound['mol2_handle_ref'] = handle_id
#                count += 1

               
               
#        if count:
#            message = 'Successfully fetched {} Mol2 files from Staging Path'.format(count)


        ## Create Extended Report
        
        output_files = self._generate_output_file_list(self.scratch)


        report_params = {'workspace_id': params['workspace_id'],
                         'objects_created': [],
                         'report_object_name': 'kb_deseq2_report_' + str(uuid.uuid4())}

        report = KBaseReport(self.callback_url)
        
        report_info = report.create_extended_report(report_params)

        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }

        return [output]

        #END run_CompMolNWChem
        
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
