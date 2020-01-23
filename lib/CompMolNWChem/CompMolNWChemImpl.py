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
from installed_clients.CompoundSetUtils import CompoundSetUtils
import compound_parsing as parse
import export as ex
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
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    def _check_param(in_params, req_param, opt_param=list()):
        """
        Check if each of the params in the list are in the input params
        """
        for param in req_param:
            if param not in in_params:
                raise ValueError('{} parameter is required'.format(param))
        defined_param = set(req_param+opt_param)
        for param in in_params:
            if param not in defined_param:
                logging.warning("Received unexpected parameter {}".format(param))

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

    def _export_compound_set(self, ref, file_type):
        logging.info("Exporting {} as {}".format(ref, file_type))
        compoundset = self.dfu.get_objects(
            {'object_refs': [ref]}
        )['data'][0]['data']
        temp_dir = "{}/{}".format(self.scratch, uuid.uuid4())
        os.mkdir(temp_dir)
        out_dir = "{}/{}".format(temp_dir, compoundset['name'])
        os.mkdir(out_dir)
        target = "{}/{}.{}".format(out_dir, compoundset['name'], file_type)
        if file_type == 'tsv':
            parse.write_tsv(compoundset, target)
        elif file_type == 'sdf':
            parse.write_sdf(compoundset, target)
        else:
            raise ValueError("Bad file_type: {}".format(file_type))
        handle = self.dfu.package_for_download(
            {'file_path': out_dir, 'ws_refs': [ref]})
        output = {'shock_id': handle['shock_id']}
        return output

    def _fetch_mol2_files(self, ref):
        compoundset_obj = self.dfu.get_objects(
            {'object_refs': [ref]}
        )['data'][0]
        compoundset_info = compoundset_obj['info']
        compoundset = compoundset_obj['data']
        temp_dir = "{}/{}".format(self.scratch, uuid.uuid4())
        os.mkdir(temp_dir)

        compounds = compoundset.get('compounds')

        mol2_files = []
        comp_id_mol2_file_name_map = {}
        for compound in compounds:
            mol2_handle_ref = compound.get('mol2_handle_ref')

            if mol2_handle_ref:
                mol2_file_path = self.dfu.shock_to_file(
                                            {'handle_id': mol2_handle_ref,
                                             'file_path': temp_dir}).get('file_path')
                mol2_files.append(mol2_file_path)
                comp_id_mol2_file_name_map[compound['id']] = os.path.basename(mol2_file_path)

        packed_mol2_files_path = None
        if mol2_files:
            packed_mol2_files_path = os.path.join(temp_dir, compoundset_info[1] + '_mol2_files.zip')
            with zipfile.ZipFile(packed_mol2_files_path, 'w') as zipMe:
                for file in mol2_files:
                    zipMe.write(file, compress_type=zipfile.ZIP_DEFLATED)

        return packed_mol2_files_path, comp_id_mol2_file_name_map

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

        sim_dir = '~/../simulation'
        os.system('ls')
        import pandas as pd
        
        # Read inputs from .tsv file

        df = CompoundSetUtils.compound_set_from_file(self,ctx,params['Input_File'])
        
        #df = pd.read_csv(params['Input_File'], sep ='\t')
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

        #for j in range(length):
        #    os.chdir('./'+ids[j]+'/dft')
        #    os.system('ls')

        #    with open(ids[j]+'_Mulliken.mol2') as IN:
                
            
        #        AllChem.MolToSmiles(IN,True)
            #xx = com._make_compound_info(open(ids[j]+'_Mulliken.mol2'))
            #print(xx)
                
            #os.chdir('../..')

        #
        #xx = ex._save_to_ws_and_report(self,params['workspace_name'],'~/../../simulation','compoundset')
        
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
