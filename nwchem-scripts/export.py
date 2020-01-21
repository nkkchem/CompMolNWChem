# -*- coding: utf-8 -*-
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
import compound_parsing as parse

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
