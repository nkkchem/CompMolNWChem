# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid
import zipfile
import copy
from subprocess import Popen, PIPE

#import CompoundSetUtils.compound_parsing as parse
#import CompoundSetUtils.zinc_db_util as zinc_db_util
#from installed_clients.DataFileUtilClient import DataFileUtil
#from installed_clients.KBaseReportClient import KBaseReport

import '/Users/mcna892/bin/CompMolNWChem/lib/installed_clients/DataFileUtilClient.py' 

import compound_parsing as parse
#END_HEADER

class CompoundSetUtils:
    
    def compound_set_to_file(self, ctx, params):
        """
        CompoundSetToFile
        string compound_set_name
        string output_format
        :param params: instance of type "compoundset_download_params" ->
           structure: parameter "compound_set_ref" of String, parameter
           "output_format" of String
        :returns: instance of type "compoundset_download_results" ->
           structure: parameter "file_path" of String, parameter
           "packed_mol2_files_path" of String, parameter
           "comp_id_mol2_file_name_map" of mapping from String to String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN compound_set_to_file
        self._check_param(params, ['compound_set_ref', 'output_format'])
        ret = self.dfu.get_objects(
            {'object_refs': [params['compound_set_ref']]}
        )['data'][0]
        compoundset = ret['data']
        ext = params['output_format']
        out = f"{self.scratch}/{uuid.uuid4()}"
        os.mkdir(out)
        out += f"/{compoundset['name']}"
        if ext == 'sdf':
            outfile_path = parse.write_sdf(compoundset, out)
        elif ext == 'tsv':
            outfile_path = parse.write_tsv(compoundset, out)
        else:
            outfile_path = parse.write_mol_dir(compoundset, out, ext)

        packed_mol2_files_path, comp_id_mol2_file_name_map = self._fetch_mol2_files(
                                                                    params['compound_set_ref'])

        output = {'file_path': outfile_path, 'packed_mol2_files_path': packed_mol2_files_path,
                  'comp_id_mol2_file_name_map': comp_id_mol2_file_name_map}

        #END compound_set_to_file

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method compound_set_to_file return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

CompoundSetUtils.compound_set_to_file(self,ctx,'/Users/mcna892/bin/CompMolNWChem/test/test_compounds.tsv')
