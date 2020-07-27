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
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
import compound_parsing as parse
#END_HEADER


class CompoundSetUtils:
    '''
    Module Name:
    CompoundSetUtils
    Module Description:
    A KBase module: CompoundSetUtils
    Contains tools for import & export of compound sets
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "2.1.2"
    GIT_URL = "https://github.com/Tianhao-Gu/CompoundSetUtils.git"
    GIT_COMMIT_HASH = "12e1f23022354f475d7ceb3631913956eb5831a7"

    #BEGIN_CLASS_HEADER
    @staticmethod
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

    def _covert_mol2_files_to_pdbqt(self, ref):
        compoundset_obj = self.dfu.get_objects(
            {'object_refs': [ref]}
        )['data'][0]
        compoundset_info = compoundset_obj['info']
        compoundset = compoundset_obj['data']
        mol2_temp_dir = "{}/{}".format(self.scratch, uuid.uuid4())
        os.mkdir(mol2_temp_dir)
        pdbqt_temp_dir = "{}/{}".format(self.scratch, uuid.uuid4())
        os.mkdir(pdbqt_temp_dir)

        compounds = compoundset.get('compounds')

        pdbqt_files = []
        comp_id_pdbqt_file_name_map = {}
        for compound in compounds:
            mol2_handle_ref = compound.get('mol2_handle_ref')

            if mol2_handle_ref:
                mol2_file_path = self.dfu.shock_to_file(
                                            {'handle_id': mol2_handle_ref,
                                             'file_path': mol2_temp_dir}).get('file_path')
                pdbqt_file_path = os.path.join(pdbqt_temp_dir, compound['id'] + '.pdbqt')

                command = ['obabel', '-i', 'mol2', mol2_file_path, '-o', 'pdbqt', '-O', pdbqt_file_path]
                process = Popen(command, stdout=PIPE, stderr=PIPE)
                stdout, stderr = process.communicate()

                if 'converted' in str(stderr) and 'molecule' in str(stderr):
                    logging.info('Successfully converted Mol2 to pdbqt format: {}'.format(
                                                            os.path.basename(mol2_file_path)))
                    pdbqt_files.append(pdbqt_file_path)
                    comp_id_pdbqt_file_name_map[compound['id']] = os.path.basename(
                                                                            pdbqt_file_path)
                else:
                    logging.warning('Cannot convert Mol2 file to pdbqt format: {}'.format(
                                                            os.path.basename(mol2_file_path)))
                    logging.warning(stderr)

        packed_pdbqt_files_path = None
        if pdbqt_files:
            packed_pdbqt_files_path = os.path.join(pdbqt_temp_dir, compoundset_info[1] + '_pdbqt_files.zip')
            with zipfile.ZipFile(packed_pdbqt_files_path, 'w') as zipMe:
                for file in pdbqt_files:
                    zipMe.write(file, compress_type=zipfile.ZIP_DEFLATED)

        return packed_pdbqt_files_path, comp_id_pdbqt_file_name_map
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.scratch = config['scratch']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        #END_CONSTRUCTOR
        pass


    def compound_set_from_file(self, ctx, params):
        """
        CompoundSetFromFile
        string staging_file_path
        :param params: instance of type "compoundset_upload_params" ->
           structure: parameter "workspace_id" of String, parameter
           "staging_file_path" of String, parameter "compound_set_name" of
           String, parameter "mol2_staging_file_path" of String
        :returns: instance of type "compoundset_upload_results" -> structure:
           parameter "report_name" of String, parameter "report_ref" of
           String, parameter "compoundset_ref" of type "obj_ref"
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN compound_set_from_file
        self._check_param(params, ['workspace_id', 'staging_file_path',
                                   'compound_set_name'], opt_param=['mol2_staging_file_path'])
        scratch_file_path = self.dfu.download_staging_file(
            {'staging_file_subdir_path': params['staging_file_path']}
        ).get('copy_file_path')
        # I probably should be uploading the raw files to shock

        mol2_staging_file_path = params.get('mol2_staging_file_path')

        mol2_file_dir = None
        if mol2_staging_file_path:
            mol2_scratch_file_path = self.dfu.download_staging_file(
                    {'staging_file_subdir_path': mol2_staging_file_path}
                ).get('copy_file_path')

            try:
                logging.info("start unpacking mol2 file")
                mol2_file_path_out = self.dfu.unpack_file(
                                            {'file_path': mol2_scratch_file_path})['file_path']
                mol2_file_dir = os.path.dirname(mol2_file_path_out)
            except Exception:
                raise ValueError('Cannot unpack mol2 file: {}'.format(
                                                        os.path.basename(mol2_file_path_out)))

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

        compoundset = {
            'id': params['compound_set_name'],
            'name': params['compound_set_name'],
            'description': 'Compound Set produced from %s' % file_name,
            'compounds': compounds,
        }

        output = self._save_to_ws_and_report(params['workspace_id'],
                                             params['staging_file_path'],
                                             compoundset)
        #END compound_set_from_file

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method compound_set_from_file return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

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

    def compound_set_from_model(self, ctx, params):
        """
        CompoundSetFromModel
        required:
        string workspace_id
        string model_ref
        string compound_set_name
        :param params: instance of type "compoundset_from_model_params" ->
           structure: parameter "workspace_id" of String, parameter
           "model_ref" of String, parameter "compound_set_name" of String
        :returns: instance of type "compoundset_upload_results" -> structure:
           parameter "report_name" of String, parameter "report_ref" of
           String, parameter "compoundset_ref" of type "obj_ref"
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN compound_set_from_model
        self._check_param(params, ['workspace_id', 'model_ref',
                                   'compound_set_name'])
        model = self.dfu.get_objects(
            {'object_refs': [params['model_ref']]}
        )['data'][0]['data']
        compounds, undef = parse.parse_model(model)
        compoundset = {
            'id': params['compound_set_name'],
            'name': params['compound_set_name'],
            'description': 'Compound Set produced from %s, a metabolic model'
                           % model['id'],
            'compounds': compounds,
        }

        output = self._save_to_ws_and_report(params['workspace_id'],
                                             model['name'], compoundset)
        #END compound_set_from_model

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method compound_set_from_model return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def export_compoundset_as_tsv(self, ctx, params):
        """
        :param params: instance of type "ExportParams" (input and output
           structure functions for standard downloaders) -> structure:
           parameter "input_ref" of String
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN export_compoundset_as_tsv
        output = self._export_compound_set(params['input_ref'], 'tsv')
        #END export_compoundset_as_tsv

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method export_compoundset_as_tsv return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def export_compoundset_as_sdf(self, ctx, params):
        """
        :param params: instance of type "ExportParams" (input and output
           structure functions for standard downloaders) -> structure:
           parameter "input_ref" of String
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN export_compoundset_as_sdf
        output = self._export_compound_set(params['input_ref'], 'sdf')
        #END export_compoundset_as_sdf

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method export_compoundset_as_sdf return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def export_compoundset_mol2_files(self, ctx, params):
        """
        :param params: instance of type "ExportParams" (input and output
           structure functions for standard downloaders) -> structure:
           parameter "input_ref" of String
        :returns: instance of type "export_mol2_files_results" -> structure:
           parameter "packed_mol2_files_path" of String, parameter
           "comp_id_mol2_file_name_map" of mapping from String to String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN export_compoundset_mol2_files
        self._check_param(params, ['input_ref'])

        packed_mol2_files_path, comp_id_mol2_file_name_map = self._fetch_mol2_files(params['input_ref'])

        output = {'packed_mol2_files_path': packed_mol2_files_path,
                  'comp_id_mol2_file_name_map': comp_id_mol2_file_name_map}
        #END export_compoundset_mol2_files

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method export_compoundset_mol2_files return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def convert_compoundset_mol2_files_to_pdbqt(self, ctx, params):
        """
        :param params: instance of type "ExportParams" (input and output
           structure functions for standard downloaders) -> structure:
           parameter "input_ref" of String
        :returns: instance of type "convert_mol2_files_results" -> structure:
           parameter "packed_pdbqt_files_path" of String, parameter
           "comp_id_pdbqt_file_name_map" of mapping from String to String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN convert_compoundset_mol2_files_to_pdbqt
        self._check_param(params, ['input_ref'])

        packed_pdbqt_files_path, comp_id_pdbqt_file_name_map = self._covert_mol2_files_to_pdbqt(
                                                                            params['input_ref'])

        output = {'packed_pdbqt_files_path': packed_pdbqt_files_path,
                  'comp_id_pdbqt_file_name_map': comp_id_pdbqt_file_name_map}
        #END convert_compoundset_mol2_files_to_pdbqt

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method convert_compoundset_mol2_files_to_pdbqt return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def fetch_mol2_files_from_zinc(self, ctx, params):
        """
        :param params: instance of type "FetchZINCMol2Params" -> structure:
           parameter "workspace_id" of String, parameter "compoundset_ref" of
           type "obj_ref", parameter "over_write" of Long
        :returns: instance of type "compoundset_upload_results" -> structure:
           parameter "report_name" of String, parameter "report_ref" of
           String, parameter "compoundset_ref" of type "obj_ref"
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN fetch_mol2_files_from_zinc
        self._check_param(params, ['workspace_id', 'compoundset_ref'], opt_param=['over_write'])
        over_write = params.get('over_write', False)
        compoundset = self.dfu.get_objects(
                                {'object_refs': [params['compoundset_ref']]})['data'][0]['data']

        compoundset_copy = copy.deepcopy(compoundset)

        count = 0
        for compound in compoundset_copy.get('compounds'):
            if not compound.get('mol2_handle_ref') or over_write:
                temp_dir = os.path.join(self.scratch, str(uuid.uuid4()))
                os.mkdir(temp_dir)
                mol2_file_path = os.path.join(temp_dir, compound.get('id'))
                inchikey = compound.get('inchikey')
                if zinc_db_util.inchikey_to_mol2(inchikey, mol2_file_path):
                    handle_id = self.dfu.file_to_shock({'file_path': mol2_file_path,
                                                        'make_handle': True})['handle']['hid']
                    compound['mol2_handle_ref'] = handle_id
                    compound['mol2_source'] = 'ZINC15'
                    count += 1
                else:
                    logging.warning('Cannot find Mol2 file from ZINC for {}'.format(inchikey))

        if count:
            message = 'Successfully fetched {} Mol2 files from ZINC database'.format(count)
        else:
            message = 'Fetched 0 Mol2 files from ZINC database. The CompoundSet object remains unchanged.'

        output = self._save_to_ws_and_report(
                    params['workspace_id'], '', compoundset_copy,
                    message=message)

        #END fetch_mol2_files_from_zinc

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method fetch_mol2_files_from_zinc return value ' +
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
