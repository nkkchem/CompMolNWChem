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
