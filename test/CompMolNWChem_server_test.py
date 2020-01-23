# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from unittest.mock import patch

from CompMolNWChem.CompMolNWChemImpl import CompMolNWChem
from CompMolNWChem.CompMolNWChemServer import MethodContext
from CompMolNWChem.authclient import KBaseAuth as _KBaseAuth

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class CompMolNWChemTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('CompMolNWChem'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'CompMolNWChem',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = CompMolNWChem(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsId(self):
        if hasattr(self.__class__, 'wsId'):
            return self.__class__.wsId
        suffix = int(time.time() * 1000)
        wsName = "test_CompoundSetUtils_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsId = ret[0]
        return ret[0]

    #def getImpl(self):
    #    return self.__class__.serviceImpl
    
   # @staticmethod
   # def fake_staging_download(params):
    #    scratch = '/kb/module/work/tmp/'
     #   inpath = params['staging_file_subdir_path']
      #  shutil.copy('/kb/module/test/'+inpath, scratch+inpath)
       # return {'copy_file_path': scratch+inpath}
    
    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa

    #@patch.object(DataFileUtil, "download_staging_file",
    #              new=fake_staging_download)
    #def test_compound_set_from_file_tsv(self):
    #    params = {'workspace_id': self.getWsId(),
    #              'staging_file_path': 'test_compounds.tsv',
    #              'compound_set_name': 'tsv_set_1',
    #              'mol2_staging_file_path': 'mol2_files_missing_comp.zip'}
    #    ret = self.getImpl().compound_set_from_file(self.getContext(), params)[0]
    #    assert ret and ('report_name' in ret)

    def test_your_method(self):

        ret = self.serviceImpl.run_CompMolNWChem(self.ctx, {'workspace_name': self.wsName,'workspace_id':self.getWsId(),
                                                                 'Input_File':'test_compounds.tsv','calculation_type':'energy'})

        print("Output")
        print (ret)
