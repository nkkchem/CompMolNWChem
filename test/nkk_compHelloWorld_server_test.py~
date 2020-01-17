# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from nkk_compHelloWorld.nkk_compHelloWorldImpl import nkk_compHelloWorld
from nkk_compHelloWorld.nkk_compHelloWorldServer import MethodContext
from nkk_compHelloWorld.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class nkk_compHelloWorldTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('nkk_compHelloWorld'):
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
                            {'service': 'nkk_compHelloWorld',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = nkk_compHelloWorld(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "nkk_compHelloWorld_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_your_method(self):

        ret = self.serviceImpl.run_nkk_compHelloWorld(self.ctx, {'workspace_name': self.wsName,
                                                                 'InChI_key': 'CKFGINPQOCXMAZ-UHFFFAOYSA-N','InChIes':'InChI=1S/H2O/h1H2'})
        print("Output")
        print (ret)
# def test_special_smiles(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
#        smiles_string = 'N[C@H](C)C(=O)O'
#        ret = self.serviceImpl.run_nkk_compHelloWorld(self.ctx, {'workspace_name': self.wsName,
#                                                             'smiles_string': smiles_string})
#        print("Output")
#        print (ret)
