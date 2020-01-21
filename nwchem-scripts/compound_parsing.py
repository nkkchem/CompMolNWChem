import csv
import json
import logging
import os
from rdkit.Chem import AllChem, Descriptors
from installed_clients.DataFileUtilClient import DataFileUtil


def _make_compound_info(mol_object):

    return {
        'smiles': AllChem.MolToSmiles(mol_object, True),
        'inchikey': AllChem.InchiToInchiKey(AllChem.MolToInchi(mol_object)),
        'mass': Descriptors.MolWt(mol_object),
        'exactmass': AllChem.CalcExactMolWt(mol_object),
        'formula': AllChem.CalcMolFormula(mol_object),
        'charge': AllChem.GetFormalCharge(mol_object),
        'fingerprints': {
            'maccs': dict([(str(x), 1) for x in AllChem.GetMACCSKeysFingerprint(mol_object).GetOnBits()]),
            'rdkit': dict([(str(x), 1) for x in AllChem.RDKFingerprint(mol_object).GetOnBits()]),
        },
        'dblinks': {},
    }


def read_tsv(file_path, structure_field='structure',
             inchi_path='/kb/module/data/Inchikey_IDs.json', mol2_file_dir=None,
             callback_url=None):

    inchi_dict = json.load(open(inchi_path))
    cols_to_copy = {'name': str, 'deltag': float, 'deltagerr': float}
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    compounds = []
    w = csv.DictReader(open(file_path), dialect='excel-tab')
    for i, line in enumerate(w):
        user_id = line.get('id')
        mol2_source = line.get('mol2_source')
        handle_id = None
        if user_id and mol2_file_dir:
            if not mol2_source:
                raise ValueError('Please indicate mol2 file source in TSV file')
            mol2_file_path = None
            for root, dirs, files in os.walk(mol2_file_dir):
                for file in files:
                    if os.path.splitext(file)[0] == user_id:
                        logging.info('Found a matching mol2 file {} for compound {}'.format(str(file), user_id))
                        mol2_file_path = os.path.join(root, str(file))

            if mol2_file_path:
                dfu = DataFileUtil(callback_url)
                handle_id = dfu.file_to_shock({'file_path': mol2_file_path,
                                               'make_handle': True})['handle']['hid']
            else:
                logging.warning('Unable to find a matching mol2 file for compound: {}'.format(user_id))

        mol = None
        # Generate Mol object from InChI code if present
        if 'structure' in line:
            if "InChI=" in line['structure']:
                mol = AllChem.MolFromInchi(line['structure'])
            # Otherwise generate Mol object from SMILES string
            else:
                mol = AllChem.MolFromSmiles(line['structure'])
        elif 'smiles' in line:
            mol = AllChem.MolFromSmiles(line['smiles'])
        if not mol:
            logging.warning("Unable to Parse %s" % line[structure_field])
            continue
        comp = _make_compound_info(mol)

        if comp['inchikey'] in inchi_dict:
            comp['kb_id'] = inchi_dict[comp['inchikey']]
        else:
            comp['kb_id'] = '%s_%s' % (file_name, i + 1)

        if user_id:
            comp['id'] = user_id
        else:
            comp['id'] = comp['kb_id']

        for col in cols_to_copy:
            if col in line and line[col]:
                comp[col] = cols_to_copy[col](line[col])

        if handle_id:
            comp['mol2_handle_ref'] = handle_id
            comp['mol2_source'] = mol2_source

        compounds.append(comp)

    return compounds


def read_sdf(file_path, inchi_path='/kb/module/data/Inchikey_IDs.json', mol2_file_dir=None,
             callback_url=None):

    inchi_dict = json.load(open(inchi_path))
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    sdf = AllChem.SDMolSupplier(file_path.encode('ascii', 'ignore'))
    compounds = []
    for i, mol in enumerate(sdf):
        user_id = mol.GetPropsAsDict().get('id')
        print('Found compound ID: {}'.format(user_id))
        handle_id = None
        if user_id and mol2_file_dir:
            mol2_file_path = None
            for root, dirs, files in os.walk(mol2_file_dir):
                for file in files:
                    if os.path.splitext(file)[0] == user_id:
                        logging.info('Found a matching mol2 file {} for compound {}'.format(str(file), user_id))
                        mol2_file_path = os.path.join(root, str(file))

            if mol2_file_path:
                dfu = DataFileUtil(callback_url)
                handle_id = dfu.file_to_shock({'file_path': mol2_file_path,
                                               'make_handle': True})['handle']['hid']
            else:
                logging.warning('Unable to find a matching mol2 file for compound: {}'.format(user_id))

        comp = _make_compound_info(mol)
        comp['name'] = mol.GetProp("_Name")
        comp['mol'] = AllChem.MolToMolBlock(mol)
        if comp['inchikey'] in inchi_dict:
            comp['kb_id'] = inchi_dict[comp['inchikey']]
        else:
            comp['kb_id'] = '%s_%s' % (file_name, i + 1)

        if user_id:
            comp['id'] = user_id
        else:
            comp['id'] = comp['kb_id']

        if handle_id:
            comp['mol2_handle_ref'] = handle_id
            comp['mol2_source'] = 'user uploaded'

        compounds.append(comp)
    return compounds


def parse_model(model, struct_path='/kb/module/data/Compound_Structures.json'):
    # At the moment this function relies on cached data. when I get a chance
    # I'll figure how to pull the biochemistry object from the workspace so it
    # stays up to date
    struct_dict = json.load(open(struct_path))
    compounds = []
    undef_structures = []
    for model_comp in model['modelcompounds']:
        cid = model_comp['id'].split('_')[0]
        if cid not in struct_dict:
            undef_structures.append(cid)
            continue
        set_comp = {'id': cid,
                    'name': model_comp['name'],
                    'compound_ref': model_comp['compound_ref'],
                    'modelcompound_ref': model_comp['id'],
                    }
        mol = AllChem.MolFromInchi(str(struct_dict[cid]))
        set_comp.update(_make_compound_info(mol))
        compounds.append(set_comp)
    return compounds, undef_structures


def write_tsv(compound_set, outfile_path):
    cols = ['id', 'kb_id', 'name', 'smiles', 'inchikey', 'charge', 'formula', 'mass',
            'exactmass', 'compound_ref', 'modelcompound_ref', 'deltag',
            'deltagerr']
    outfile_path + ".tsv"
    writer = csv.DictWriter(open(outfile_path, 'w'), cols, dialect='excel-tab',
                            extrasaction='ignore')
    writer.writeheader()
    for compound in compound_set['compounds']:
        writer.writerow(compound)
    return outfile_path


def write_sdf(compound_set, outfile_path):
    no_export = {'smiles', 'fingerprints', 'dblinks'}
    outfile_path += ".sdf"
    writer = AllChem.SDWriter(open(outfile_path, 'w'))
    for compound in compound_set['compounds']:
        mol = _get_mol_from_compound(compound)
        for prop, val in compound.items():
            if prop in no_export:
                continue
            mol.SetProp(str(prop), str(val))
        writer.write(mol)
    return outfile_path


def write_mol_dir(compound_set, outfile_path, out_type='mol'):
    os.mkdir(outfile_path)
    for compound in compound_set['compounds']:
        mol = _get_mol_from_compound(compound)
        if out_type == 'mol':
            AllChem.MolToMolFile(mol, f'{outfile_path}/{compound["id"]}.mol')
        elif out_type == 'pdb':
            AllChem.MolToPDBFile(mol, f'{outfile_path}/{compound["id"]}.pdb')
        else:
            ValueError('Invalid output_format. Expects tsv, sdf, mol, or pdb')
    return outfile_path


def _get_mol_from_compound(compound):
    if compound.get('mol'):
        mol = AllChem.MolFromMolBlock(compound['mol'])
    else:
        mol = AllChem.MolFromSmiles(compound['smiles'])
        _calc_3d_coord(mol)
    return mol


def _calc_3d_coord(mol):
    AllChem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useRandomCoords=True)
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except ValueError:
        logging.warning("Unable to make 3d cords.")
    AllChem.RemoveHs(mol)
