from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit
from rdkit.Chem import Draw
import os
import sys
import csv
from pybel import *
#from csv import DictReader

def inchiTo3D(args):
    """Reads the InChI string from the .inchi file, which is in a {inchi-key} molecule folder.
       Converts the InChI string to .mol, .mol2, .xyz, and .png files
    """
    with open(args['inchifile']) as f:
        smiles = f.readline()    # read the InChI string from the .inchi file
    m4 = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(m4)
    m5 = Chem.AddHs(m4)
    if m5.GetNumAtoms() > 110:
        AllChem.EmbedMolecule(m5, useRandomCoords=True)
    else:
        AllChem.EmbedMolecule(m5)

    AllChem.MMFFOptimizeMolecule(m5)
    ii = Chem.MolToMolBlock(m5).splitlines()

# create .xyz file    
    with open(args['xyzfile'],'w') as f:
        f.write(str(m5.GetNumAtoms()))
        f.write('\n\n')
        for i in range(4, 4+m5.GetNumAtoms()):
            jj=' '.join((ii[i].split()))
            kk=jj.split()
            f.write("{}\t{:>8}  {:>8}  {:>8}\n".format(kk[3], kk[0], kk[1], kk[2]))


#create .charge file
    with open(args['chargefile'],'w') as ff:
         ff.write(str(args['chargefile']).split('/')[0]+'\t')
         ff.write(str(rdkit.Chem.rdmolops.GetFormalCharge(m5)) + '\t')
         ff.write(str(m5.GetNumAtoms()))

# create .mol file

    with open(args['molfile'],'w') as f:
        m4 = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(m4)
        m5 = Chem.AddHs(m4)
    #
        if m5.GetNumAtoms() > 110:
            AllChem.EmbedMolecule(m5, useRandomCoords=True)
        else:
            AllChem.EmbedMolecule(m5)

        AllChem.MMFFOptimizeMolecule(m5)
        f.write(Chem.MolToMolBlock(m5))

# create .png file
        Draw.MolToFile(m5, args['pngfile'])


# create .mol2 file  
#    mol22 = readfile('mol', args['molfile']).__next__()
#    output = Outputfile('mol2',args['mol2file'])
#    output.write(mol22)

if __name__ == '__main__':
    args={}
    print("check")
    args['inchifile'] = sys.argv[1]
    args['xyzfile'] = sys.argv[2]
    args['molfile'] = sys.argv[3]
#    args['mol2file'] = sys.argv[4]
    args['pngfile'] = sys.argv[4]
    args['chargefile'] = sys.argv[5]

    print(args)
    inchiTo3D(args)
