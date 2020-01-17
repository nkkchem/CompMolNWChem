"""
Script to read INCHIKEY and INCHI-STRING to generate .xyz, .mol, .nw and .sbatch
files and to submit jobs

Needs test.csv file with data

--> excute these two first in shell to activate rdkit

conda create -c rdkit -n my-rdkit-env rdkit
conda activate my-rdkit-env

"""

from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit
from rdkit.Chem import Draw
import os
import sys
import subprocess
from pybel import *
from csv import DictReader

# Generate dft folder from inchi string

#os.system('rm -rf Z* A* K* D* F* P* W* X*')

def inchi_to_dft(InChI_key,InChIes):
    os.chdir('../../../simulation')
    dd = dict(zip(InChI_key, InChIes))
    #print('dd:',dd)
    #print('dd_items:',dd.items())
    #dd = dict(InChI_key=InChIes)
    #print(dd)
    for key, value in dd.items():
        os.mkdir(key)
        os.chdir(key)
        #print('key:',key)
        #print('value:',value)
        m4 = Chem.inchi.MolFromInchi(value)
        #print('m4:',m4)
        AllChem.Compute2DCoords(m4)
        m5 = Chem.AddHs(m4)
        #print('m5',m5)
        if m5.GetNumAtoms() > 110:
            AllChem.EmbedMolecule(m5, useRandomCoords=True)
        else:
            AllChem.EmbedMolecule(m5)
    
        AllChem.MMFFOptimizeMolecule(m5)
        ii = Chem.MolToMolBlock(m5).splitlines()
    
        os.mkdir('dft')   # create dft files
        os.chdir('dft')
    
        f=open(str(key)+".xyz", 'w')
        f.write(str(m5.GetNumAtoms()))
        f.write('\n\n')
        for i in range(4, 4+m5.GetNumAtoms()):
            jj=' '.join((ii[i].split()))
            kk=jj.split()
            f.write("{}\t{:>8}  {:>8}  {:>8}\n".format(kk[3], kk[0], kk[1], kk[2]))
        f.close()

    #
    #   write .nw file
    #
        nw=open(str(key)+".nw", 'w')
        nw.write("title " +'"'+key+'"'+'\n')
        nw.write('start '+key+'\n\n')
        nw.write('memory global 1600 mb heap 100 mb stack 600 mb'+'\n\n')
        nw.write('permanent_dir ' + os.getcwd()+'\n')
        nw.write('#scratch_dir /scratch'+'\n\n')
        nw.write('echo'+'\n')
        nw.write('print low'+'\n\n')
        nw.write('charge ' + str(rdkit.Chem.rdmolops.GetFormalCharge(m5))+'\n')
        nw.write('geometry noautoz noautosym'+'\n')
        nw.write('load '+ os.getcwd()+'/'+key+'.xyz'+'\n')
        nw.write('end'+'\n')
        nw.write('basis'+'\n')
        nw.write('* library 6-31G*'+'\n')
        nw.write('end'+'\n\n')
        nw.write('driver'+'\n')
        nw.write(' maxiter 50'+'\n')
        nw.write('xyz FXEUKVKGTKDDIQ-UWVGGRQHSA-M_geom'+'\n')
        nw.write('end'+'\n\n')
        nw.write('set lindep:n_dep 0'+'\n\n')
        nw.write('dft'+'\n')
        nw.write('  maxiter 100'+'\n')
        nw.write('  xc b3lyp'+'\n')
        #nw.write('  disp vdw 3'+'\n')
        nw.write('  mulliken'+ '\n')
        nw.write('  print "mulliken ao"'+'\n')
        nw.write('  print "final vectors analysis"'+'\n')
        nw.write('end'+'\n\n')
        nw.write('task dft optimize ignore'+'\n')
        #nw.write('task dft freq'+'\n\n')
        #nw.write('cosmo'+'\n')
        #nw.write(' dielec 80.4'+'\n')
        #nw.write(' lineq  0'+'\n')
        #nw.write('end'+'\n\n')
        #nw.write('task dft energy'+'\n\n')
        nw.write('property'+'\n')
        nw.write(' dipole'+'\n')
        nw.write('end'+'\n\n')
        nw.write('task dft property'+'\n')
        nw.close()

        nwfile = str(key)+".nw"
        #print('nwfile:',nwfile)
        os.system("mpirun -np 2 --allow-run-as-root nwchem *.nw > nwchem.out 2>error")
        #os.system('ls')
        #os.system('tail -20 *.out')
        cwd = os.getcwd()
        #print('After Submission:',cwd)
        os.chdir('../..')







        
# Reads InChI and InChI-key from .csv file
#with open("test.csv") as f:
#     InChI_key = [row["InChI-Key"].split('InChIKey=')[1] for row in DictReader(f)] 
#InChI_key = ['FFQKYPRQEYGKAF-UHFFFAOYSA-L']
#with open("test.csv") as f:
#     InChIes   = [row["InChI"] for row in DictReader(f)]
#InChIes   = ["InChI=1S/CH4NO5P/c2-1(3)7-8(4,5)6/h(H2,2,3)(H2,4,5,6)/p-2"]

#inchi_to_dft(InChI_key, InChIes)


