import os
import sys


symbol_to_atomic_num = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C':6, 'N': 7, 'O':8, 'F':9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl':17,'Ar': 18, 'K': 19, 'Ca':20, 'Sc': 21, 'Ti': 22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga': 31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36, 'Rb':37, 'Sr':38, 'Y':39, 'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49, 'Sn':50}

def dft_nw_file(InChI_key):
    dft_dir    =  InChI_key + '/' + 'dft' 
    nwfile     = InChI_key +'.nw'
    chargefile = InChI_key +'.charge'
    xyzfile    = InChI_key +'.xyz'

    os.chdir(dft_dir)
    path    =  os.getcwd()
    with open(nwfile, 'r') as file:
         lines = file.read()
         lines = lines.replace('INCHIKEY', InChI_key)
         lines = lines.replace('PATH', path)
#
         with open(chargefile, 'r') as chgfile:
              for i in chgfile.readlines():
                  CHARGE = i.split('\t')[1]
#
         lines = lines.replace('CHARGE', CHARGE)

#

    with open(xyzfile, 'r') as xyz:
         natom = xyz.readline()
         comment = xyz.readline()
         sum_elctrons = 0
         for i in range(0, int(natom.rstrip())):
             symbol = xyz.readline().split()[0]
             sum_elctrons += float(symbol_to_atomic_num[symbol])

         total_elctrons = sum_elctrons + float(CHARGE)
         if total_elctrons % 2 == 0:
             lines = lines.replace('multiplicity', '3')
         else:
             lines = lines.replace('multiplicity', '2 \n  odft')

#

    with open(nwfile, 'w') as file:
         file.write(lines)
    os.chdir('../../')

def dft_sbatch_file(InChI_key):
    sbatchFile = InChI_key +'.sbatch' 
    dft_dir  =  InChI_key + '/' + 'dft'

    os.chdir(dft_dir)
    path    =  os.getcwd()
    with open(sbatchFile, 'r') as file:
         lines = file.read()
         lines = lines.replace('INCHIKEY', InChI_key)
         lines = lines.replace('PATH', path)
    with open(sbatchFile, 'w') as file:
         file.write(lines)
    os.chdir('../../')

if __name__ == '__main__':
    InChI_key = sys.argv[1]
    dft_nw_file(InChI_key)
    dft_sbatch_file(InChI_key)


