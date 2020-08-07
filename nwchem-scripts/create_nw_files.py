import os
import sys

def dft_nw_file(InChI_key):
    dft_dir    =  InChI_key + '/' + 'dft' 
    nwfile     = InChI_key +'.nw'
    chargefile = InChI_key +'.charge'

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
         lines = lines.replace('CHARGE', '-1')

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




