import sys
import os
import re
import mmap
import glob
from csv import DictReader
import openbabel
from openbabel import *

# This program is for collecting data from nwchem output files: Neeraj Kumar's group working with Yarrowia Molecules
# ***** This code is not optimized *****

# Core Functions
#################################################################
#gets the number of atoms      could grab the last xyz file ********
def getNumberOfAtoms(open_file):
    #gets the line with "Charge" as this only occurs twice, and the first time is above the list of atoms
   
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "Charge" not in line:
        line = open_file.readline()
    open_file.readline() # moves to the horizontal --- ------
    line = open_file.readline() # moves to the first atom

    #counts the number of atoms
    atomCount = 0;
    while not line.isspace():
        line = open_file.readline() # moves to the next atom
        atomCount+=1
    return atomCount

#################################################################
#gets internal energy at 0 K *******
def getInternalEnergy0K(open_file):
    open_file.seek(0,os.SEEK_SET) # goes to end of the file
    lineList = open_file.readlines()
    for line in reversed(lineList):
        if "Total DFT energy" in line:
#        if "sol phase energy" in line:
            break
    result = re.findall("[+-]?\d+\.\d+",line)
    return result[0]

######## *********
def getMullikenCharge(open_file, nAtoms):
    open_file.seek(0,os.SEEK_SET)
    lineList = open_file.readlines()
    atom_coords = []

    xyz_ind =  len(lineList) - 1 - lineList[::-1].index('                         Geometry "geometry" -> "geometry"\n')
    atoms = nAtoms #int(getNumberOfAtoms(inchi_key))
    for i in range(xyz_ind+7, xyz_ind+7+atoms):
        atom_coords.append(lineList[i].split())

    ind = len(lineList) - 1 - lineList[::-1].index('      Total Density - Mulliken Population Analysis\n')
    atoms = nAtoms #int(getNumberOfAtoms(inchi_key))
    for i in range(ind+5, ind+5+atoms):
        atom_coords.append(lineList[i].split()[3])

    output = "%s\t%12s\t%12s\t%12s\t%12s" % (atom_coords[0][1], atom_coords[0][3], atom_coords[0][4], atom_coords[0][5],str(float(atom_coords[0][2]) - float(atom_coords[atoms+0]))+'\n')
    for l in range(1, atoms):
        if l < atoms-1:
           output += "%s\t%12s\t%12s\t%12s\t%8s" % (atom_coords[l][1], atom_coords[l][3], atom_coords[l][4], atom_coords[l][5], str(float(atom_coords[l][2]) - float(atom_coords[atoms+l]))+'\n')
        else:
           output += "%s\t%12s\t%12s\t%12s\t%8s" % (atom_coords[l][1], atom_coords[l][3], atom_coords[l][4], atom_coords[l][5], str(float(atom_coords[l][2]) - float(atom_coords[atoms+l])))
    return output

def calculate(InChI_key):
    os.chdir('../..')
    cwd = os.getcwd()

    dft_dir    =  InChI_key + '/' + 'dft'
    prop_file = InChI_key + '_prop.xyz'  


    matches = []
    os.chdir(dft_dir)
    pattern = InChI_key+'_nwchem.out'
    index=pattern.split('.')[1]
    matches.append(pattern)
    for j in matches:
        with open(j,'r') as f:

            outputData = ("%s%s\t%s" % (str(nAtoms)+'\n',InChI_key,E0K))
            outputData += "\n" + getMullikenCharge(f, nAtoms) + "\n"# + getFrequencies(f)
            with open(prop_file, 'w') as output:
                output.write(outputData)
            
#
#    Convert propert included xyz file to mol2 file
#
            xyz_all = [i for i in os.listdir('./') if i.endswith('_prop.xyz')]
            #print(xyz_all)
            for i in xyz_all:
                ii = open(i, 'r')
                ii.readline()
                inchikey = ii.readline().split('\t')[0]
                os.system('babel -i xyz %s  -o mol2 %s' %(i, inchikey.rstrip()+'_prop.mol2'))

#
#   Append Mulliken charges to mol2 file
#
            for xyz in xyz_all:
                key  = xyz.split('.')[0]
                mol2 = str(key) + '.mol2'

                xyzfile = open(xyz, 'r')
                mol2file = open(mol2, 'r')

                num_lines = sum(1 for line in open(mol2))
                natoms = xyzfile.readline().rstrip()
                xyzfile.readline()

                popn =  []
                for i in range(int(natoms)):
                    xyz = xyzfile.readline().rstrip()
                    xy = xyz.strip().split('\t')
                    abc = [x for x in xy if x]
                    popn.append(abc[4])

                mol2_key = key.split('_')[0]
                mol2_updated = mol2_key +'_Mulliken.mol2'
                ff = open(mol2_updated, 'w+')

                for i in range(7):
                    abc = mol2file.readline()
                    ff.write(abc)

                for i in range(7, 7+int(natoms)):
                    geometry = mol2file.readline().rstrip()
                    ff.write(str((' '.join(geometry.split(' ')[:-1])) + popn[i-7] + '\n'))

                for i in range(7+int(natoms), int(num_lines)):
                    cde = mol2file.readline()
                    ff.write(cde)

                ff.close()
            os.chdir('../../')

