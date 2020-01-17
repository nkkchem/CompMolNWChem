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
#gets the Consecutive, 1-based integer identifier of molecule      job output number
def get1BasedIdentifier(fileName):
    result = fileName.split('.')
    return(result[2])

#################################################################
# gets the dipole moment in Debye
def getDipoleMoment(open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "Dipole moment" and "Debye" not in line:
        line = open_file.readline()
    result = re.findall("[+-]?\d+\.\d+",line)
    return (result[0])

#################################################################
#gets the energy of the Highest Occupied Molecular Orbital (HOMO) and the Lowest Unoccupied Molecular Orbital (LUMO)
def getHOMO_LUMO(open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    lineList = open_file.readlines()
    for line in reversed(lineList):
        if "Occ=0" in line:
            LUMO = re.findall("[+-]?\d+\.\d+[D+-]*\d+",line)
            LUMO = [i.replace('D', 'E') for i in LUMO]
#            LUMO = re.findall("[+-]?\d+\.\d+",line) # will rewrite this one each time Occ=0 is found until Occ=2 is found
        elif "Occ=2" in line:
            HOMO = re.findall("[+-]?\d+\.\d+[D+-]*\d+",line)
            HOMO = [i.replace('D', 'E') for i in HOMO]
#            HOMO = re.findall("[+-]?\d+\.\d+",line)
            break
    return (float(HOMO[1]), float(LUMO[1]))

#################################################################
#gets the difference between LUMO and HOMO
def getGap(LUMO,HOMO):
    return LUMO - HOMO

#################################################################
#gets zero point vibrational energy
def getZeroPointVibrationEnergy(correction, open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "Total Entropy" not in line:
        line = open_file.readline()
    open_file.readline()
    open_file.readline()
    line = open_file.readline()
    result = re.findall("[+-]?\d+\.\d+",line)
    return (float(result[0]) / (1000*627.509469)) + correction

#################################################################
#gets internal energy at 298 K
def getInternalEnergy(internal,Thermalcorrection):
    return internal + float(Thermalcorrection)

def ThermalCorrectionEnergy(open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "Thermal correction to Energy" not in line:
        line = open_file.readline()
    correction = re.findall("[+-]?\d+\.\d+",line)
    return correction[1]
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

#################################################################
#gets Enthalpy at 298.15 K. Does this by taking Internal Energy and adding it to the thermal correction to Enthalpy
def getEnthalpy(internalEnergy, open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "Thermal correction to Enthalpy" not in line:
        line = open_file.readline()
    correction = re.findall("[+-]?\d+\.\d+",line)
    return float(correction[1]) + internalEnergy        # [0]= kcal  [1]= au

#################################################################
# gets free energy at 298.15 K. free energy (G) is Correction + Enthalpy*627.509469 - Temperature*(Total Entropy / (1000))
# this gives results in kcal
def getFreeEnergy(internalEnergy,open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "Thermal correction to Enthalpy" not in line:
        line = open_file.readline()
    correction = re.findall("[+-]?\d+\.\d+",line)
    while "Total Entropy" not in line:
        line = open_file.readline()
    totalEntropy = re.findall("[+-]?\d+\.\d+",line)
    return float(correction[0]) + internalEnergy*627.509469 - (298.15*(float(totalEntropy[0]) / (1000)))        # [0]= kcal  [1]= au
#################################################################
#gets the solvation energy
def getSolvationEnergy(open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "solvation energy" not in line:
        line = open_file.readline()
    result = re.findall("[+-]?\d+\.\d+",line)
    return float(result[0]) * 627.509469       # 627.509469 is a conversion factor

#################################################################
#gets heat capacity at 298.15 K
def getHeatCapacity(open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "heat capacity" not in line:
        line = open_file.readline()
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

#################################################################
#gets the frequencies (3na-5 or 3na-6)
def getFrequencies(open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "Normal Eigenvalue" not in line:
        line = open_file.readline()
    open_file.readline()
    open_file.readline()
    line = open_file.readline()
    output = ""
    while "-" not in line:
        result = [x.strip() for x in line.split()]
        output += result[1] + " \t"
        line = open_file.readline() # moves to the next atom
    return output
#################################################################
#gets the SMILES
def getSMILES(open_file):
    open_file.seek(0, os.SEEK_SET) # goes to beginning of the file

#################################################################
#gets the InChI
def getInChI(open_file):
    open_file.seek(0 ,os.SEEK_SET) # goes to beginning of the file

#================================================================#
#Support Functionality
def getZeroPointCorrection(open_file):
    open_file.seek(0 ,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "Zero-Point" not in line:
        line = open_file.readline()
    result = re.findall("[+-]?\d+\.\d+",line)
    return result[1]        # [0]= kcal  [1]= au

def getRotationalConstants(open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "Rotational Constants" not in line:
        line = open_file.readline()
    line = open_file.readline()
    A_line = open_file.readline()
    A = re.findall('\d.\d*', A_line)[0]
    B_line = open_file.readline()
    B = re.findall('\d.\d*', B_line)[0]
    C_line = open_file.readline()
    C = re.findall('\d.\d*', C_line)[0]
    return (float(A), float(B), float(C))
#################################################################

def calculate(InChI_key):
    os.chdir('../..')
    cwd = os.getcwd()

    dft_dir    =  InChI_key + '/' + 'dft'
    prop_file = InChI_key + '_prop.xyz'  


    matches = []
    os.chdir(dft_dir)
    pattern = 'nwchem.out'
    #pattern = InChI_key + '.out'
    index=pattern.split('.')[1]
    matches.append(pattern)
    for j in matches:
        with open(j,'r') as f:
#            outputData = ("nAtoms\tindex\tmu\tHOMO\tLUMO\tgap\tZPVE\tE\tE0K\tH\tG\tS\tCv\n")

            #zero = float(getZeroPointCorrection(f))
            #nAtoms = getNumberOfAtoms(f)
            #A,B,C=getRotationalConstants(f)
            #mu = getDipoleMoment(f)
            #HOMO, LUMO = getHOMO_LUMO(f)
            #gap = getGap(LUMO,HOMO)
            #ZPVE = getZeroPointVibrationEnergy(zero,f)
            #E0K = getInternalEnergy0K(f)
            #internalEnergy = float(E0K)
            #ThermalCorr = ThermalCorrectionEnergy(f)
            #E = getInternalEnergy(internalEnergy,ThermalCorr)
            #H = getEnthalpy(internalEnergy,f)
            #G = getFreeEnergy(internalEnergy,f)
            #tEsolvation = getSolvationEnergy(f)
            #cV = getHeatCapacity(f)
            #geom = getMullikenCharge(f, nAtoms)
            outputData = ("%s%s\t%s" % (str(nAtoms)+'\n',InChI_key,E0K))
            outputData += "\n" + getMullikenCharge(f, nAtoms) + "\n"# + getFrequencies(f)
            with open(prop_file, 'w') as output:
                output.write(outputData)
            
#
#    Convert propert included xyz file to mol2 file
#
            xyz_all = [i for i in os.listdir('./') if i.endswith('_prop.xyz')]
            print(xyz_all)
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
            #os.system('ls')
            #os.system('pwd')
            #os.system('more *_Mulliken.mol2')
            os.chdir('../../')
            


#if __name__ == '__main__':
#    InChI_key = sys.argv[1]
#    calculate(InChI_key)

#    with open("test.csv") as f:
#         InChI_key = [row["InChI-Key"].split('InChIKey=')[1] for row in DictReader(f)]

#    for inchi_key in InChI_key:
#        calculate(inchi_key)


