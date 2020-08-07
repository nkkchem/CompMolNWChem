import sys
import os
import re
import mmap
import glob
import openbabel as obabel

# This program is for collecting data from nwchem output files: Neeraj Kumar's group working with Yarrowia Molecules
# ***** This code is not optimized *****

# Core Functions
#################################################################
#gets the number of atoms      could grab the last xyz file
def getNumberOfAtoms(inchi_key):
#    os.chdir(inchi_key)
    xyz_list = []
    xyz_list += [each for each in os.listdir('./') if each.endswith('.xyz')]

    start_geom = sorted(xyz_list)[0]

    with open(start_geom, 'r') as geom_file:
         lines = geom_file.readlines()
         atomCount = lines[0].rstrip()
#    os.chdir('../')
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
#gets internal energy at 0 K
def getInternalEnergy0K(open_file):
    open_file.seek(0,os.SEEK_SET) # goes to end of the file
    lineList = open_file.readlines()
    for line in reversed(lineList):
#        if "Total DFT energy" in line:
        if "total free energy in solvent including G(SMD-CDS)" in line:
            break
    result = re.findall("[+-]?\d+\.\d+",line)
    return result[0]

#################################################################
#gets Enthalpy at 298.15 K. Does this by taking Internal Energy and adding it to the thermal correction to Enthalpy
def getEnthalpy(internalEnergy, ThermalCorrection):
    return float(ThermalCorrection) + internalEnergy        # [0]= kcal  [1]= au

#################################################################
# gets free energy at 298.15 K. free energy (G) is Correction + Enthalpy*627.509469 - Temperature*(Total Entropy / (1000))
# this gives results in kcal
def getFreeEnergy(internalEnergy, ThermalCorrection, totalEntropy):
    return float(ThermalCorrection) + internalEnergy*627.509469 - (298.15*(float(totalEntropy[0]) / (1000)))        # [0]= kcal  [1]= au

#################################################################
#gets the solvation energy
def getSolvationEnergy(open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "delta internal energy" not in line:
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

#################################################################
#gets the Element type, coordinate (x,y,z)(Angstrom), and Mulliken partial charge (e) of the atom)
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
           output += "%s\t%12s\t%12s\t%12s\t%8s" % (atom_coords[l][1], atom_coords[l][3], atom_coords[l][4], atom_coords[l][5], str(float(atom_coords[l][2]) - float(atom_coords[atoms+l]))
+'\n')
        else:
           output += "%s\t%12s\t%12s\t%12s\t%8s" % (atom_coords[l][1], atom_coords[l][3], atom_coords[l][4], atom_coords[l][5], str(float(atom_coords[l][2]) - float(atom_coords[atoms+l]))
)
    return output

#################################################################
#gets the frequencies (3na-5 or 3na-6)
def getFrequencies(open_file, nAtoms):
    open_file.seek(0,os.SEEK_SET)
    lineList = open_file.readlines()
    output = ""

    ind = len(lineList) - 1 - lineList[::-1].index(' Normal Eigenvalue ||           Projected Infra Red Intensities\n')
    atoms = nAtoms #int(getNumberOfAtoms(inchi_key))
    for i in range(ind+3, ind+3+3*atoms):
        output += lineList[i].split()[1] + '\t'
    return output
#################################################################
#gets the starting and optimized SMILES
def getSMILES(inchi_key):
    dft = inchi_key + '/dft'
#    os.chdir(dft)
    xyz_files = []
    xyz_files += [each for each in os.listdir('./') if each.endswith('.xyz')]
#
    converged_geom = sorted(xyz_files)[-2]
    os.system('obabel %s -O converged.smiles' % (converged_geom))
#
    file = open('converged.smiles', 'r')
    smiles_lines = file.readlines()
    converged_smiles = smiles_lines[0].split('\t')[0]
#    os.chdir('../../')
    return ('%s' %(converged_smiles))

#################################################################
#gets the starting and optimized InChI
def getInChI(inchi_key):
    dft = inchi_key + '/dft'
#    os.chdir(dft)
    xyz_files = []
    xyz_files += [each for each in os.listdir('./') if each.endswith('.xyz')]
#
    converged_geom = sorted(xyz_files)[-2]
    os.system('obabel %s -O converged.inchi' % (converged_geom))
#
    file = open('converged.inchi', 'r')
    inchi_lines = file.readlines()
    converged_inchi = inchi_lines[0].rstrip()
#    os.chdir('../../')
    return (('%s' %(converged_inchi)))

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

def ThermalCorrectionEnthalpy(open_file):
    open_file.seek(0,os.SEEK_SET) # goes to beginning of the file
    line = open_file.readline()
    while "Thermal correction to Enthalpy" not in line:
        line = open_file.readline()
    correction = re.findall("[+-]?\d+\.\d+",line)
    while "Total Entropy" not in line:
        line = open_file.readline()
    totalEntropy = re.findall("[+-]?\d+\.\d+",line)
    return correction[0], totalEntropy

#################################################################
# call all functions to extract properties to a file

def calculate(InChI_key):
    dft_dir    =  InChI_key + '/' + 'dft'

    matches = []
    os.chdir(dft_dir)
    pattern = InChI_key + '.out'
    property_output_file = InChI_key + '_properties.dat'
    index=pattern.split('.')[1]
    matches.append(pattern)
    for j in matches:
#
        with open(j,'r') as f:
            outputData = ("nAtoms\tindex\tmu\tHOMO\tLUMO\tgap\tZPVE\tE\tE0K\tH\tG\tS\tCv\n")

            nAtoms = int(getNumberOfAtoms(InChI_key))

            zero = float(getZeroPointCorrection(f))
            A,B,C=getRotationalConstants(f)
            ZPVE = getZeroPointVibrationEnergy(zero,f)
            ThermalCorr = ThermalCorrectionEnergy(f)
            EnthalpyCorr, Entropy = ThermalCorrectionEnthalpy(f)
            cV = getHeatCapacity(f)
            freq = getFrequencies(f, nAtoms)
            Mulliken = getMullikenCharge(f, nAtoms)

            mu = getDipoleMoment(f)
            HOMO, LUMO = getHOMO_LUMO(f)
            gap = getGap(LUMO,HOMO)
            E0K = getInternalEnergy0K(f)
            internalEnergy = float(E0K)
            E = getInternalEnergy(internalEnergy,ThermalCorr)
            H = getEnthalpy(internalEnergy,ThermalCorr)
            G = getFreeEnergy(internalEnergy, EnthalpyCorr, Entropy)/627.509469
            print(G)
            tEsolvation = getSolvationEnergy(f)
            SMILES = getSMILES(InChI_key)
            InChI = getInChI(InChI_key)
            print('\n')

            outputData = ("%s\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (nAtoms, InChI_key, index, A, B, C, mu, HOMO, LUMO, gap, ZPVE, E0K, E, H, G, cV))
            outputData += "\n" + Mulliken + "\n" + freq + '\n'
            outputData += SMILES + "\n"
            outputData += InChI
            with open(property_output_file, 'w') as output:
                output.write(outputData)
            os.chdir('../../')


if __name__ == '__main__':
    InChI_key = sys.argv[1]
    calculate(InChI_key)
