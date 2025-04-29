import numpy as np
import pandas as pd
import sys
import argparse
import os
import re


parser = argparse.ArgumentParser(prog='vib_analysis', 
                                 description='Extracting vibrational data from ORCA, GAUSSIAN and TURBOMOLE output files and analyzing it.\n\
                                 Program already creates output files so there\'s no need to pipe the output to a file.\n',
                                 usage='%(prog)s filename [options]',
                                 epilog='Enjoy the program! :)')

parser.add_argument('filename', type=str, help='Name of the file to be analyzed. It can be either ORCA, GAUSSIAN or TURBOMOLE (aoforce) output file. \
                    When analyzing TURBOMOLE output file, the control file must be specified with the \'-c\' argument.')
parser.add_argument('-c', '--control', type=str, help='Name of the control file in TM vibrational analysis.')
parser.add_argument('-r', '--range', default=[1000,2000], type=float, nargs=2, help='Frequency range to be included in the output. Default is 1000-2000 cm^-1')
parser.add_argument('-s', '--scaling', default=1.00, type=float, help='Scaling factor for the frequencies. Default is 1.00')
parser.add_argument('-x', '--indices', type=int, nargs='*', help='List of atoms to be included in the functional group contribution analysis.\
                    When not specified, the program will try to guess the functional groups based on the atom types and proximity.\
                    The supported groups are CO, NH2, and NH.')
parser.add_argument('-int', '--intensities', action='store_true', help='Calculates Raman intensities based on the Raman activities and save them to a file.')
parser.add_argument('-xf', '--indices-file', type=str, help='File containing the list of atoms to be included in the functional group contribution analysis. \
                    Each group should be in a separate line. The atoms in the group should be separated by spaces.')
parser.add_argument('-p', '--print', action='store_true', help='Prints the contributions of functional groups to the console.')
parser.add_argument('-P', '--print-all', action='store_true', help='Prints the contributions of the functional groups and individual atoms to the console.')
parser.add_argument('-T', type=float, default=298.15, help="Temperature in Kelvin (default: 298.15)")
parser.add_argument('-WL', type=float, default=18797, help="Excitation wavelength in cm^-1 (default: 18797 cm^-1 / 532 nm)")
parser.add_argument('-sr', '--sort', type=int, nargs=1, help='Sorts the output frequencies by the contributions of inputed group (given as a index of a group).')

args = parser.parse_args()

filename = args.filename
scaling = args.scaling
min_max = args.range
temp = args.T
wl = args.WL
if args.sort != None:
    sort_idx = args.sort[0] - 1

# Constants for Raman intensity calculation
h = 6.62607015e-34  # Planck's constant, in J·s
c = 2.99792458e10   # Speed of light, in cm/s
k = 1.380649e-23    # Boltzmann constant, in J/K
C = 1e-12           # Scaling constant (arbitrary)

atom2mass = {
    'H'  : 1,     # Hydrogen
    'He' : 4,     # Helium
    'Li' : 7,     # Lithium
    'Be' : 9,     # Beryllium
    'B'  : 11,    # Boron
    'C'  : 12,    # Carbon
    'N'  : 14,    # Nitrogen
    'O'  : 16,    # Oxygen
    'F'  : 18,    # Fluorine
    'Ne' : 20,    # Neon
    'Na' : 23,    # Sodium
    'Mg' : 24,    # Magnesium
    'Al' : 27,    # Aluminum
    'Si' : 28,    # Silicon
    'P'  : 31,    # Phosphorus
    'S'  : 32,    # Sulfur
    'Cl' : 35,    # Chlorine
    'Ar' : 40,    # Argon
    'K'  : 39,    # Potassium
    'Ca' : 40,    # Calcium
    'Sc' : 45,    # Scandium
    'Ti' : 47,    # Titanium
    'V'  : 51,    # Vanadium
    'Cr' : 52,    # Chromium
    'Mn' : 55,    # Manganese
    'Fe' : 56,    # Iron
    'Ni' : 59,    # Nickel
    'Cu' : 63,    # Copper
    'Zn' : 65,    # Zinc
    'Ga' : 69,    # Gallium
    'Ge' : 73,    # Germanium
    'As' : 75,    # Arsenic
    'Se' : 79,    # Selenium
    'Br' : 80,    # Bromine
    'Kr' : 84,    # Krypton
    'Rb' : 85,    # Rubidium
    'Sr' : 87,    # Strontium
    'Y'  : 89,    # Yttrium
    'Zr' : 91,    # Zirconium
    'Nb' : 93,    # Niobium
    'Mo' : 96,    # Molybdenum
    'Tc' : 98,    # Technetium
    'Ru' : 101,   # Ruthenium
    'Rh' : 103,   # Rhodium
    'Pd' : 106,   # Palladium
    'Ag' : 108,   # Silver
    'Cd' : 112,   # Cadmium
    'In' : 115,   # Indium
    'Sn' : 119,   # Tin
    'Sb' : 122,   # Antimony
    'I'  : 127,   # Iodine
    'Xe' : 131,   # Xenon
    'Cs' : 133,   # Cesium
    'Ba' : 137,   # Barium
    'La' : 138,   # Lanthanum
    'Ce' : 140,   # Cerium
    'Pr' : 141,   # Praseodymium
    'Nd' : 144,   # Neodymium
    'Pm' : 145,   # Promethium
    'Sm' : 150,   # Samarium
    'Eu' : 152,   # Europium
    'Gd' : 157,   # Gadolinium
    'Tb' : 159,   # Terbium
    'Dy' : 162,   # Dysprosium
    'Ho' : 164,   # Holmium
    'Er' : 167,   # Erbium
    'Tm' : 169,   # Thulium
    'Yb' : 173,   # Ytterbium
    'Lu' : 175,   # Lutetium
    'Hf' : 178,   # Hafnium
    'Ta' : 180,   # Tantalum
    'W'  : 184,   # Tungsten
    'Re' : 186,   # Rhenium
    'Os' : 190,   # Osmium
    'Ir' : 192,   # Iridium
    'Pt' : 195,   # Platinum
    'Au' : 197,   # Gold
    'Hg' : 201,   # Mercury
    'Tl' : 204,   # Thallium
    'Pb' : 207,   # Lead
    'Bi' : 209,   # Bismuth
    'Po' : 209,   # Polonium
    'At' : 210,   # Astatine
    'Rn' : 222,   # Radon
}


num2mass = {
    '1'  : 1,     # Hydrogen
    '2'  : 4,     # Helium
    '3'  : 7,     # Lithium
    '4'  : 9,     # Beryllium
    '5'  : 11,    # Boron
    '6'  : 12,    # Carbon
    '7'  : 14,    # Nitrogen
    '8'  : 16,    # Oxygen
    '9'  : 18,    # Fluorine
    '10' : 20,    # Neon
    '11' : 23,    # Sodium
    '12' : 24,    # Magnesium
    '13' : 27,    # Aluminum
    '14' : 28,    # Silicon
    '15' : 31,    # Phosphorus
    '16' : 32,    # Sulfur
    '17' : 35,    # Chlorine
    '18' : 40,    # Argon
    '19' : 39,    # Potassium
    '20' : 40,    # Calcium
    '21' : 45,    # Scandium
    '22' : 47,    # Titanium
    '23' : 51,    # Vanadium
    '24' : 52,    # Chromium
    '25' : 55,    # Manganese
    '26' : 56,    # Iron
    '27' : 59,    # Nickel
    '28' : 63,    # Copper
    '29' : 65,    # Zinc
    '30' : 69,    # Gallium
    '31' : 73,    # Germanium
    '32' : 75,    # Arsenic
    '33' : 79,    # Selenium
    '34' : 80,    # Bromine
    '35' : 84,    # Krypton
    '36' : 85,    # Rubidium
    '37' : 87,    # Strontium
    '38' : 89,    # Yttrium
    '39' : 91,    # Zirconium
    '40' : 93,    # Niobium
    '41' : 96,    # Molybdenum
    '42' : 98,    # Technetium
    '43' : 101,   # Ruthenium
    '44' : 103,   # Rhodium
    '45' : 106,   # Palladium
    '46' : 108,   # Silver
    '47' : 112,   # Cadmium
    '48' : 115,   # Indium
    '49' : 119,   # Tin
    '50' : 122,   # Antimony
    '51' : 127,   # Iodine
    '52' : 131,   # Xenon
    '53' : 133,   # Cesium
    '54' : 137,   # Barium
    '55' : 138,   # Lanthanum
    '56' : 140,   # Cerium
    '57' : 141,   # Praseodymium
    '58' : 144,   # Neodymium
    '59' : 145,   # Promethium
    '60' : 150,   # Samarium
    '61' : 152,   # Europium
    '62' : 157,   # Gadolinium
    '63' : 159,   # Terbium
    '64' : 162,   # Dysprosium
    '65' : 164,   # Holmium
    '66' : 167,   # Erbium
    '67' : 169,   # Thulium
    '68' : 173,   # Ytterbium
    '69' : 175,   # Lutetium
    '70' : 178,   # Hafnium
    '71' : 180,   # Tantalum
    '72' : 184,   # Tungsten
    '73' : 186,   # Rhenium
    '74' : 190,   # Osmium
    '75' : 192,   # Iridium
    '76' : 195,   # Platinum
    '77' : 197,   # Gold
    '78' : 201,   # Mercury
    '79' : 204,   # Thallium
    '80' : 207,   # Lead
    '81' : 209,   # Bismuth
    '82' : 209,   # Polonium
    '83' : 210,   # Astatine
    '84' : 222,   # Radon
}


num2atom = {
    '1'  : 'H',    # Hydrogen
    '2'  : 'He',   # Helium
    '3'  : 'Li',   # Lithium
    '4'  : 'Be',   # Beryllium
    '5'  : 'B',    # Boron
    '6'  : 'C',    # Carbon
    '7'  : 'N',    # Nitrogen
    '8'  : 'O',    # Oxygen
    '9'  : 'F',    # Fluorine
    '10' : 'Ne',   # Neon
    '11' : 'Na',   # Sodium
    '12' : 'Mg',   # Magnesium
    '13' : 'Al',   # Aluminum
    '14' : 'Si',   # Silicon
    '15' : 'P',    # Phosphorus
    '16' : 'S',    # Sulfur
    '17' : 'Cl',   # Chlorine
    '18' : 'Ar',   # Argon
    '19' : 'K',    # Potassium
    '20' : 'Ca',   # Calcium
    '21' : 'Sc',   # Scandium
    '22' : 'Ti',   # Titanium
    '23' : 'V',    # Vanadium
    '24' : 'Cr',   # Chromium
    '25' : 'Mn',   # Manganese
    '26' : 'Fe',   # Iron
    '27' : 'Ni',   # Nickel
    '28' : 'Cu',   # Copper
    '29' : 'Zn',   # Zinc
    '30' : 'Ga',   # Gallium
    '31' : 'Ge',   # Germanium
    '32' : 'As',   # Arsenic
    '33' : 'Se',   # Selenium
    '34' : 'Br',   # Bromine
    '35' : 'Kr',   # Krypton
    '36' : 'Rb',   # Rubidium
    '37' : 'Sr',   # Strontium
    '38' : 'Y',    # Yttrium
    '39' : 'Zr',   # Zirconium
    '40' : 'Nb',   # Niobium
    '41' : 'Mo',   # Molybdenum
    '42' : 'Tc',   # Technetium
    '43' : 'Ru',   # Ruthenium
    '44' : 'Rh',   # Rhodium
    '45' : 'Pd',   # Palladium
    '46' : 'Ag',   # Silver
    '47' : 'Cd',   # Cadmium
    '48' : 'In',   # Indium
    '49' : 'Sn',   # Tin
    '50' : 'Sb',   # Antimony
    '51' : 'I',    # Iodine
    '52' : 'Xe',   # Xenon
    '53' : 'Cs',   # Cesium
    '54' : 'Ba',   # Barium
    '55' : 'La',   # Lanthanum
    '56' : 'Ce',   # Cerium
    '57' : 'Pr',   # Praseodymium
    '58' : 'Nd',   # Neodymium
    '59' : 'Pm',   # Promethium
    '60' : 'Sm',   # Samarium
    '61' : 'Eu',   # Europium
    '62' : 'Gd',   # Gadolinium
    '63' : 'Tb',   # Terbium
    '64' : 'Dy',   # Dysprosium
    '65' : 'Ho',   # Holmium
    '66' : 'Er',   # Erbium
    '67' : 'Tm',   # Thulium
    '68' : 'Yb',   # Ytterbium
    '69' : 'Lu',   # Lutetium
    '70' : 'Hf',   # Hafnium
    '71' : 'Ta',   # Tantalum
    '72' : 'W',    # Tungsten
    '73' : 'Re',   # Rhenium
    '74' : 'Os',   # Osmium
    '75' : 'Ir',   # Iridium
    '76' : 'Pt',   # Platinum
    '77' : 'Au',   # Gold
    '78' : 'Hg',   # Mercury
    '79' : 'Tl',   # Thallium
    '80' : 'Pb',   # Lead
    '81' : 'Bi',   # Bismuth
    '82' : 'Po',   # Polonium
    '83' : 'At',   # Astatine
    '84' : 'Rn',   # Radon
}



atomnames = []
atommasses = []
coordinates = []
frequencies = []
intensities = []
normalmodes = []

def calculate_raman_intensity(frequencies, activities, excitation_freq, temperature):
    """
    Convert Raman activity to intensity using theoretical formula.
    
    Parameters:
        frequencies (array): Vibrational frequencies (cm^-1).
        activities (array): Raman activities.
        excitation_freq (float): Excitation frequency (cm^-1).
        temperature (float): Temperature in Kelvin.
    
    Returns:
        array: Raman intensities.
    """
    exp_factor = np.exp(-h * c * frequencies / (k * temperature))
    B_i = 1 - exp_factor
    non_zero = frequencies != 0
    intensities = np.zeros_like(frequencies)
    intensities[non_zero] = C * ((excitation_freq - frequencies[non_zero])**4) * (1 / frequencies[non_zero]) * (1 / B_i[non_zero]) * activities[non_zero]
    return intensities

'''Detecting if data comes from ORCA or GAUSSIAN'''
with open(filename, 'r') as in_file:
    for line in in_file:
        if 'ORCA' in line:
            which = 'ORCA'
            break
        elif 'TURBOMOLE' in line:
            which = 'TURBOMOLE'
            break
        elif 'Entering Gaussian System' in line:
            which = 'Gaussian'
            break

'''Grabing data from ORCA output file'''
if which == 'ORCA':
    with open(filename, 'r') as in_file:
        write = False
        for line in in_file:
            if write == True and len(line.split()) > 1:
                if line.split()[0] in atom2mass.keys():
                    atomnames.append(line.split()[0])
                    coordinates.append([float(x) for x in line.split()[1:]])
                    atommasses.append(atom2mass[line.split()[0]])
                    atommasses.append(atom2mass[line.split()[0]])
                    atommasses.append(atom2mass[line.split()[0]])
            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                write = True
                continue
            if 'CARTESIAN COORDINATES (A.U.)' in line:
                break
    
    
    with open(filename, 'r') as in_file:
        write = False
        for line in in_file:
            if write == True and len(line.split()) > 1:
                if 'cm**-1' in line:
                    frequencies.append(float(line.split()[1]))
            if 'VIBRATIONAL FREQUENCIES' in line:
                write = True
                continue
            if 'NORMAL MODES' in line:
                break
    

    for x in range(3*len(atomnames)):
        normalmodes.append([])
    with open(filename, 'r') as in_file:
        write = False
        for line in in_file:
            if write == True and '.' in line:
                if line.split()[0] in [str(x) for x in range(0,3*len(atomnames))]:
                    idx = int(line.split()[0])
                    for x in line.split()[1:]:
                        normalmodes[idx].append(float(x))
            if 'NORMAL MODES' in line:
                write = True
                continue
            if 'IR SPECTRUM' in line:
                break
    normalmodes = np.array(normalmodes).T
            

    for i in range(0,6):
        intensities.append(float(0))
    with open(filename, 'r') as in_file:
        write = False
        for line in in_file:
            if write == True and len(line.split()) > 1:
                if ':' in line:
                    intensities.append(float(line.split()[2]))
            if 'RAMAN SPECTRUM' in line:
                write = True
                continue
            if 'THERMOCHEMISTRY' in line:
                break


    while len(intensities) != len(frequencies):
        intensities.insert(0, float(0))


    '''Calculating Raman intensities'''
    frequencies = np.array(frequencies)
    intensities = np.array(intensities)

    if args.intensities:
        intensities = calculate_raman_intensity(frequencies, intensities, wl, temp)
            

'''TURBOMOLE data extraction'''
if which =='TURBOMOLE':
    with open(filename, 'r') as in_file:
        write = False
        for line in in_file:
            if write == True and len(line.split()) > 1:
                if line.split()[3].capitalize() in atom2mass.keys():
                    atomnames.append(line.split()[3].capitalize())
                    coordinates.append([float(x)*0.529177249 for x in line.split()[0:3]])
                    atommasses.append(atom2mass[line.split()[3].capitalize()])
                    atommasses.append(atom2mass[line.split()[3].capitalize()])
                    atommasses.append(atom2mass[line.split()[3].capitalize()])
            if 'Atomic coordinate, charge and isotop information' in line:
                write = True
                continue
            if 'center of nuclear mass' in line:
                break

    modelines=[]

    with open(filename, 'r') as in_file:
        write = False
        for line in in_file:
            if write == True and len(line.split()) > 1:
                if line.split()[2] in ['x','y','z']:
                    modelines.append(line.split()[3:])
                if line.split()[0] in ['x','y','z']:
                    modelines.append(line.split()[1:])
            if 'NORMAL MODES and VIBRATIONAL FREQUENCIES' in line:
                write = True
                continue
            if '**************************************************************' in line:
                break


    for x in range(3*len(atomnames)):
        normalmodes.append([])
    for idx, line in enumerate(modelines):
        for x in line:
            normalmodes[idx % (3*len(atomnames))].append(float(x))
    
    normalmodes = np.array(normalmodes).T

    if not args.control:  
        print('Error: TM vibrational analysis requires a control file to be specified.')
        sys.exit(1)
    else:
        with open(args.control, 'r') as in_file:
            write = False
            for line in in_file:
                if write == True and len(line.split()) == 8 and line.split()[0].isdigit():
                    intensities.append(float(line.split()[6].replace('D', 'e')) + float(line.split()[7].replace('D', 'e')))
                    frequencies.append(float(line.split()[2]))
                if 'raman spectrum' in line:
                    write = True
                    continue
                if 'end' == line:
                    break


    if len(frequencies) > len(normalmodes)-5:
        print('Error: Number of frequencies and normal modes does not match.\nCheck if the control file and the output file are consistent.')
        sys.exit(1)

    while len(frequencies) != len(normalmodes):
        frequencies.insert(0, 0)
        intensities.insert(0, 0)

    frequencies = np.array(frequencies)
    intensities = np.array(intensities)

    

'''Grabing data from GAUSSIAN output file'''
if which == 'Gaussian':
    with open(filename, 'r') as in_file:
        write = False
        for line in in_file:
            if write == True and len(line.split()) > 1:
                if line.split()[0] in [str(x) for x in range(1,100)]:
                    atomnames.append(num2atom[line.split()[1]])
                    coordinates.append([float(x) for x in line.split()[3:]])
                    atommasses.append(num2mass[line.split()[1]])
                    atommasses.append(num2mass[line.split()[1]])
                    atommasses.append(num2mass[line.split()[1]])
            if 'Input orientation:' in line:
                write = True
                continue
            if 'Distance matrix (angstroms)' in line:
                break
            
    with open(filename, 'r') as in_file:
        for line in in_file:
            if 'Frequencies --' in line:
                frequencies.append(float(line.split()[-3]))
                frequencies.append(float(line.split()[-2]))
                frequencies.append(float(line.split()[-1]))
            if '- Thermochemistry -' in line:
                break

    with open(filename, 'r') as in_file:
        for line in in_file:
            if 'Raman Activ --' in line:
                intensities.append(float(line.split()[-3]))
                intensities.append(float(line.split()[-2]))
                intensities.append(float(line.split()[-1]))
            if '- Thermochemistry -' in line:
                break    
        

    with open(filename, 'r') as in_file:
        write = False
        temp1 = []
        temp2 = []
        temp3 = []
        for line in in_file:
            if write == True:
                if len(line.split()) > 3:
                    if line.split()[0] in [str(x) for x in range(1,100)]:
                        temp1.append(float(line.split()[2]))
                        temp1.append(float(line.split()[3]))
                        temp1.append(float(line.split()[4]))
                        temp2.append(float(line.split()[5]))
                        temp2.append(float(line.split()[6]))
                        temp2.append(float(line.split()[7]))
                        temp3.append(float(line.split()[8]))
                        temp3.append(float(line.split()[9]))
                        temp3.append(float(line.split()[10]))
            if 'Frequencies --' in line:
                if temp1 != []:
                    normalmodes.append(temp1)
                    normalmodes.append(temp2)
                    normalmodes.append(temp3)
                temp1 = []
                temp2 = []
                temp3 = []
                write = True
            if '- Thermochemistry -' in line:
                break
    normalmodes = np.array(normalmodes)


if 3*len(atomnames) != len(normalmodes):
    print('Error: Number of atoms and normal modes does not match.\nPossibly the file contains unidentified atom types.')
    sys.exit(1)
    

'''Normalizing intensities'''
intensities = 10000*(intensities - np.min(intensities))/(np.max(intensities) - np.min(intensities))


'''Reformatting data'''
mwdispls = []
mwdispls3D = []
for nm in normalmodes:
    temp_displs = np.multiply(np.sqrt(atommasses[0:len(nm)]),nm)
    mwdispls.append(temp_displs)
    mwdispls3D.append(np.reshape(temp_displs,(-1,3)))


'''Normalizing mass-weighted displacements'''
coordcontributions = []
for nm in mwdispls:
    total = np.sum(nm**2)
    if total > 0:
        normnm = np.divide(nm**2,total)
    else:
        normnm = np.zeros(3*len(atomnames))
    coordcontributions.append(normnm)

'''Calculating contributions per atom'''
atomcontributions =[]
for nm in coordcontributions:
    temp = np.reshape(nm,(-1,3))
    temp = np.sum(temp,axis=1)
    atomcontributions.append(temp)

'''Output contributions per atom'''
with open('Z_OUT_' + filename, 'w') as out_file:
    string_list = ['idx', 'v']
    out_file.write(''.join(x.rjust(8) for x in string_list))
    out_file.write(''.join(x.rjust(6) for x in atomnames) + '\n')
    out_file.write(''.join(x.rjust(8) for x in [' ',' ']))
    out_file.write(''.join(str(x).rjust(6) for x in range(1,len(atomnames)+1)) + '\n')
    for idx, freq, nm in zip(range(len(frequencies)),frequencies, atomcontributions):
        if min_max[0] <= float(freq)*scaling <= min_max[1]:
            pass
        else: 
            continue
        out_file.write(str(idx).rjust(8))
        out_file.write(str(round(scaling*freq,1)).rjust(8))
        out_file.write(''.join(str(round(100*x,1)).rjust(6) for x in nm))
        out_file.write('\n')

df = pd.read_csv('Z_OUT_' + filename, delim_whitespace=True)
df.to_csv('Z_OUT_' + filename + '.csv', index=False)
        
'''Estimating contributions of CO stretch, NH bend, and NH2 bend
   - these are estimates based on total displacements on the relevant atoms'''
modes_CO = []
modes_NH = []
modes_NH2 = []
tmp_atoms = []


if args.indices_file:
    try:
        if args.indices_file not in os.listdir():
            raise FileNotFoundError('Error: Indices file not found.')
        with open(args.indices_file, 'r') as file:
            idx_file = file.readlines()

            ## Remove comments from the file
            for ln_num, line in enumerate(idx_file):
                if '#' in line:
                    comm_idx = None
                    for idx, words in enumerate(line):
                        if '#' in words:
                            comm_idx = idx
                            break
                    idx_file[ln_num] = line[:comm_idx]

            #allows naming the groups for easier output reading
            alias = []
            for ln_num, line in enumerate(idx_file):
                patern = r'([\s\d]+)([\'\"])(\w*)([\'\"])'
                if re.search(patern, line):
                    alias.append((re.search(patern, line).group(3)).rjust(12))
                    idx_file[ln_num] = re.search(patern, line).group(1)
                else:
                    if len(line.split()) > 4:
                        alias.append('...'.rjust(12))
                    else:
                        alias.append(''.join([atomnames[int(atom)-1] for atom in line.split()]).rjust(12))

            idx_file = [list(map(int, x.split())) for x in idx_file]
            idx_file = list(map(lambda x: [y-1 for y in x], idx_file))
        pop_line = []
        for idx, line in enumerate(idx_file):
#            if len(line) == 1:
#                print(f'Error: Not enough indices in the line: {line}. Specify at least two indices for a group.')
#            elif len(line) == 0:
#                idx_file.pop(idx)
#            for i in line:
#                if i not in range(len(atomnames)):
#                    print(f'Error: Index {i+1} out of range.')
#                    sys.exit(1)
#                i_xyz = coordinates[i]
#                for j in line:
#                    if j not in range(len(atomnames)):
#                        print(f'Error: Index {j+1} out of range.')
#                        sys.exit(1)
#                    if i == j:
#                        continue
#                    j_xyz = coordinates[j]
#                    ij_vec = np.subtract(j_xyz, i_xyz)
#                    if np.sqrt(np.dot(ij_vec,ij_vec)) >= 1.5:
#                        dist = np.sqrt(np.dot(ij_vec,ij_vec))
#                        print(f'Warning: Atoms {i+1} and {j+1} are too far apart ({dist}) to be considered bonded.')
#                        continue
            con = [[] for x in range(len(line))]
            if len(line) == 0:
                pop_line.append(idx)
                continue
            for i_idx, i in enumerate(line):
                if i not in range(len(atomnames)):
                    print(f'Error: Index {i+1} out of range.')
                    sys.exit(1)
                i_xyz = coordinates[i]
                for j in line:
                    if j not in range(len(atomnames)):
                        print(f'Error: Index {j+1} out of range.')
                        sys.exit(1)
                    if i == j:
                        con[i_idx].append(10)
                        continue
                    j_xyz = coordinates[j]
                    ij_vec = np.subtract(j_xyz, i_xyz)
                    con[i_idx].append(np.sqrt(np.dot(ij_vec,ij_vec)))
            con = np.array(con)
            for idx, col in enumerate(np.all(con > 1.5, axis=0)):
                if col:
                    min_dist = np.min(con[idx])
                    print(f'!!!Warning: Atom {line[idx]+1} might be too far apart from the rest of the specified functional group ({min_dist:.4f}A from the closest atom).')
        pop_line.sort(reverse=True)
        for line in pop_line:
            idx_file.pop(line)
            alias.pop(line)
            
    except FileNotFoundError as fnfe:
        print(fnfe)
        sys.exit(1)
    except IOError as ie:
        print(ie)
        sys.exit(1)
    

    
elif args.indices:
    for idx in args.indices:
        if idx not in range(len(atomnames)+1):
            print(f'Error: Index {idx} out of range.')
            sys.exit(1)
        tmp_atoms.append(idx - 1)

    print(f'\nUsing user-defined indices for functional group analysis: {args.indices}\n')
else:
    tmp_atoms = range(len(atomnames))
    print('\nNo user-defined indices for functional group analysis. Attempting to guess the functional groups based on atom types and proximity.\n')

if args.indices_file:
    groups_contributions = []
    for displs in mwdispls3D:
        total = np.sum(displs**2)
        temp_group = []
        for group in idx_file:
            temp_atom = []
            for atom in group:
                temp_atom.append(displs[atom])
            if total == 0 or total == np.nan:
                temp_group.append(0)
            else:
                temp_group.append(np.sum(np.array(temp_atom)**2)/total)
        groups_contributions.append(temp_group)
      

    '''Output digested contributions and input for plot_Raman_specs.py'''
    with open('Z_MODES_' + filename, 'w') as out_file:
        numbering = []
        for line in idx_file:
            if len(line) > 4:
                numbering.append(f'{line[0]+1},...,{line[-1]+1}'.rjust(12))
            else:
                numbering.append(','.join([str(num+1) for num in line]).rjust(12))
        string_list = ['idx', 'v']
        if args.intensities:
            string_list.append('Int')
        out_file.write(''.join(x.rjust(8) for x in string_list))
        out_file.write(''.join(name for name in alias) + '\n')
        out_file.write(''.join(x.rjust(8) for x in [' ',' ']))
        if args.intensities:
            out_file.write(' '.rjust(9))
        for num in numbering:
            out_file.write(num)
        out_file.write('\n')
        iterable = zip(range(len(frequencies)),frequencies, intensities, groups_contributions)
        if args.sort == None:
            pass
        elif sort_idx+1 <= len(groups_contributions[0]):
            iterable = sorted(iterable, key=lambda x: x[3][sort_idx], reverse=True)
        elif sort_idx+1 > len(groups_contributions[0]):
            print('\nThe group index for sorting is out of range!\n'.upper())
        for idx, freq, intens, group in iterable:
            if min_max[0] <= float(freq)*scaling <= min_max[1]:
                pass
            else: 
                continue
            
            out_file.write(str(idx+1).rjust(8))
            out_file.write(str(round(scaling*freq,1)).rjust(8))
            if args.intensities:
                out_file.write(str(round(intens,2)).rjust(8))
            out_file.write(''.join(str(round(100*x,1)).rjust(12) for x in group))
            out_file.write('\n')

    df = pd.read_csv('Z_MODES_' + filename, delim_whitespace=True)
    df.to_csv('Z_MODES_' + filename + '.csv', index=False)


elif args.indices:
    groups_contributions = []
    for displs in mwdispls3D:
        total = np.sum(displs**2)
        temp_group = []
        temp_atom = []
        for atom in tmp_atoms:
            temp_atom.append(displs[atom])
        if total == 0 or total == np.nan:
            temp_group.append(0)
        else:
            temp_group.append(np.sum(np.array(temp_atom)**2)/total)
        groups_contributions.append(temp_group)


    '''Output digested contributions and input for plot_Raman_specs.py'''
    with open('Z_MODES_' + filename, 'w') as out_file:
        string_list = ['idx', 'v']
        if args.intensities:
            string_list.append('Int')
        out_file.write(''.join(x.rjust(8) for x in string_list))
        out_file.write(''.join([atomnames[atom] for atom in tmp_atoms]).rjust(10) + '\n')
        out_file.write(''.join(x.rjust(8) for x in [' ',' ']))
        if args.intensities:
            out_file.write(' '.rjust(8))
        out_file.write(','.join(str(x+1) for x in tmp_atoms).rjust(10) + '\n')
        for idx, freq in zip(range(len(frequencies)), frequencies):
            if min_max[0] <= float(freq)*scaling <= min_max[1]:
                pass
            else: 
                continue
            
            out_file.write(str(idx+1).rjust(8))
            out_file.write(str(round(scaling*freq,1)).rjust(8))
            if args.intensities:
                out_file.write(str(round(intensities[idx],2)).rjust(8))
            out_file.write(''.join(str(round(100*x,1)).rjust(10) for x in groups_contributions[idx]))
            out_file.write('\n')

    df = pd.read_csv('Z_MODES_' + filename, delim_whitespace=True)
    df.to_csv('Z_MODES_' + filename + '.csv', index=False)

else:
    for idx in tmp_atoms:
        temp_CO = []
        if atomnames[idx] == 'O':
            O_idx = idx
            O_xyz = coordinates[idx]
            temp_CO.append(idx)
            for jdx in tmp_atoms:
                if jdx != idx and atomnames[jdx] == 'C':
                    C_idx = jdx
                    C_xyz = coordinates[jdx]
                    CO_vec = np.subtract(C_xyz, O_xyz)
                    if np.sqrt(np.dot(CO_vec,CO_vec)) <= 1.5:
                        temp_CO.append(jdx)
            if len(temp_CO) == 2:
                modes_CO.append(temp_CO)


        temp_NHx = []
        if atomnames[idx] == 'N':
            N_idx = idx
            N_xyz = coordinates[idx]
            temp_NHx.append(idx)
            for jdx in tmp_atoms:
                if jdx != idx and atomnames[jdx] == 'H':
                    H_idx = jdx
                    H_xyz = coordinates[jdx]
                    NH_vec = np.subtract(N_xyz, H_xyz)
                    if np.sqrt(np.dot(NH_vec,NH_vec)) <= 1.3:
                        temp_NHx.append(jdx)
            if len(temp_NHx) == 3:
                for jdx in tmp_atoms:
                    if jdx != idx and atomnames[jdx] == 'C':
                        C_idx = jdx
                        C_xyz = coordinates[jdx]
                        NC_vec = np.subtract(N_xyz, C_xyz)
                        if np.sqrt(np.dot(NC_vec,NC_vec)) <= 1.5:
                            temp_NHx.append(jdx)
                if len(temp_NHx) == 4:
                    modes_NH2.append(temp_NHx)            
            elif len(temp_NHx) == 2:
                for jdx in tmp_atoms:
                    if jdx != idx and atomnames[jdx] == 'C':
                        C_idx = jdx
                        C_xyz = coordinates[jdx]
                        NC_vec = np.subtract(N_xyz, C_xyz)
                        if np.sqrt(np.dot(NC_vec,NC_vec)) <= 1.5:
                            temp_NHx.append(jdx)
                if len(temp_NHx) == 4:
                    modes_NH.append(temp_NHx)                    

    CO = []
    for displs in mwdispls3D:
        total = np.sum(displs**2)
        temp_CO = []
        for mode in modes_CO:      
            C_displ = displs[mode[1]]
            O_displ = displs[mode[0]]
            if total == 0:
                contrib = 0
            else:
                contrib = (np.sum(C_displ**2) + np.sum(O_displ**2))/total
            temp_CO.append(contrib)
        CO.append(temp_CO)

    NH2 = []
    for displs in mwdispls3D:
        total = np.sum(displs**2)
        temp_NH2 = []
        for mode in modes_NH2:
            N_displ = displs[mode[0]]
            H1_displ = displs[mode[1]]
            H2_displ = displs[mode[2]]
            if total == 0:
                contrib = 0
            else:
                contrib = (np.sum(N_displ**2) + np.sum(H1_displ**2) + np.sum(H2_displ**2))/total 
            temp_NH2.append(contrib)
        NH2.append(temp_NH2)

    NH = []
    for displs in mwdispls3D:
        total = np.sum(displs**2)
        temp_NH = []
        for mode in modes_NH:
            N_displ = displs[mode[0]]
            H_displ = displs[mode[1]]
            if total == 0:
                contrib = 0
            else:
                contrib = (np.sum(N_displ**2) + np.sum(H_displ**2))/total 
            temp_NH.append(contrib)
        NH.append(temp_NH)


    '''Output digested contributions and input for plot_Raman_specs.py'''
    with open('Z_MODES_' + filename, 'w') as out_file:
        string_list = ['idx', 'v']
        if args.intensities:
            string_list.append('Int')
        out_file.write(''.join(x.rjust(8) for x in string_list))
        out_file.write(''.join(('CO').rjust(10) for x in modes_CO))
        out_file.write(''.join(('NH2').rjust(10) for x in modes_NH2))
        out_file.write(''.join(('NH').rjust(10) for x in modes_NH) + '\n')
        if args.intensities:
            out_file.write(' '.rjust(9))
        out_file.write(''.join(x.rjust(8) for x in [' ',' ']))
        out_file.write(''.join((','.join(str(y+1) for y in x)).rjust(10) for x in modes_CO))
        out_file.write(''.join((','.join(str(y+1) for y in x[0:3])).rjust(10) for x in modes_NH2))
        out_file.write(''.join((','.join(str(y+1) for y in x[0:2])).rjust(10) for x in modes_NH) + '\n')

        contr = [x+y+z for x, y, z in zip(CO, NH2, NH)]
        iterable = zip(range(len(frequencies)),frequencies, intensities, contr)
        if args.sort == None:
            pass
        elif sort_idx+1 <= len(contr[0]):
            iterable = sorted(iterable, key=lambda x: x[3][sort_idx], reverse=True)
        elif sort_idx+1 > len(contr[0]):
            print('\nThe group index for sorting is out of range!\n'.upper())


        for idx, freq, intens, group in iterable:
            if min_max[0] <= float(freq)*scaling <= min_max[1]:
                pass
            else: 
                continue
            
            out_file.write(str(idx+1).rjust(8))
            out_file.write(str(round(scaling*freq,1)).rjust(8))
            if args.intensities:
                out_file.write(str(round(intens,2)).rjust(8))
            out_file.write(''.join(str(round(100*x,1)).rjust(10) for x in group))
            out_file.write('\n')

    df = pd.read_csv('Z_MODES_' + filename, delim_whitespace=True)
    df.to_csv('Z_MODES_' + filename + '.csv', index=False)


if args.print:
    print('\n\t*** Group contributions ***\n'.upper())
    with open('Z_MODES_' + filename, 'r') as out_file:
        for line in out_file:
            print(line)
elif args.print_all:
    print('\n\t*** Group contributions ***\n'.upper())
    with open('Z_MODES_' + filename, 'r') as out_file:
        for line in out_file:
            print(line)
    print('\n\n')
    print('\n\t*** Atom contributions ***\n'.upper())
    with open('Z_OUT_' + filename, 'r') as out_file:
        for line in out_file:
            print(line)


'''Graveyard'''
# CO_str = []
# for displs in mwdispls3D:
#     total = np.sum(displs**2)
#     temp_CO = []
#     for mode in modes_CO:
#         C_xyz = coordinates[mode[1]]
#         O_xyz = coordinates[mode[0]]
#         CO_vec = np.subtract(C_xyz,O_xyz)
#         CO_vec = CO_vec / np.linalg.norm(CO_vec)
        
#         C_displ = np.dot(displs[mode[1]],CO_vec)
#         O_displ = np.dot(displs[mode[0]],CO_vec)
#         trans = (C_displ + O_displ)**2
#         stretch = (C_displ - O_displ)**2
#         #Optional weighting - how much of it is stretching:
#         correct = stretch/(trans+stretch)
#         contrib = (C_displ**2 + O_displ**2)/total * correct
#         temp_CO.append(round(100*contrib,1))
#     CO_str.append(temp_CO)
        
# NH2_sciss = []
# NH2_wagg = []
# for displs in mwdispls3D:
#     total = np.sum(displs**2)
#     temp_NH2s = []
#     temp_NH2w = []
#     for mode in modes_NH2:
#         N_xyz = coordinates[mode[0]]
#         H1_xyz = coordinates[mode[1]]
#         H2_xyz = coordinates[mode[2]]
#         C_xyz = coordinates[mode[3]]
        
#         NH1_vec = np.subtract(H1_xyz,N_xyz)
#         NH1_vec = NH1_vec / np.linalg.norm(NH1_vec)
#         NH2_vec = np.subtract(H2_xyz,N_xyz)
#         NH2_vec = NH2_vec / np.linalg.norm(NH2_vec)
#         bisec_vec = np.add(np.subtract(H1_xyz,N_xyz), np.subtract(H2_xyz,N_xyz))
#         bisec_vec = bisec_vec / np.linalg.norm(bisec_vec)
#         HH_vec = np.subtract(H1_xyz,H2_xyz)
#         HH_vec = HH_vec / np.linalg.norm(HH_vec)
        
#         NH1_vec = np.cross(NH1_vec,bisec_vec)
#         NH1_vec = NH1_vec / np.linalg.norm(NH1_vec)
#         NH1_vec = np.cross(NH1_vec,bisec_vec)
#         NH1_vec = NH1_vec / np.linalg.norm(NH1_vec)
        
#         NH2_vec = np.cross(NH2_vec,bisec_vec)
#         NH2_vec = NH2_vec / np.linalg.norm(NH2_vec)
#         NH2_vec = np.cross(NH2_vec,bisec_vec)
#         NH2_vec = NH2_vec / np.linalg.norm(NH2_vec)
        
#         #scissoring:
#         N_z = np.dot(displs[mode[0]],bisec_vec)
#         NH1 = np.dot(displs[mode[1]],NH1_vec)
#         NH2 = np.dot(displs[mode[2]],NH2_vec)  
#         C_z = np.dot(displs[mode[3]],bisec_vec)
#         H1_x = np.dot(displs[mode[1]],HH_vec)
#         H2_x = np.dot(displs[mode[2]],HH_vec)
#         frac1 = (N_z - C_z)**2 / ((N_z - C_z)**2 + (N_z + C_z)**2)
#         frac2 = (H1_x - H2_x)**2 / ((H1_x - H2_x)**2 + (H1_x + H2_x)**2)
#         contrib = (N_z**2 + C_z**2)/total * frac1 + (NH1**2 + NH2**2)/total * frac2
#         temp_NH2s.append(round(100*contrib,1))
        
#         #wagging
#         N_x = np.dot(displs[mode[0]],HH_vec)
#         frac = (H1_x + H2_x)**2 / ((H1_x - H2_x)**2 + (H1_x + H2_x)**2)
#         contrib = (N_x**2)/total + (NH1**2 + NH2**2)/total * frac #* (asym/(sym+asym))
#         temp_NH2w.append(round(100*contrib,1))
#     NH2_sciss.append(temp_NH2s)
#     NH2_wagg.append(temp_NH2w)

# for nm in mwdispls:
#     total = np.sum(nm**2)
#     if total > 0:
#         normnm = np.divide(nm**2,total)
#     else:
#         normnm = np.zeros(3*len(atomnames))
#     coordcontributions.append(normnm)
        
        
        
        
        
        
        
        