"""
@file:      read_structures.py
@author:    Michele Galasso
@brief:     Reads the necessary information from individuals-like and gatheredPOSCARS-like files.
"""

from pymatgen.core.structure import Structure

from .iterator_poscar_file import iterator_poscar_file


def read_structures(individuals, gatheredPOSCARS, remove_hydrogens=False, fixcomp=False):
    iterator = iterator_poscar_file(gatheredPOSCARS)
    uspexpy = False                 # suppose that we are reading MATLAB USPEX
    data = []
    with open(individuals, 'r') as f:
        while True:
            current_line = f.readline()
            if current_line == '':
                break
            if '+----' in current_line:
                uspexpy = True      # we are reading pythonic USPEX
                continue

            if uspexpy:
                values = current_line.split('|')
                if 'ID' in values[1]:
                    headers = []
                    for value in values:
                        headers.append(value.strip())
                else:
                    string = next(iterator)
                    ID = string.split()[0].strip('EA')

                    if ID == values[1].strip():
                        structure = Structure.from_str(string, fmt='poscar')
                        pmg_composition = structure.composition

                        if remove_hydrogens:
                            structure.remove_species(['H'])

                        if len(structure.frac_coords) != 0:
                            enthalpy_position = headers.index('Enthalpy (eV)')

                            if fixcomp:
                                enthalpy = float(values[enthalpy_position])
                                fitness = float(values[enthalpy_position]) / structure.num_sites
                            else:
                                fitness_position = headers.index('Enthalpy above CH (eV/block)')

                                enthalpy = float(values[enthalpy_position])
                                fitness = float(values[fitness_position])

                            data.append((structure, ID, enthalpy, fitness, pmg_composition))
                    else:
                        raise IOError('Structures in {} do not match data in {}'.format(gatheredPOSCARS, individuals))
            else:
                values = current_line.split()
                if values[0].isdigit():
                    composition_string = current_line[current_line.find('[') + 1:current_line.find(']')]
                    composition = [int(x) for x in composition_string.split()]

                    string = next(iterator)
                    ID = string.split()[0].strip('EA')

                    if fixcomp:
                        ID_position = 1
                    else:
                        ID_position = 0

                    if ID == values[ID_position]:
                        structure = Structure.from_str(string, fmt='poscar')
                        pmg_composition = structure.composition

                        if remove_hydrogens:
                            structure.remove_species(['H'])

                        if len(structure.frac_coords) != 0:
                            if fixcomp:
                                enthalpy = float(values[5 + len(composition)])
                                fitness = float(values[5 + len(composition)]) / structure.num_sites
                            else:
                                enthalpy = float(values[3 + len(composition)])
                                fitness = float(values[5 + len(composition)])

                            data.append((structure, ID, enthalpy, fitness, pmg_composition))
                    else:
                        raise IOError('Structures in {} do not match data in {}'.format(gatheredPOSCARS, individuals))
    return data
