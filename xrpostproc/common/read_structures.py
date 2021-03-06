"""
@file:      read_structures.py
@author:    Michele Galasso
@brief:     Reads the necessary information from individuals-like and gatheredPOSCARS-like files.
"""

from pymatgen.core.structure import Structure

from .iterator_poscar_file import iterator_poscar_file


def read_structures(individuals, gatheredPOSCARS, remove_hydrogens=False, fixcomp=False):
    iterator = iterator_poscar_file(gatheredPOSCARS)
    data = []
    with open(individuals, 'r') as f:
        while True:
            current_line = f.readline()
            if current_line == '':
                break
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
