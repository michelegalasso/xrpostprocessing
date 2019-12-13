"""
@file:      read_convex_hull.py
@author:    Michele Galasso
@brief:     Reads the necessary information from extended_convex_hull and extended_convex_hull_POSCARS.
"""

from pymatgen.core.structure import Structure

from .iterator_poscar_file import iterator_poscar_file


def read_convex_hull(extended_convex_hull, extended_convex_hull_POSCARS, remove_hydrogens=False):
    iterator = iterator_poscar_file(extended_convex_hull_POSCARS)
    data = []
    with open(extended_convex_hull, 'r') as f:
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
                if ID == values[0]:
                    structure = Structure.from_str(string, fmt='poscar')
                    pmg_composition = structure.composition

                    if remove_hydrogens:
                        structure.remove_species(['H'])

                    if len(structure.frac_coords) != 0:
                        enthalpy = float(values[3 + len(composition)])
                        fitness = float(values[5 + len(composition)])
                        data.append((structure, ID, enthalpy, fitness, pmg_composition))
                else:
                    raise IOError('Structures in {} do not match data in {}'.format(extended_convex_hull_POSCARS,
                                                                                    extended_convex_hull))
    return data
