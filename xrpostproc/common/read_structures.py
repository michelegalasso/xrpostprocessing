"""
@file:      read_structures.py
@author:    Michele Galasso
@brief:     Reads the necessary information from individuals-like and gatheredPOSCARS-like files.
"""

from pymatgen.core.structure import Structure

from .iterator_poscar_file import iterator_poscar_file


def read_structures(gatheredPOSCARS, remove_hydrogens=False):
    # initialization
    data = []

    # read data
    for string in iterator_poscar_file(gatheredPOSCARS):
        ID = string.split()[0].strip('EA')

        structure = Structure.from_str(string, fmt='poscar')
        pmg_composition = structure.composition

        if remove_hydrogens:
            structure.remove_species(['H'])

        if len(structure.frac_coords) != 0:
            data.append((structure, ID, pmg_composition))
    return data
