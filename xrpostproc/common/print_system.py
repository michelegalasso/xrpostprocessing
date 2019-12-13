"""
@file:      print_system.py
@author:    Michele Galasso
@brief:     Function for creating cif files.
"""

import os
import spglib


def print_system(i, structure, enthalpy, fitness, pressure):
    dtset = spglib.get_symmetry_dataset((structure.lattice.matrix, structure.frac_coords, structure.atomic_numbers),
                                        symprec=0.2)
    filename = '{}_{}_{}_{}_{}_{}.cif'.format(i + 1, fitness, enthalpy, structure.composition.reduced_formula,
                                              pressure, dtset['number'])
    structure.to(filename=os.path.join('results', filename), symprec=0.2)
