"""
@file:      create_cif.py
@author:    Michele Galasso
@brief:     Function for creating cif files.
"""

import os
import spglib


def create_cif(i, structure, ID, enthalpy, fitness, pmg_composition, pressure):
    dtset = spglib.get_symmetry_dataset((structure.lattice.matrix, structure.frac_coords, structure.atomic_numbers),
                                        symprec=0.2)
    filename = '{}_EA{}_{}_{}_{}_{}_spg{}.cif'.format(i + 1, ID, fitness, enthalpy, pmg_composition.iupac_formula,
                                                      pressure, dtset['number'])
    structure.to(filename=os.path.join('results', filename), symprec=0.2)
