"""
@file:      fixcomp_sublattice_split_CIFs.py
@author:    Michele Galasso
@brief:     Script for splitting varcomp USPEX results into multiple CIF files without hydrogens.
"""

import os
import sys

from xrpostproc.common.read_structures import read_structures
from xrpostproc.common.create_cif import create_cif


if len(sys.argv) not in [4, 5]:
    print('ERROR: wrong number of arguments.')
    print('Usage: python fixcomp_sublattice_split_CIFs.py [extended_convex_hull] [extended_convex_hull_POSCARS] [pressure]')
    print('Example: python fixcomp_sublattice_split_CIFs.py extended_convex_hull extended_convex_hull_POSCARS 50GPa')
    sys.exit()

# read arguments
extended_convex_hull = sys.argv[1]
extended_convex_hull_POSCARS = sys.argv[2]
pressure = sys.argv[3]

if len(sys.argv) == 4:
    n_returned_structs = 10
else:
    n_returned_structs = int(sys.argv[4])

# read files
data = read_structures(extended_convex_hull, extended_convex_hull_POSCARS, remove_hydrogens=True, fixcomp=True)

# sort according to enthalpy
data.sort(key=lambda tup: tup[2])
systems_to_print = []
compositions = {}

os.mkdir('results')
for structure, ID, enthalpy, fitness, pmg_composition in data:
    if pmg_composition.reduced_formula in compositions.keys():
        if compositions[pmg_composition.reduced_formula] < n_returned_structs:
            compositions[pmg_composition.reduced_formula] += 1
            systems_to_print.append((structure, ID, enthalpy, fitness, pmg_composition))
    else:
        compositions[pmg_composition.reduced_formula] = 1
        systems_to_print.append((structure, ID, enthalpy, fitness, pmg_composition))

# sort according to fitness
systems_to_print.sort(key=lambda tup: tup[3])

for i, (structure, ID, enthalpy, fitness, pmg_composition) in enumerate(systems_to_print):
    create_cif(i, structure, ID, enthalpy, fitness, pmg_composition, pressure)

print('Done.')
