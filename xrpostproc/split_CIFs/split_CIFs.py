"""
@file:      split_CIFs.py
@author:    Michele Galasso
@brief:     Script for splitting USPEX results into multiple CIF files.
"""

import os
import sys

from xrpostproc.common.read_convex_hull import read_convex_hull
from xrpostproc.common.create_cif import create_cif


if len(sys.argv) != 4:
    print('ERROR: wrong number of arguments.')
    print('Usage: python split_CIFs.py [extended_convex_hull] [extended_convex_hull_POSCARS] [pressure]')
    print('Example: python split_CIFs.py extended_convex_hull extended_convex_hull_POSCARS 50GPa')
    sys.exit()

# read arguments
extended_convex_hull = sys.argv[1]
extended_convex_hull_POSCARS = sys.argv[2]
pressure = sys.argv[3]

# read files
data = read_convex_hull(extended_convex_hull, extended_convex_hull_POSCARS)

# sort according to enthalpy
data.sort(key=lambda tup: tup[2])
systems_to_print = []
compositions = {}

os.mkdir('results')
for structure, ID, enthalpy, fitness in data:
    if structure.composition.reduced_formula in compositions.keys():
        if compositions[structure.composition.reduced_formula] < 5:
            compositions[structure.composition.reduced_formula] += 1
            systems_to_print.append((structure, ID, enthalpy, fitness))
    else:
        compositions[structure.composition.reduced_formula] = 1
        systems_to_print.append((structure, ID, enthalpy, fitness))

# sort according to fitness
systems_to_print.sort(key=lambda tup: tup[3])

for i, (structure, ID, enthalpy, fitness) in enumerate(systems_to_print):
    create_cif(i, structure, ID, enthalpy, fitness, pressure)

print('Done.')
