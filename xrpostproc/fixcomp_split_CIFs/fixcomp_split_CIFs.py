"""
@file:      fixcomp_split_CIFs.py
@author:    Michele Galasso
@brief:     Script for splitting fixcomp USPEX results into multiple CIF files.
"""

import os
import sys

from xrpostproc.common.read_structures import read_structures
from xrpostproc.common.create_cif import create_cif


if len(sys.argv) != 4:
    print('ERROR: wrong number of arguments.')
    print('Usage: python fixcomp_split_CIFs.py [Individuals] [gatheredPOSCARS] [pressure]')
    print('Example: python fixcomp_split_CIFs.py Individuals gatheredPOSCARS 50GPa')
    sys.exit()

# read arguments
individuals = sys.argv[1]
gatheredPOSCARS = sys.argv[2]
pressure = sys.argv[3]

# read files
data = read_structures(individuals, gatheredPOSCARS, fixcomp=True)

# sort according to fitness
data.sort(key=lambda tup: tup[3])

os.mkdir('results')
for i, (structure, ID, enthalpy, fitness, pmg_composition) in enumerate(data):
    create_cif(i, structure, ID, enthalpy, fitness, pmg_composition, pressure)

print('Done.')
