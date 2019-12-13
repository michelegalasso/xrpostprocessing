import os
import sys
import spglib

from pymatgen.core.structure import Structure

from xrpostproc.common.read_convex_hull import read_convex_hull
from xrpostproc.common.create_cif import create_cif


def print_system(i, structure, ID, enthalpy, fitness, full_composition, pressure):
    dtset = spglib.get_symmetry_dataset((structure.lattice.matrix, structure.frac_coords, structure.atomic_numbers),
                                        symprec=0.2)
    filename = '{}_{}_{}_{}_{}_{}.cif'.format(i + 1, fitness, enthalpy, full_composition, pressure, dtset['number'])
    structure.to(filename=os.path.join('results', filename), symprec=0.2)


if len(sys.argv) != 4:
    print('ERROR: wrong number of arguments.')
    print('Usage: python sublattice_split_CIFs.py [extended_convex_hull] [extended_convex_hull_POSCARS] [pressure]')
    print('Example: python sublattice_split_CIFs.py extended_convex_hull extended_convex_hull_POSCARS 50GPa')
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
for structure, ID, enthalpy, fitness, pmg_composition in data:
    if pmg_composition.reduced_formula in compositions.keys():
        if compositions[pmg_composition.reduced_formula] < 5:
            compositions[pmg_composition.reduced_formula] += 1
            systems_to_print.append((structure, ID, enthalpy, fitness, pmg_composition))
    else:
        compositions[pmg_composition.reduced_formula] = 1
        systems_to_print.append((structure, ID, enthalpy, fitness, pmg_composition))

# sort according to fitness
systems_to_print.sort(key=lambda tup: tup[3])

for i, (structure, ID, enthalpy, fitness, pmg_composition) in enumerate(systems_to_print):
    print_system(i, structure, ID, enthalpy, fitness, pmg_composition, pressure)

print('Done.')
