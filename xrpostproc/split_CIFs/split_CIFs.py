import os
import sys

from pymatgen.core.structure import Structure

from xrpostproc.common.iterator_poscar_file import iterator_poscar_file
from xrpostproc.common.print_system import print_system


if len(sys.argv) != 4:
    print('ERROR: wrong number of arguments.')
    print('Usage: python split_CIFs.py [extended_convex_hull] [extended_convex_hull_POSCARS] [pressure]')
    print('Example: python split_CIFs.py extended_convex_hull extended_convex_hull_POSCARS 50GPa')
    sys.exit()

# read files
iterator = iterator_poscar_file(sys.argv[2])
data = []
with open(sys.argv[1], 'r') as f:
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
                enthalpy = float(values[3 + len(composition)])
                fitness = float(values[5 + len(composition)])
                data.append((structure, enthalpy, fitness))
            else:
                raise IOError('Structures in {} do not match data in {}'.format(sys.argv[2], sys.argv[1]))

# sort according to enthalpy
data.sort(key=lambda tup: tup[1])
pressure = sys.argv[3]
systems_to_print = []
compositions = {}

os.mkdir('results')
for structure, enthalpy, fitness in data:
    if structure.composition.reduced_formula in compositions.keys():
        if compositions[structure.composition.reduced_formula] < 5:
            compositions[structure.composition.reduced_formula] += 1
            systems_to_print.append((structure, enthalpy, fitness))
    else:
        compositions[structure.composition.reduced_formula] = 1
        systems_to_print.append((structure, enthalpy, fitness))

# sort according to fitness
systems_to_print.sort(key=lambda tup: tup[2])

for i, (structure, enthalpy, fitness) in enumerate(systems_to_print):
    print_system(i, structure, enthalpy, fitness, pressure)

print('Done.')
