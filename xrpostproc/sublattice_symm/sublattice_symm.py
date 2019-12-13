import os
import sys
import spglib

from pymatgen.core.structure import Structure


def iterator_poscar_file(file_name):
    with open(file_name, 'r') as f:
        while True:
            current_line = f.readline()
            if current_line == '':
                return
            lines = [current_line]
            for i in range(6):
                lines.append(f.readline())
            atom_nums = lines[-1].split()
            num_atoms = 0
            for i in range(len(atom_nums)):
                num_atoms += int(atom_nums[i])

            for i in range(num_atoms + 1):
                lines.append(f.readline())
            total_line = ''
            for line in lines:
                total_line += line
            yield total_line


def print_system(i, structure, enthalpy, fitness, full_composition, pressure):
    dtset = spglib.get_symmetry_dataset((structure.lattice.matrix, structure.frac_coords, structure.atomic_numbers),
                                        symprec=0.2)
    filename = '{}_{}_{}_{}_{}_{}.cif'.format(i + 1, fitness, enthalpy, full_composition, pressure, dtset['number'])
    structure.to(filename=os.path.join('results', filename), symprec=0.2)


if len(sys.argv) != 4:
    print('ERROR: wrong number of arguments.')
    print('Usage: python3 split_CIFs.py [Individuals] [gatheredPOSCARS] [pressure]')
    print('Example: python3 split_CIFs.py Individuals gatheredPOSCARS 50GPa')
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
                full_composition = structure.composition.reduced_formula
                structure.remove_species(['H'])

                if len(structure.frac_coords) != 0:
                    enthalpy = float(values[3 + len(composition)])
                    fitness = float(values[5 + len(composition)])
                    data.append((structure, enthalpy, fitness, full_composition))
            else:
                raise IOError('Structures in {} do not match data in {}'.format(sys.argv[2], sys.argv[1]))

# sort according to enthalpy
data.sort(key=lambda tup: tup[1])
pressure = sys.argv[3]
systems_to_print = []
compositions = {}

os.mkdir('results')
for structure, enthalpy, fitness, full_composition in data:
    if full_composition in compositions.keys():
        if compositions[full_composition] < 5:
            compositions[full_composition] += 1
            systems_to_print.append((structure, enthalpy, fitness, full_composition))
    else:
        compositions[full_composition] = 1
        systems_to_print.append((structure, enthalpy, fitness, full_composition))

# sort according to fitness
systems_to_print.sort(key=lambda tup: tup[2])

for i, (structure, enthalpy, fitness, full_composition) in enumerate(systems_to_print):
    print_system(i, structure, enthalpy, fitness, full_composition, pressure)

print('Done.')
