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


if len(sys.argv) != 4:
    print('ERROR: wrong number of arguments.')
    print('Usage: python3 fixComp.py [file with enthalpies] [gatheredPOSCARS] [pressure]')
    print('Example: python3 fix_comp.py Individuals gatheredPOSCARS 50GPa')
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
            if ID == values[1]:
                structure = Structure.from_str(string, fmt='poscar')
                realFitness = float(values[5 + len(composition)]) / structure.num_sites
                data.append((ID, realFitness, structure))
            else:
                raise IOError('Structures in {} do not match data in {}'.format(sys.argv[2], sys.argv[1]))

# sort according to real_fitness
data.sort(key=lambda tup: tup[1])
pressure = sys.argv[3]

os.mkdir('results')
for i, (ID, realFitness, structure) in enumerate(data):
    nospaceIupacFormula = structure.composition.iupac_formula.replace(' ', '')
    symm = spglib.get_symmetry_dataset((structure.lattice.matrix, structure.frac_coords, structure.atomic_numbers),
                                       symprec=0.2)['number']
    filename = '{}_{:f}_{}_{}_ID{}_{}.cif'.format(i + 1, realFitness, nospaceIupacFormula, symm, ID, pressure)
    structure.to(filename=os.path.join('results', filename), symprec=0.2)

print('Done.')
