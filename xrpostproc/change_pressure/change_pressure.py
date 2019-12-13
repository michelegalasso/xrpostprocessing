import os
import sys
import numpy as np

from tqdm import tqdm
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

if len(sys.argv) != 3:
    print('ERROR: wrong number of arguments.')
    print('Usage: python3 change_pressure.py [initial_pressure] [final_pressure]')
    print('Example: python3 change_pressure.py 50 58')
    sys.exit()

initial_pressure = int(sys.argv[1])
final_pressure = int(sys.argv[2])
deltaP = final_pressure - initial_pressure
k = (300 / (150 + (22500 + 300 * deltaP) ** (1 / 2))) ** (1 / 3)

cif_files = [f for f in os.listdir('.') if '.cif' in f]

for f in tqdm(cif_files):
    structure = Structure.from_file(f)
    structure.lattice = Lattice(np.diag([k, k, k]) @ structure.lattice.matrix)

    filename = f.strip('.cif') + '_to' + str(final_pressure) + 'GPa' + '.cif'
    try:
        structure.to(filename=filename, symprec=0.2)
    except TypeError:
        print('WARNING: {} failed to be created.'.format(filename))

print('Done.')
