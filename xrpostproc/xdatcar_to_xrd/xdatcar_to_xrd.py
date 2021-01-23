"""
@file:      xdatcar_to_xrd.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      5 January 2021
@brief:     Script which reads XDATCAR and generates spectra.
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core.structure import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator

# from xrpostproc.common.iterator_poscar_file import iterator_poscar_file


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
            total_line = ""
            for line in lines:
                total_line += line
            yield total_line[0:-1]


# increase font in plots
plt.rcParams.update({'font.size': 22})

parser = argparse.ArgumentParser()
parser.add_argument('-w', '--wavelength', help='wavelength of the inciding radiation', type=float)
parser.add_argument('-i', '--input', help='name of the input XDATCAR file', type=str)
parser.add_argument('-m', '--minangle', help='minimum angle of the generated XRD spectra', type=float)
parser.add_argument('-M', '--maxangle', help='maximum angle of the generated XRD spectra', type=float)

args = parser.parse_args()

os.mkdir('results')
calculator = XRDCalculator(args.wavelength)

for i, poscar_string in enumerate(iterator_poscar_file(args.input)):
    # symmetrize with tolerance 0.2
    tmp = Structure.from_str(poscar_string, fmt='poscar')
    cif_string = tmp.to(fmt='cif', symprec=0.2)
    structure = Structure.from_str(cif_string, fmt='cif')

    # get pattern
    pattern = calculator.get_pattern(structure, two_theta_range=(args.minangle, args.maxangle))

    # generate picture
    plt.figure(figsize=(16, 9))
    plt.stem(pattern.x, pattern.y, 'r', markerfmt='None', basefmt='None', use_line_collection=True)
    plt.xlim(args.minangle, args.maxangle)
    plt.savefig(f'results/{i:03d}.png')
    plt.close()

    # generate spectrum file
    with open(f'results/{i:03d}.txt', 'w') as f:
        for angle, intensity in zip(pattern.x, pattern.y):
            f.write(f'{angle:5.2f}  {intensity:10.4f}\n')

total_x = []
total_y = []

for filename in os.listdir('results'):
    if filename.endswith('txt'):
        with open(os.path.join('results', filename), 'r') as f:
            for line in f:
                content = line.split()
                angle = float(content[0])
                intensity = float(content[1])

                if angle in total_x:
                    ind = total_x.index(angle)
                    total_y[ind] += intensity
                else:
                    total_x.append(angle)
                    total_y.append(intensity)

# sort avg spectrum with increasing angle
total_x = np.array(total_x)
total_y = np.array(total_y)

indices = np.argsort(total_x)
avg_x = total_x[indices]
avg_y = total_y[indices] / max(total_y) * 100

# generate avg picture
plt.figure(figsize=(16, 9))
plt.stem(avg_x, avg_y, 'r', markerfmt='None', basefmt='None', use_line_collection=True)
plt.savefig('results/AVERAGE.png')
plt.close()

# generate avg spectrum file
with open('results/AVERAGE.txt', 'w') as f:
    for angle, intensity in zip(avg_x, avg_y):
        f.write(f'{angle:5.2f}  {intensity:10.4f}\n')
