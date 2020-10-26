"""
@file:      pareto_picture.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      25 March 2020
@brief:     Script which produces an svg picture with Pareto fronts.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from xrpostproc.common.iterator_poscar_file import iterator_poscar_file
from xrpostproc.common.plot_against_experiment import plot_against_experiment

plt.rcParams.update({'font.size': 22})

# input parameters
wavelength = 1.5406
spectrum_file = 'subtracted_KTaWO6_4GPa.txt'
goodStructures = 'goodStructures_4gpa'
goodStructures_POSCARS = 'goodStructures_POSCARS_4gpa'
folder_name = 'pareto_fronts_4gpa'
picture_name = 'pareto_4gpa.png'

poscars_iterator = iterator_poscar_file(goodStructures_POSCARS)

# in this folder I will save the produced pictures
os.mkdir(folder_name)

p_frontsX, p_frontsY = [[]], [[]]
with open(goodStructures, 'r') as f:
    for line in f:
        values = [value.strip() for value in line.split(('|'))]
        if (len(values) > 1) and (values[1] != 'ID'):
            ID = 'EA' + values[1]
            rank = int(values[2])
            enthalpy = float(values[5])
            xraydistance = float(values[6])

            if len(p_frontsX) < rank + 1:
                p_frontsX.append([])
                p_frontsY.append([])
            p_frontsY[rank].append(enthalpy)
            p_frontsX[rank].append(xraydistance)

            poscar_string = next(poscars_iterator)
            if ID == poscar_string.split()[0]:
                if rank < 3:
                    plot_against_experiment(poscar_string, spectrum_file, wavelength, folder_name, rank + 1,
                                            xraydistance, enthalpy)
            else:
                raise ValueError('Structure IDs in the two input files do not match.')

plt.figure(figsize=(16, 9))
ordinal = lambda n: "%d%s" % (n, "tsnrhtdd"[(n / 10 % 10 != 1)*(n % 10 < 4) * n % 10::4])

for i in range(3):
    indices = np.argsort(p_frontsX[i])
    p_frontX = np.array(p_frontsX[i])[indices]
    p_frontY = np.array(p_frontsY[i])[indices]
    plt.plot(p_frontX, p_frontY, 'o-', label='{} Pareto front'.format(ordinal(i + 1)))

X, Y = [], []
for i, (frontX, frontY) in enumerate(zip(p_frontsX, p_frontsY)):
    if i >= 3:
        for xraydistance, enthalpy in zip(frontX, frontY):
            if True:
                X.append(xraydistance)
                Y.append(enthalpy)

plt.plot(X, Y, 'o')
plt.ylabel('Enthalpy (eV/f.u.)')
plt.xlabel('Distance between calculated and experimental X-ray spectrum')
# plt.xlim(0.132, 0.155)
# plt.ylim(-0.1, 2.3)
plt.legend()

# plt.show()
plt.savefig(picture_name)
plt.close()
