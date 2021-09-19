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
work_dir = 'KTaWO6_0gpa'
spectrum_file = 'subtracted_KTaWO6_0GPa.txt'        # None for single crystal

if spectrum_file:
    poscars_iterator = iterator_poscar_file(os.path.join(work_dir, 'goodStructures_POSCARS'))

    # in this folder I will save the produced pictures
    os.mkdir(os.path.join(work_dir, 'pareto_fronts'))

p_frontsX, p_frontsY = [[]], [[]]
with open(os.path.join(work_dir, 'goodStructures'), 'r') as f:
    for line in f:
        values = [value.strip() for value in line.split(('|'))]
        if (len(values) > 1):
            if (values[1] == 'ID'):
                xrd_index = values.index('xraydistance')
            else:
                ID = 'EA' + values[1]
                rank = int(values[2])
                enthalpy = float(values[5])
                xraydistance = float(values[xrd_index])

                if len(p_frontsX) < rank + 1:
                    p_frontsX.append([])
                    p_frontsY.append([])
                p_frontsY[rank].append(enthalpy)
                p_frontsX[rank].append(xraydistance)

                if spectrum_file:
                    poscar_string = next(poscars_iterator)
                    if ID == poscar_string.split()[0]:
                        if rank < 3:
                            plot_against_experiment(poscar_string, os.path.join(work_dir, spectrum_file), wavelength,
                                                    os.path.join(work_dir, 'pareto_fronts'), rank + 1,
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
plt.ylabel('Enthalpy (eV)')
plt.xlabel('$F$')
# plt.xlim(0, 4)
# plt.ylim(-53, -20)
plt.legend(loc='upper right', framealpha=1)

# plt.show()
plt.savefig(os.path.join(work_dir, 'pareto.png'))
plt.close()
