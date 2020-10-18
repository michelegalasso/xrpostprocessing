"""
@file:      pareto_picture.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      25 March 2020
@brief:     Script which produces an svg picture with Pareto fronts.
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 22})


p_frontsX, p_frontsY = [[]], [[]]
with open('goodStructures_KTaWO6_4gpa', 'r') as f:
    for line in f:
        values = [value.strip() for value in line.split(('|'))]
        if (len(values) > 1) and (values[1] != 'ID'):
            rank = int(values[2])
            enthalpy = float(values[5])
            xraydistance = float(values[6])

            if len(p_frontsX) < rank + 1:
                p_frontsX.append([])
                p_frontsY.append([])
            p_frontsY[rank].append(enthalpy)
            p_frontsX[rank].append(xraydistance)

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

plt.show()
# plt.savefig('pareto_4gpa.png')
