"""
@file:      pareto_picture.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      25 March 2020
@brief:     Script which produces an svg picture with Pareto fronts.
"""

import matplotlib.pyplot as plt

from xrpostproc.pareto_picture.pareto_front import pareto_front


X = []
Y = []
p_frontX, p_frontY = [], []
with open('structures.txt', 'r') as f:
    for current_line in f.readlines():
        ID = int(current_line.split()[0])
        enthalpy = float(current_line.split()[1])
        xraydistance = float(current_line.split()[2])
        if xraydistance < 25 and enthalpy < 0.4:
            Y.append(enthalpy)
            X.append(xraydistance)

plt.figure(figsize=(16, 9))
plt.plot(X, Y, 'o')

ordinal = lambda n: "%d%s" % (n, "tsnrhtdd"[(n / 10 % 10 != 1)*(n % 10 < 4) * n % 10::4])

for i in range(3):
    if len(X) != 0:
        p_frontX, p_frontY, X, Y = pareto_front(X, Y, maxX=False, maxY=False)
        plt.plot(p_frontX, p_frontY, 'o-', label='{} Pareto front'.format(ordinal(i + 1)))

        # print('{} FRONT'.format(ordinal(i + 1)))
        # for a, b in zip(p_frontX, p_frontY):
        #     print('{}    {}'.format(a, b))
        # print('\n')

plt.ylabel('Enthalpy of formation (eV/atom)')
plt.xlabel('Distance between calculated and experimental X-ray spectrum')
plt.xlim(0, 26)
plt.ylim(-0.12, 0.42)
plt.legend()

# plt.show()
plt.savefig('pareto.png', dpi=300)
