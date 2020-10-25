"""
@file:      plot_against_experiment.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      25 March 2020
@brief:     Script which plots theoretical peaks against experimental spectra.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core.structure import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator

from xrpostproc.common.iterator_poscar_file import iterator_poscar_file


# load experimental spectrum
exp_spectrum = np.loadtxt('XRD_135gpa_pure_BaH12.txt')
exp_angles = exp_spectrum[:, 0]
exp_intensities = exp_spectrum[:, 1] / max(exp_spectrum[:, 1]) * 100

# in this folder I will save the produced pictures
os.mkdir('plots')

for poscar in iterator_poscar_file('goodStructures_POSCARS'):
    ID = poscar.split()[0]

    tmp = Structure.from_str(poscar, fmt='POSCAR')
    string = tmp.to(fmt='cif', symprec=0.2)
    structure = Structure.from_str(string, fmt='cif')

    calculator = XRDCalculator(wavelength=0.6199)
    th_spectrum = calculator.get_pattern(structure, two_theta_range=(min(exp_angles), max(exp_angles)))
    th_angles = th_spectrum.x
    th_intensities = th_spectrum.y

    plt.rcParams.update({'font.size': 22})
    plt.figure(figsize=(16, 9))
    plt.plot(exp_angles, exp_intensities, label='experimental')
    plt.stem(th_angles, th_intensities, 'r', markerfmt='None', basefmt='None', use_line_collection=True,
             label='{} predicted'.format(ID))
    plt.ylabel('2θ')
    plt.ylabel('Intensity')
    plt.legend()

    plt.savefig('plots/{}.png'.format(ID))
    # plt.show()
