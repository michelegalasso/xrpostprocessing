"""
@file:      plot_against_experiment.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      25 March 2020
@brief:     Script which plots theoretical peaks against experimental spectra.
"""

import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core.structure import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator


def plot_against_experiment(poscar_string, exp_spectrum_file, wavelength, folder_name=None, front=None,
                            xraydistance=None, enthalpy=None):
    # get structure ID
    ID = poscar_string.split()[0]

    # save poscar file
    if (folder_name is not None) and (front is not None) and (xraydistance is not None) and (enthalpy is not None):
        with open('{}/{}_{:6.3f}_{:6.3f}_{}.vasp'.format(folder_name, front, xraydistance, enthalpy, ID), 'w') as f:
            f.write(poscar_string)

    # load experimental spectrum
    exp_spectrum = np.loadtxt(exp_spectrum_file)
    exp_angles = exp_spectrum[:, 0]
    exp_intensities = exp_spectrum[:, 1] / max(exp_spectrum[:, 1]) * 100

    tmp = Structure.from_str(poscar_string, fmt='POSCAR')
    string = tmp.to(fmt='cif', symprec=0.2)
    structure = Structure.from_str(string, fmt='cif')

    calculator = XRDCalculator(wavelength)
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

    # plt.show()
    if (folder_name is not None) and (front is not None) and (xraydistance is not None) and (enthalpy is not None):
        plt.savefig('{}/{}_{:6.3f}_{:6.3f}_{}.png'.format(folder_name, front, xraydistance, enthalpy, ID))
    else:
        plt.savefig('{}.png'.format(ID))
    plt.close()
