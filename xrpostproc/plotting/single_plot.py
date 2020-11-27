"""
@file:      single_plot.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      25 March 2020
@brief:     Script which plots one structure against experiment.
"""

from xrpostproc.common.plot_against_experiment import plot_against_experiment

# input parameters
wavelength = 1.5406
poscar_file = 'EA2763.vasp'
exp_spectrum_file = 'subtracted_KTaWO6_4GPa.txt'

with open(poscar_file, 'r') as f:
    poscar_string = f.read()

plot_against_experiment(poscar_string, exp_spectrum_file, wavelength)
