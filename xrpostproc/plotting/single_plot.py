"""
@file:      single_plot.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      25 March 2020
@brief:     Script which plots one structure against experiment.
"""

from xrpostproc.common.plot_against_experiment import plot_against_experiment

# input parameters
wavelength = 0.2952
poscar_file = 'chlorine-31-gpa-refined.vasp'
exp_spectrum_file = 'XRD_30gpa_NaCl.txt'

with open(poscar_file, 'r') as f:
    poscar_string = f.read()

plot_against_experiment(poscar_string, exp_spectrum_file, wavelength)
