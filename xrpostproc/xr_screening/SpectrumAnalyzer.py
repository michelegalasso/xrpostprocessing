"""
@file:      SpectrumAnalyzer.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      4 October 2019
@brief:     Class which implements the fitness function for X-ray optimization.
"""

import os
import spglib
import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.analysis.diffraction.xrd import XRDCalculator

from ..common.read_structures import read_structures


class SpectrumAnalyzer(object):
    def __init__(self, exp_pressure: int, th_pressure: int, spectrum_starts: float, spectrum_ends: float,
                 wavelength: float, sigma: float, spectrum_file: str, extended_convex_hull: str,
                 extended_convex_hull_POSCARS: str):
        self.exp_pressure = exp_pressure
        self.th_pressure = th_pressure
        self.deltaP = exp_pressure - th_pressure
        self.spectrum_starts = spectrum_starts
        self.spectrum_ends = spectrum_ends
        self.wavelength = wavelength
        self.sigma = sigma
        self.spectrum = np.loadtxt(spectrum_file)
        self.extended_convex_hull = extended_convex_hull
        self.extended_convex_hull_POSCARS = extended_convex_hull_POSCARS

    def run(self, match_tol: float, factors : list, individuals: bool = False):
        exp_angles = self.spectrum[:, 0]
        exp_intensities = self.spectrum[:, 1]
        exp_intensities = exp_intensities / exp_intensities.max() * 100

        # pressure correction
        k = (300 / (150 + (22500 + 300 * self.deltaP) ** (1 / 2))) ** (1 / 3)

        # initialize calculator
        calculator = XRDCalculator(wavelength=self.wavelength)

        # read files
        data = read_structures(self.extended_convex_hull, self.extended_convex_hull_POSCARS, fixcomp=individuals)

        os.mkdir('results')
        print('Processing structures..')
        for i, (tmp, ID, enth, fit, pmg_comp) in enumerate(data):
            string = tmp.to(fmt='cif', symprec=0.2)
            structure = Structure.from_str(string, fmt='cif')
            structure.lattice = Lattice(np.diag([k, k, k]) @ structure.lattice.matrix)

            # ignore pure hydrogen
            if structure.composition.chemical_system == 'H':
                continue

            fitness = 0
            pattern = calculator.get_pattern(structure, two_theta_range=(self.spectrum_starts, self.spectrum_ends))
            th_angles = pattern.x
            th_intensities = pattern.y

            # match corresponding peaks
            exp_matches, th_matches = [], []
            for exp_index, (exp_angle, exp_intensity) in enumerate(zip(exp_angles, exp_intensities)):
                partial = 0
                counter = 0
                for th_index, (th_angle, th_intensity) in enumerate(zip(th_angles, th_intensities)):
                    if np.abs(th_angle - exp_angle) < match_tol:
                        counter += 1
                        exp_matches.append(exp_index)
                        th_matches.append(th_index)
                        partial += ((exp_intensity - th_intensity) / 100) ** 2 * (exp_intensity / 100) ** 2
                # average out in the case when multiple theoretical peaks match to the same experimental peak
                if partial:
                    fitness += partial / counter

            exp_intensities_rest = np.delete(exp_intensities, exp_matches)
            th_intensities_rest = np.delete(th_intensities, th_matches)

            # experimental rest
            for intensity in exp_intensities_rest:
                fitness += (intensity / 100) ** 2

            # theoretical rest
            for intensity in th_intensities_rest:
                fitness += (intensity / 100) ** 2

            # write cif
            composition = structure.composition.iupac_formula.replace(' ', '')
            dtset = spglib.get_symmetry_dataset((structure.lattice.matrix, structure.frac_coords,
                                                 structure.atomic_numbers), symprec=0.2)
            if dtset is not None:
                filename = '{:08.4f}_EA{}_{}_{}_{}GPa_spg{}'.format(fitness, ID, fit, composition, self.th_pressure,
                                                                    dtset['number'])
                structure.to(filename=os.path.join('results', filename + '.cif'), symprec=0.2)
            else:
                filename = '{:08.4f}_EA{}_{}_{}_{}GPa_spgND'.format(fitness, ID, fit, composition, self.th_pressure)

            # write txt
            # with open('structures.txt', 'a') as f:
            #     f.write('{:>6}    {:8.4f}    {:8.4f}\n'.format(ID, fit, fitness))

            # write png
            plt.figure(figsize=(16, 9))
            x = np.arange(self.spectrum_starts, self.spectrum_ends, 0.01)
            y_exp = np.zeros_like(x)

            for angle, intensity in zip(exp_angles, exp_intensities):
                y_exp += intensity * np.exp(-((x - angle) ** 2) / (2 * self.sigma ** 2))

            plt.plot(x, y_exp, label='experimental')

            y_th = np.zeros_like(x)
            for angle, intensity in zip(th_angles, th_intensities):
                y_th += intensity * np.exp(-((x - angle) ** 2) / (2 * self.sigma ** 2))

            plt.plot(x, y_th, label='calculated')

            plt.legend()
            plt.savefig(os.path.join('results', filename + '.png'))
            plt.close()

            print(i + 1)
