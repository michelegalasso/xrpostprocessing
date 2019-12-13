"""
@file:      SpectrumAnalyzer.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      4 October 2019
@brief:     Class which implements the fitness function for X-ray optimization.
"""

import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from .iterator_poscar_file import iterator_poscar_file


class SpectrumAnalyzer(object):
    def __init__(self, exp_pressure: int, deltaP: int, spectrum_starts: int, spectrum_ends: int, wavelength: float,
                 sigma: float, spectrum_file: str, gatheredPOSCARS: str):
        self.exp_pressure = exp_pressure
        self.deltaP = deltaP
        self.spectrum_starts = spectrum_starts
        self.spectrum_ends = spectrum_ends
        self.wavelength = wavelength
        self.sigma = sigma
        self.spectrum = np.loadtxt(spectrum_file)
        self.gatheredPOSCARS = gatheredPOSCARS

    def run(self, match_tol: float, factors : list):
        amplitude = self.spectrum_ends - self.spectrum_starts
        exp_angles = self.spectrum[:, 0]
        exp_intensities = self.spectrum[:, 1]
        exp_intensities = exp_intensities / exp_intensities.max() * 100

        # pressure correction
        k = (300 / (150 + (22500 + 300 * self.deltaP) ** (1 / 2))) ** (1 / 3)

        # initialize calculator
        calculator = XRDCalculator(wavelength=self.wavelength)

        print('Analyzing structures..')
        for i, string in enumerate(iterator_poscar_file(self.gatheredPOSCARS)):
            ID = string.split()[0]
            tmp = Structure.from_str(string, fmt='poscar')
            string = tmp.to(fmt='cif', symprec=0.2)
            structure = Structure.from_str(string, fmt='cif')
            structure.lattice = Lattice(np.diag([k, k, k]) @ structure.lattice.matrix)
            composition = structure.composition.iupac_formula.replace(' ', '')

            # ignore pure hydrogen
            if structure.composition.chemical_system == 'H':
                continue

            fitness = 0
            pattern = calculator.get_pattern(structure)
            th_angles = pattern.x[np.logical_and(pattern.x < self.spectrum_ends, pattern.x > self.spectrum_starts)]
            th_intensities = pattern.y[np.logical_and(pattern.x < self.spectrum_ends, pattern.x > self.spectrum_starts)]

            # match corresponding peaks
            exp_matches, th_matches = [], []
            for exp_index, (exp_angle, exp_intensity) in enumerate(zip(exp_angles, exp_intensities)):
                partial = 0
                counter = 0
                for th_index, (th_angle, th_intensity) in enumerate(zip(th_angles, th_intensities)):
                    if np.abs(th_angle - exp_angle) < match_tol and th_intensity > 1:
                        counter += 1
                        exp_matches.append(exp_index)
                        th_matches.append(th_index)
                        factor = self.choose_factor(exp_intensity, factors)
                        partial += factor * np.abs(exp_angle - th_angle) ** 2 / amplitude ** 2
                        partial += factor * np.abs(exp_intensity - th_intensity) ** 2 / 100 ** 2
                # average out in the case when multiple theoretical peaks match to the same experimental peak
                if partial:
                    fitness += partial / counter

            exp_angles_rest = np.delete(exp_angles, exp_matches)
            exp_intensities_rest = np.delete(exp_intensities, exp_matches)
            th_angles_rest = np.delete(th_angles, th_matches)
            th_intensities_rest = np.delete(th_intensities, th_matches)

            # experimental rest
            for angle, intensity in zip(exp_angles_rest, exp_intensities_rest):
                factor = self.choose_factor(intensity, factors)
                fitness += factor * angle ** 2 / amplitude ** 2
                fitness += factor * intensity ** 2 / 100 ** 2

            # theoretical rest
            for angle, intensity in zip(th_angles_rest, th_intensities_rest):
                factor = self.choose_factor(intensity, factors)
                fitness += factor * angle ** 2 / amplitude ** 2
                fitness += factor * intensity ** 2 / 100 ** 2

            # write output
            analyzer = SpacegroupAnalyzer(structure, symprec=0.2)
            try:
                spg = analyzer.get_space_group_number()
                structure.to(filename='{:08.4f}_{}_{}_spg{}_{}GPa.cif'.format(fitness, composition, ID, spg,
                                                                              self.exp_pressure), symprec=0.2)
            except:
                spg = 'ND'
            self.plot_results(exp_angles, exp_intensities, th_angles, th_intensities, ID, fitness, composition, spg)

            print(i + 1)

    def plot_results(self, exp_angles, exp_intensities, th_angles, th_intensities, ID, fitness, composition, spg):
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
        plt.savefig('{:08.4f}_{}_{}_spg{}_{}GPa.png'.format(fitness, composition, ID, spg, self.exp_pressure))
        plt.close()

    @staticmethod
    def choose_factor(intensity, choices):
        if intensity > 90:
            factor = choices[0]
        elif 50 < intensity <= 90:
            factor = choices[1]
        elif 10 < intensity <= 50:
            factor = choices[2]
        elif 1 < intensity <= 10:
            factor = choices[3]
        else:
            factor = choices[4]
        return factor
