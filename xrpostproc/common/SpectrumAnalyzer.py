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

from .get_reflections import get_reflections
from .iterator_poscar_file import iterator_poscar_file
from .read_structures import read_structures


class SpectrumAnalyzer(object):
    def __init__(self, exp_pressure: int, th_pressure: int, extended_convex_hull: str,
                 extended_convex_hull_POSCARS: str, spectrum_file: str = None,
                 spectrum_starts: float = None, spectrum_ends: float = None,
                 wavelength: float = None, sigma: float = None, hkl_file: str = None,
                 min_d_spacing: float = None):
        """
        Initializes the class.

        Args:
            exp_pressure:
            th_pressure:
            extended_convex_hull:
            extended_convex_hull_POSCARS:
            spectrum_file:
            spectrum_starts:
            spectrum_ends:
            wavelength:
            sigma:
            hkl_file:
        """
        self.exp_pressure = exp_pressure
        self.th_pressure = th_pressure
        self.deltaP = exp_pressure - th_pressure
        self.extended_convex_hull = extended_convex_hull
        self.extended_convex_hull_POSCARS = extended_convex_hull_POSCARS

        if spectrum_file is not None:
            self.spectrum_file = spectrum_file
            self.spectrum_starts = spectrum_starts
            self.spectrum_ends = spectrum_ends
            self.wavelength = wavelength
            self.sigma = sigma

            for param in ['spectrum_starts', 'spectrum_ends', 'wavelength', 'sigma']:
                if getattr(self, param) is None:
                    raise ValueError(f'Powder mode requires parameter {param}.')
            self.mode = 'powder'
        else:
            self.hkl_file = hkl_file
            self.min_d_spacing = min_d_spacing

            for param in ['hkl_file', 'min_d_spacing']:
                if getattr(self, param) is None:
                    raise ValueError(f'SCXRD mode requires parameter {param}.')
            self.mode = 'scxrd'

    def run(self, match_tol: float = None, individuals: bool = False):
        """
        Creates pictures and cif files in results folder.

        Args:
            match_tol:
            individuals:

        Returns:
            None
        """
        os.mkdir('results')
        print('Processing structures..')

        # pressure correction
        k = (300 / (150 + (22500 + 300 * self.deltaP) ** (1 / 2))) ** (1 / 3)

        if self.mode == 'powder':
            if match_tol is None:
                raise ValueError('Powder mode requires parameter match_tol.')

            spectrum = np.loadtxt(self.spectrum_file)
            exp_angles = spectrum[:, 0]
            exp_intensities = spectrum[:, 1]
            exp_intensities = exp_intensities / exp_intensities.max() * 100

            # initialize calculator
            calculator = XRDCalculator(wavelength=self.wavelength)

            # read files
            data = read_structures(self.extended_convex_hull, self.extended_convex_hull_POSCARS, fixcomp=individuals)

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
        else:
            for poscar_string in iterator_poscar_file(self.extended_convex_hull_POSCARS):
                structure = Structure.from_str(poscar_string, fmt='poscar')
                # structure.lattice = Lattice(np.diag([k, k, k]) @ structure.lattice.matrix)

                reflections = get_reflections(structure, self.min_d_spacing)

                # TODO: compute weighted R-factor based on obtained reflections
