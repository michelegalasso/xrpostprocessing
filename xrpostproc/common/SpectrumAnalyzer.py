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
    def __init__(self, exp_pressure: int = None, th_pressure: int = None, k_coeffs: np.ndarray = None,
                 extended_convex_hull: str = None, extended_convex_hull_POSCARS: str = None,
                 spectrum_file: str = None, spectrum_starts: float = None, spectrum_ends: float = None,
                 wavelength: float = None, sigma: float = None, B0: float = 300,
                 dB0: float = 3, hkl_file: str = None):
        """
        Initializes the class.

        Args:
            exp_pressure:
            th_pressure:
            k_coeffs:
            extended_convex_hull:
            extended_convex_hull_POSCARS:
            spectrum_file:
            spectrum_starts:
            spectrum_ends:
            wavelength:
            sigma:
            hkl_file:
            B0:
            dB0:
        """
        self.extended_convex_hull_POSCARS = extended_convex_hull_POSCARS

        if spectrum_file is not None:
            self.th_pressure = th_pressure

            if exp_pressure is not None and k_coeffs is None:
                self.exp_pressure = exp_pressure
                self.k_coeffs = None
            elif exp_pressure is None and k_coeffs is not None:
                self.exp_pressure = None
                self.k_coeffs = k_coeffs
            else:
                raise ValueError("Please give either exp_pressure or k_coeffs.")

            self.extended_convex_hull = extended_convex_hull
            self.spectrum_file = spectrum_file
            self.spectrum_starts = spectrum_starts
            self.spectrum_ends = spectrum_ends
            self.wavelength = wavelength
            self.sigma = sigma
            self.B0 = B0
            self.dB0 = dB0

            for param in ['spectrum_starts', 'spectrum_ends', 'wavelength', 'sigma']:
                if getattr(self, param) is None:
                    raise ValueError(f'Powder mode requires parameter {param}.')
            self.mode = 'powder'
        else:
            self.hkl_file = hkl_file

            for param in ['hkl_file']:
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
        print('Processing structures..')

        if self.mode == 'powder':
            os.mkdir('results')
            if match_tol is None:
                raise ValueError('Powder mode requires parameter match_tol.')

            k_coeffs = []
            if self.k_coeffs is None:
                deltaP = self.exp_pressure - self.th_pressure
                num = (self.B0 * (self.dB0 - 1))
                den = (self.B0 * (self.dB0 - 2)) + np.sqrt(self.B0 ** 2 + 2 * num * deltaP)
                k_coeffs.append((num / den) ** (1 / 3))
            else:
                k_coeffs = self.k_coeffs

            spectrum = np.loadtxt(self.spectrum_file)
            exp_angles = spectrum[:, 0]
            exp_intensities = spectrum[:, 1]
            exp_intensities = exp_intensities / exp_intensities.max() * 100

            # initialize calculator
            calculator = XRDCalculator(wavelength=self.wavelength)

            # read files
            data = read_structures(self.extended_convex_hull, self.extended_convex_hull_POSCARS, fixcomp=individuals)

            for i, (tmp, ID, enth, fit, pmg_comp) in enumerate(data):
                # structure symmetrization
                string = tmp.to(fmt='cif', symprec=0.2)

                # structure scaling
                fitnesses = []
                for k in k_coeffs:
                    structure = Structure.from_str(string, fmt='cif')
                    structure.lattice = Lattice(np.diag([k, k, k]) @ structure.lattice.matrix)

                    # ignore pure hydrogen
                    if structure.composition.chemical_system == 'H':
                        continue

                    fitness = 0
                    pattern = calculator.get_pattern(structure=structure,
                                                     two_theta_range=(self.spectrum_starts, self.spectrum_ends))
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

                    fitnesses.append(fitness)

                # get structure with minimum fitness
                structure = Structure.from_str(string, fmt='cif')

                # ignore pure hydrogen
                if structure.composition.chemical_system == 'H':
                    continue

                # retrieve minimum fitness value and corresponding pressure
                j = np.argmin(fitnesses)
                fitness = fitnesses[j]

                # scale the structure
                structure.lattice = Lattice(np.diag([k_coeffs[j], k_coeffs[j], k_coeffs[j]]) @ structure.lattice.matrix)
                pattern = calculator.get_pattern(structure=structure,
                                                 two_theta_range=(self.spectrum_starts, self.spectrum_ends))

                # angles and intensities to plot
                th_angles = pattern.x
                th_intensities = pattern.y

                # write cif
                composition = structure.composition.iupac_formula.replace(' ', '')
                dtset = spglib.get_symmetry_dataset((structure.lattice.matrix, structure.frac_coords,
                                                     structure.atomic_numbers), symprec=0.2)
                if dtset is not None:
                    filename = '{:08.4f}_EA{}_{}_{}_{:.3f}_{}GPa_spg{}'.format(fitness, ID, fit, composition,
                                                                           k_coeffs[j], self.th_pressure, dtset.number)
                    structure.to(filename=os.path.join('results', filename + '.cif'), symprec=0.2)
                else:
                    filename = '{:08.4f}_EA{}_{}_{}_{:.3f}_{}GPa_spgND'.format(fitness, ID, fit, composition,
                                                                           k_coeffs[j], self.th_pressure)

                # write png
                plt.figure(figsize=(16, 9))
                x = np.arange(self.spectrum_starts, self.spectrum_ends, 0.01)
                y_exp = np.zeros_like(x)

                for angle, intensity in zip(exp_angles, exp_intensities):
                    y_exp += intensity * np.exp(-((x - angle) ** 2) / (2 * self.sigma ** 2))

                plt.plot(x, y_exp, label='experimental')

                plt.stem(th_angles, th_intensities, 'r', markerfmt='None', basefmt='None',
                         label='{} predicted'.format(ID))

                plt.legend()
                plt.savefig(os.path.join('results', filename + '.png'))
                plt.close()

                if len(k_coeffs) > 1:
                    # write png
                    plt.figure(figsize=(16, 9))
                    plt.plot(k_coeffs, fitnesses, marker=".")
                    plt.xlabel("k coefficient")
                    plt.ylabel("fitness")
                    plt.savefig(os.path.join('results', filename + '_k.png'))

                print(i + 1)
        else:
            results = {}
            for poscar_string in iterator_poscar_file(self.extended_convex_hull_POSCARS):
                ID = poscar_string.split()[0]
                structure = Structure.from_str(poscar_string, fmt='poscar')

                # read hkl file
                exp_reflections = []
                min_d_spacing = 10000
                with open(self.hkl_file, 'r') as f:
                    for line in f:
                        values = line.split()

                        # the hkl file terminates with all zeros
                        if values == ['0', '0', '0', '0.00', '0.00']:
                            break

                        i_hkl = float(values[3])
                        if i_hkl > 0:
                            hkl = (int(values[0]), int(values[1]), int(values[2]))
                            sigma_hkl = float(values[4])

                            d_spacing = structure.lattice.d_hkl(list(hkl))
                            if d_spacing < min_d_spacing:
                                min_d_spacing = d_spacing

                            exp_reflections.append([i_hkl, hkl, sigma_hkl])

                min_d_spacing = np.floor(min_d_spacing * 1000) / 1000
                th_reflections = get_reflections(structure, min_d_spacing)

                exp_reflections = np.array(exp_reflections, dtype=object)
                th_reflections = np.array(th_reflections, dtype=object)

                # scale theoretical intensities according to experimental maximum
                th_reflections[:, 0] = th_reflections[:, 0] / max(th_reflections[:, 0]) * max(exp_reflections[:, 0])

                numerator = 0
                denominator = 0
                for i_hkl, hkl, sigma_hkl in exp_reflections:
                    th_hkls = [r[1] for r in th_reflections]

                    if hkl in th_hkls:
                        index = th_hkls.index(hkl)
                    elif (-hkl[0], -hkl[1], -hkl[2]) in th_hkls:
                        index = th_hkls.index((-hkl[0], -hkl[1], -hkl[2]))
                    else:
                        raise ValueError(f'Experimental reflection {hkl} with intensity {i_hkl} not found in theory.')

                    numerator += (1 / sigma_hkl ** 2) * (i_hkl - th_reflections[index][0]) ** 2
                    denominator += (1 / sigma_hkl ** 2) * i_hkl ** 2

                wR = np.sqrt(numerator / denominator)
                results[ID] = wR

                print(ID)

            results = {k: v for k, v in sorted(results.items(), key=lambda item: item[1])}
            with open('results.txt', 'w') as f:
                for ID, wR in results.items():
                    f.write(f'{ID:6s}  {wR:6.4f}\n')
