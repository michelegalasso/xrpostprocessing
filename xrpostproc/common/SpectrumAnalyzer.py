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
                 wavelength: float = None, sigma: float = None, hkl_file: str = None):
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

        # pressure correction
        k = (300 / (150 + (22500 + 300 * self.deltaP) ** (1 / 2))) ** (1 / 3)

        if self.mode == 'powder':
            os.mkdir('results')
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
            results = {}
            for poscar_string in iterator_poscar_file(self.extended_convex_hull_POSCARS):
                ID = poscar_string.split()[0]
                structure = Structure.from_str(poscar_string, fmt='poscar')
                # structure.lattice = Lattice(np.diag([k, k, k]) @ structure.lattice.matrix)

                # read hkl file
                exp_reflections = []
                min_d_spacing = 10000
                with open(self.hkl_file, 'r') as f:
                    for line in f:
                        values = line.split()

                        # the hkl file terminates with all zeros
                        if values == ['0', '0', '0', '0.00', '0.00']:
                            break

                        hkl = (int(values[0]), int(values[1]), int(values[2]))
                        i_hkl = float(values[3])
                        sigma_hkl = float(values[4])

                        d_spacing = structure.lattice.d_hkl(list(hkl))
                        if d_spacing < min_d_spacing:
                            min_d_spacing = d_spacing

                        exp_reflections.append([i_hkl, hkl, sigma_hkl])

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
