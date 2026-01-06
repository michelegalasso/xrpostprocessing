"""
@file:      powder_xrd_screening.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      4 October 2019
@brief:     Script which uses the class SpectrumAnalyzer.
"""

import numpy as np

from xrpostproc.common.SpectrumAnalyzer import SpectrumAnalyzer


# input parameters (BaH)
gathered_POSCARS = 'UH_50_extended_convex_hull_POSCARS.txt'     # file with all structures in the POSCAR format
spectrum_file = 'UH_50_powder.txt'                              # name of the experimental spectrum file
spectrum_starts = 6.0                                           # minimum angle for theoretical spectra
spectrum_ends = 22.0                                            # maximum angle for theoretical spectra
k_coeffs = np.arange(0.91, 1.09, 0.005)                         # coefficients for isotropic cell deformation
sigma = 0.01                                                    # parameter for gaussian smearing of peaks
wavelength = 0.4130                                             # experimental wavelength
match_tol = 0.25                                                # tolerance for matching peaks in degrees
th_peaks_penalty = 5.0                                          # additional penalty for non-matching th_peaks

analyzer = SpectrumAnalyzer(gathered_POSCARS=gathered_POSCARS, spectrum_file=spectrum_file,
                            spectrum_starts=spectrum_starts, spectrum_ends=spectrum_ends, wavelength=wavelength,
                            sigma=sigma, th_peaks_penalty=th_peaks_penalty)
analyzer.run(k_coeffs=k_coeffs, match_tol=match_tol)
