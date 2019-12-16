"""
@file:      xr_screening.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      4 October 2019
@brief:     Script which uses the class SpectrumAnalyzer.
"""

from xrpostproc.xr_screening.SpectrumAnalyzer import SpectrumAnalyzer


factors = [
    5,              # fitness factor for I > 90
    1,              # fitness factor for 50 < I <= 90
    0.25,           # fitness factor for 10 < I <= 50
    0.02,           # fitness factor for 1 < I <= 10
    0               # fitness factor for I <= 1
]

# input parameters (BaH)
exp_pressure = 138      # pressure of the experimental powder
th_pressure = 150       # pressure of the USPEX calculation

spectrum_starts = 6.0           # minimum angle for theoretical spectra
spectrum_ends = 22.0            # maximum angle for theoretical spectra
wavelength = 0.6199             # experimental wavelength
sigma = 0.01                    # parameter for gaussian smearing of peaks
spectrum_file = 'spectrum.txt'                                      # name of the experimental spectrum file
extended_convex_hull = 'extended_convex_hull.txt'                   # name of the extended_convex_hull file
extended_convex_hull_POSCARS = 'extended_convex_hull_POSCARS.txt'   # name of the gatheredPOSCARS file
match_tol = 0.25                # tolerance for matching peaks in degrees

analyzer = SpectrumAnalyzer(exp_pressure, th_pressure, spectrum_starts, spectrum_ends, wavelength, sigma, spectrum_file,
                            extended_convex_hull, extended_convex_hull_POSCARS)
analyzer.run(match_tol, factors)
