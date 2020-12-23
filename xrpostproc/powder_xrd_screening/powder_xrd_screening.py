"""
@file:      powder_xrd_screening.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      4 October 2019
@brief:     Script which uses the class SpectrumAnalyzer.
"""

from xrpostproc.common.SpectrumAnalyzer import SpectrumAnalyzer


# input parameters (BaH)
exp_pressure = 138      # pressure of the experimental powder
th_pressure = 150       # pressure of the USPEX calculation

extended_convex_hull = 'extended_convex_hull.txt'                   # name of the extended_convex_hull file
extended_convex_hull_POSCARS = 'extended_convex_hull_POSCARS.txt'   # name of the gatheredPOSCARS file
spectrum_file = 'spectrum.txt'      # name of the experimental spectrum file
spectrum_starts = 6.0               # minimum angle for theoretical spectra
spectrum_ends = 22.0                # maximum angle for theoretical spectra
wavelength = 0.6199                 # experimental wavelength
sigma = 0.01                        # parameter for gaussian smearing of peaks
match_tol = 0.25                    # tolerance for matching peaks in degrees

# set to True if you are using Individuals instead of extended_convex_hull
individuals = False

analyzer = SpectrumAnalyzer(exp_pressure, th_pressure, extended_convex_hull, extended_convex_hull_POSCARS,
                            spectrum_file, spectrum_starts, spectrum_ends, wavelength, sigma)
analyzer.run(match_tol, individuals)
