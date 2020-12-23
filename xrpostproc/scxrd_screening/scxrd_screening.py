"""
@file:      scxrd_screening.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      2 November 2020
@brief:     Script which uses the class SpectrumAnalyzer.
"""

from xrpostproc.common.SpectrumAnalyzer import SpectrumAnalyzer


# input parameters (MgSiO3)
exp_pressure = 0        # pressure of the experimental powder
th_pressure = 0         # pressure of the USPEX calculation

extended_convex_hull = 'DOESNOTEXIST.txt'           # name of the extended_convex_hull file
extended_convex_hull_POSCARS = 'gatheredPOSCARD'    # name of the gatheredPOSCARS file
# spectrum_starts = 6.0               # minimum angle for theoretical spectra
# spectrum_ends = 22.0                # maximum angle for theoretical spectra
# wavelength = 0.6199                 # experimental wavelength
# sigma = 0.01                        # parameter for gaussian smearing of peaks
# spectrum_file = 'spectrum.txt'      # name of the experimental spectrum file
# match_tol = 0.25                    # tolerance for matching peaks in degrees
hkl_file = 'test_P1.hkl'            # name of the experimental hkl file

analyzer = SpectrumAnalyzer(exp_pressure, th_pressure, extended_convex_hull, extended_convex_hull_POSCARS,
                            hkl_file=hkl_file)
analyzer.run()
