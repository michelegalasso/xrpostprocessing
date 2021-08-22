"""
@file:      scxrd_screening.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      2 November 2020
@brief:     Script which uses the class SpectrumAnalyzer.
"""

from xrpostproc.common.SpectrumAnalyzer import SpectrumAnalyzer


# input parameters (MgSiO3)
gatheredPOSCARS = 'gatheredPOSCARS'                 # name of the gatheredPOSCARS file
hkl_file = 'test_P1.hkl'                            # name of the experimental hkl file

analyzer = SpectrumAnalyzer(extended_convex_hull_POSCARS=gatheredPOSCARS, hkl_file=hkl_file)
analyzer.run()
