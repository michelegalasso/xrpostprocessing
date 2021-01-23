"""
@file:      get_structure_factors.py
@author:    Michele Galasso
@contact:   m.galasso@yandex.com
@date:      3 November 2020
@brief:     Function adapted from pymatgen to get structure factors.
"""

import os
import json
import numpy as np

from math import sin, cos, asin, pi, degrees
from pymatgen.analysis.diffraction.core import AbstractDiffractionPatternCalculator
from pymatgen.analysis.diffraction.core import DiffractionPattern
from pymatgen.analysis.diffraction.core import get_unique_families


with open(os.path.join(os.path.dirname(__file__),
                       "atomic_scattering_params.json")) as f:
    ATOMIC_SCATTERING_PARAMS = json.load(f)


def get_structure_factors(structure, min_d_spacing, scaled=True):
    """
    Calculates the diffraction pattern for a structure.

    Args:
        structure (Structure): Input structure
        min_d_spacing (float): Minimum d spacing
        scaled (bool): Whether to return scaled intensities. The maximum
            peak is set to a value of 100. Defaults to True. Use False if
            you need the absolute values to combine XRD plots.

    Returns:
        (XRDPattern)
    """
    latt = structure.lattice
    is_hex = latt.is_hexagonal()

    # Obtain crystallographic reciprocal lattice points within range
    recip_latt = latt.reciprocal_lattice_crystallographic
    recip_pts = recip_latt.get_points_in_sphere(
        [[0, 0, 0]], [0, 0, 0], 1 / min_d_spacing)

    # Create a flattened array of zs, coeffs, fcoords and occus. This is
    # used to perform vectorized computation of atomic scattering factors
    # later. Note that these are not necessarily the same size as the
    # structure as each partially occupied specie occupies its own
    # position in the flattened array.
    zs = []
    coeffs = []
    fcoords = []
    occus = []
    dwfactors = []

    for site in structure:
        for sp, occu in site.species.items():
            zs.append(sp.Z)
            try:
                c = ATOMIC_SCATTERING_PARAMS[sp.symbol]
            except KeyError:
                raise ValueError("Unable to calculate XRD pattern as "
                                 "there is no scattering coefficients for"
                                 " %s." % sp.symbol)
            coeffs.append(c)
            dwfactors.append(0)
            fcoords.append(site.frac_coords)
            occus.append(occu)

    zs = np.array(zs)
    coeffs = np.array(coeffs)
    fcoords = np.array(fcoords)
    occus = np.array(occus)
    dwfactors = np.array(dwfactors)
    planes = []
    intensities = []

    for hkl, g_hkl, ind, _ in sorted(
            recip_pts, key=lambda i: (i[1], -i[0][0], -i[0][1], -i[0][2])):
        # Force miller indices to be integers.
        hkl = [int(round(i)) for i in hkl]
        if g_hkl != 0:

            # s = sin(theta) / wavelength = 1 / 2d = |ghkl| / 2 (d =
            # 1/|ghkl|)
            s = g_hkl / 2

            # Store s^2 since we are using it a few times.
            s2 = s ** 2

            # Vectorized computation of g.r for all fractional coords and
            # hkl.
            g_dot_r = np.dot(fcoords, np.transpose([hkl])).T[0]

            # Highly vectorized computation of atomic scattering factors.
            # Equivalent non-vectorized code is::
            #
            #   for site in structure:
            #      el = site.specie
            #      coeff = ATOMIC_SCATTERING_PARAMS[el.symbol]
            #      fs = el.Z - 41.78214 * s2 * sum(
            #          [d[0] * exp(-d[1] * s2) for d in coeff])
            fs = zs - 41.78214 * s2 * np.sum(
                coeffs[:, :, 0] * np.exp(-coeffs[:, :, 1] * s2), axis=1)

            dw_correction = np.exp(-dwfactors * s2)

            # Structure factor = sum of atomic scattering factors (with
            # position factor exp(2j * pi * g.r and occupancies).
            # Vectorized computation.
            f_hkl = np.sum(fs * occus * np.exp(2j * pi * g_dot_r)
                           * dw_correction)

            # Intensity for hkl is modulus square of structure factor.
            i_hkl = (f_hkl * f_hkl.conjugate()).real

            if is_hex:
                # Use Miller-Bravais indices for hexagonal lattices.
                hkl = (hkl[0], hkl[1], - hkl[0] - hkl[1], hkl[2])

            planes.append(hkl)
            intensities.append(i_hkl)

    planes = np.array(planes)
    intensities = np.array(intensities)

    indices = np.argsort(intensities)[::-1]
    planes = planes[indices]
    intensities = intensities[indices]

    # TODO: merge planes which have the same intensity

    return planes, intensities
