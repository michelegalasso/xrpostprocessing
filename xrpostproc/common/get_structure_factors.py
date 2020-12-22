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

from math import sin, cos, asin, pi, degrees, radians
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.diffraction.core import AbstractDiffractionPatternCalculator
from pymatgen.analysis.diffraction.core import DiffractionPattern
from pymatgen.analysis.diffraction.core import get_unique_families


with open(os.path.join(os.path.dirname(__file__),
                       "atomic_scattering_params.json")) as f:
    ATOMIC_SCATTERING_PARAMS = json.load(f)


def get_structure_factors(calculator, structure, scaled=True, two_theta_range=(0, 90)):
    """
    Calculates the diffraction pattern for a structure.

    Args:
        structure (Structure): Input structure
        scaled (bool): Whether to return scaled intensities. The maximum
            peak is set to a value of 100. Defaults to True. Use False if
            you need the absolute values to combine XRD plots.
        two_theta_range ([float of length 2]): Tuple for range of
            two_thetas to calculate in degrees. Defaults to (0, 90). Set to
            None if you want all diffracted beams within the limiting
            sphere of radius 2 / wavelength.

    Returns:
        (XRDPattern)
    """
    if calculator.symprec:
        finder = SpacegroupAnalyzer(structure, symprec=calculator.symprec)
        structure = finder.get_refined_structure()

    wavelength = calculator.wavelength
    latt = structure.lattice
    is_hex = latt.is_hexagonal()

    # Obtained from Bragg condition. Note that reciprocal lattice
    # vector length is 1 / d_hkl.
    min_r, max_r = (0, 2 / wavelength) if two_theta_range is None else \
        [2 * sin(radians(t / 2)) / wavelength for t in two_theta_range]

    # Obtain crystallographic reciprocal lattice points within range
    recip_latt = latt.reciprocal_lattice_crystallographic
    recip_pts = recip_latt.get_points_in_sphere(
        [[0, 0, 0]], [0, 0, 0], max_r)
    if min_r:
        recip_pts = [pt for pt in recip_pts if pt[1] >= min_r]

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
            dwfactors.append(calculator.debye_waller_factors.get(sp.symbol, 0))
            fcoords.append(site.frac_coords)
            occus.append(occu)

    zs = np.array(zs)
    coeffs = np.array(coeffs)
    fcoords = np.array(fcoords)
    occus = np.array(occus)
    dwfactors = np.array(dwfactors)
    peaks = {}
    two_thetas = []

    for hkl, g_hkl, ind, _ in sorted(
            recip_pts, key=lambda i: (i[1], -i[0][0], -i[0][1], -i[0][2])):
        # Force miller indices to be integers.
        hkl = [int(round(i)) for i in hkl]
        if g_hkl != 0:

            d_hkl = 1 / g_hkl

            # Bragg condition
            theta = asin(wavelength * g_hkl / 2)

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

            # Lorentz polarization correction for hkl
            lorentz_factor = (1 + cos(2 * theta) ** 2) / \
                             (sin(theta) ** 2 * cos(theta))

            # Intensity for hkl is modulus square of structure factor.
            i_hkl = (f_hkl * f_hkl.conjugate()).real

            two_theta = degrees(2 * theta)

            if is_hex:
                # Use Miller-Bravais indices for hexagonal lattices.
                hkl = (hkl[0], hkl[1], - hkl[0] - hkl[1], hkl[2])
            # Deal with floating point precision issues.
            ind = np.where(np.abs(np.subtract(two_thetas, two_theta)) <
                           AbstractDiffractionPatternCalculator.TWO_THETA_TOL)
            if len(ind[0]) > 0:
                peaks[two_thetas[ind[0][0]]][0] += i_hkl * lorentz_factor
                peaks[two_thetas[ind[0][0]]][1].append(tuple(hkl))
            else:
                peaks[two_theta] = [i_hkl * lorentz_factor, [tuple(hkl)],
                                    d_hkl]
                two_thetas.append(two_theta)

    # Scale intensities so that the max intensity is 100.
    max_intensity = max([v[0] for v in peaks.values()])
    x = []
    y = []
    hkls = []
    d_hkls = []
    for k in sorted(peaks.keys()):
        v = peaks[k]
        fam = get_unique_families(v[1])
        if v[0] / max_intensity * 100 > AbstractDiffractionPatternCalculator.SCALED_INTENSITY_TOL:
            x.append(k)
            y.append(v[0])
            hkls.append([{"hkl": hkl, "multiplicity": mult}
                         for hkl, mult in fam.items()])
            d_hkls.append(v[2])
    xrd = DiffractionPattern(x, y, hkls, d_hkls)
    if scaled:
        xrd.normalize(mode="max", value=100)
    return xrd
