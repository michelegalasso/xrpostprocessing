"""
@file:      scxrd_screening.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      2 November 2020
@brief:     Script which uses the class SpectrumAnalyzer.
"""

import iotbx.cif
import cctbx.crystal

from pymatgen.core.structure import Structure
from iotbx.reflection_file_reader import any_reflection_file


# import numpy as np
# hkl_data = np.loadtxt('test.hkl')
# np.savetxt('strong_reflections.hkl', hkl_data[hkl_data[:, 3] > 1000], fmt=['%4.0f', '%3.0f', '%3.0f', '%7.2f', '%7.2f'])


def compute_r_factors(fobs, fmodel, flags):
    fmodel, fobs = fmodel.common_sets(other=fobs)
    fmodel, flags = fmodel.common_sets(other=flags)
    fc_work = fmodel.select(~(flags.data()))
    fo_work = fobs.select(~(flags.data()))
    fc_test = fmodel.select(flags.data())
    fo_test = fobs.select(flags.data())
    r_work = fo_work.r1_factor(fc_work)
    r_free = fo_test.r1_factor(fc_test)
    print("r_work = %.4f" % r_work)
    print("r_free = %.4f" % r_free)
    print("")
    flags.setup_binner(n_bins=20)
    fo_work.use_binning_of(flags)
    fc_work.use_binner_of(fo_work)
    fo_test.use_binning_of(fo_work)
    fc_test.use_binning_of(fo_work)
    for i_bin in fo_work.binner().range_all():
        sel_work = fo_work.binner().selection(i_bin)
        sel_test = fo_test.binner().selection(i_bin)
        fo_work_bin = fo_work.select(sel_work)
        fc_work_bin = fc_work.select(sel_work)
        fo_test_bin = fo_test.select(sel_test)
        fc_test_bin = fc_test.select(sel_test)
        if fc_test_bin.size() == 0 : continue
        r_work_bin = fo_work_bin.r1_factor(other=fc_work_bin, assume_index_matching=True)
        r_free_bin = fo_test_bin.r1_factor(other=fc_test_bin, assume_index_matching=True)
        cc_work_bin = fo_work_bin.correlation(fc_work_bin).coefficient()
        cc_free_bin = fo_test_bin.correlation(fc_test_bin).coefficient()
        legend = flags.binner().bin_legend(i_bin, show_counts=False)
        print("%s  %8d %8d  %.4f %.4f  %.3f %.3f" % (legend, fo_work_bin.size(), fo_test_bin.size(), r_work_bin,
                                                     r_free_bin, cc_work_bin, cc_free_bin))


# read input hkl file
hkl_in = any_reflection_file(file_name='strong_reflections.hkl=intensities')

# give input unit cell and crystal symmetry
crystal_symmetry = cctbx.crystal.symmetry(unit_cell=(4.7877, 4.9480, 6.9151, 90, 90, 90), space_group_symbol='P1')

# get experimental Miller arrays
miller_arrays = hkl_in.as_miller_arrays(crystal_symmetry)
i_obs = miller_arrays[0]

# convert to amplitudes
f_obs = i_obs.f_sq_as_f()

# use any_file(structure.cif) and then get f_model somehow, maybe after refinement.
# see webpage "Tour of the cctbx"

structure = Structure.from_file('POSCAR')
cif_string = structure.to(fmt='cif')
reader = iotbx.cif.reader(input_string=cif_string)
model = reader.build_crystal_structures()['MgSiO3']
f_model = model.structure_factors(d_min=2).f_calc().as_amplitude_array()

print(f_model.show_summary())
print('\n')
print(f_obs.show_summary())
