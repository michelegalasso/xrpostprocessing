import os
import sys
import numpy as np

from tqdm import tqdm
from copy import copy
from pymatgen import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator

if len(sys.argv) < 4:
    print('ERROR: wrong number of arguments.')
    print('Usage: python3 excl.py [wavelength] [cutoff in % max] [exclusion regions in degrees]')
    print('Example: python3 exclusion.py 0.6199 25 25-28 31-32')
    sys.exit()

wl = float(sys.argv[1])
cutoff = float(sys.argv[2])
excl_regions = []
for i in range(3, len(sys.argv)):
    excl_regions.append(tuple(float(x) for x in sys.argv[i].split('-')))

cif_files = [f for f in os.listdir('.') if '.cif' in f]
to_remove = ''

for f in tqdm(cif_files):
    if to_remove in os.listdir('.'):
        os.remove(to_remove)

    structure = Structure.from_file(f)
    calculator = XRDCalculator(wavelength=wl)
    pattern = calculator.get_pattern(structure)

    excl_peaks = []
    for region in excl_regions:
        excl_peaks.extend([peak for angle, peak in zip(pattern.x, pattern.y) if region[0] < angle < region[1]])

    excl_peaks = np.array(excl_peaks)
    if not np.alltrue(excl_peaks <= cutoff):
        to_remove = copy(f)

if to_remove in os.listdir('.'):
    os.remove(to_remove)

print('Done.')
