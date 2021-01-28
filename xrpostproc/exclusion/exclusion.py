import os
import argparse
import numpy as np

from tqdm import tqdm
from copy import copy
from pymatgen import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator


example_text = '''example:

python generate_xrd.py -w 0.6199 -c 25 -l 25,31 -r 28,32 -m 5 -M 60'''

parser = argparse.ArgumentParser(epilog=example_text)
parser.add_argument('-w', '--wavelength', help='wavelength of the inciding radiation', type=float)
parser.add_argument('-m', '--minangle', help='minimum angle of the generated XRD spectra', type=float)
parser.add_argument('-M', '--maxangle', help='maximum angle of the generated XRD spectra', type=float)
parser.add_argument('-c', '--cutoff', help='cutoff in %% of the maximum peak', type=float)
parser.add_argument('-l', '--lextremes', help='comma separated left extremes of the exclusion regions', type=str)
parser.add_argument('-r', '--rextremes', help='comma separated right extremes of the exclusion regions', type=str)

args = parser.parse_args()

lextremes = args.lextremes.split(',')
rextremes = args.rextremes.split(',')
if len(lextremes) != len(rextremes):
    raise ValueError('Left and right extremes have different length.')

excl_regions = []
for l, r in zip(lextremes, rextremes):
    excl_regions.append((float(l), float(r)))

cif_files = [f for f in os.listdir('.') if f.endswith('.cif')]
to_remove = ''

for f in tqdm(cif_files):
    if to_remove in os.listdir('.'):
        os.remove(to_remove)

    structure = Structure.from_file(f)
    calculator = XRDCalculator(wavelength=args.wavelength)
    pattern = calculator.get_pattern(structure, two_theta_range=(args.minangle, args.maxangle))

    excl_peaks = []
    for region in excl_regions:
        excl_peaks.extend([peak for angle, peak in zip(pattern.x, pattern.y) if region[0] < angle < region[1]])

    excl_peaks = np.array(excl_peaks)
    if not np.alltrue(excl_peaks <= args.cutoff):
        to_remove = copy(f)

if to_remove in os.listdir('.'):
    os.remove(to_remove)

print('Done.')
