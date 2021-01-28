import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from pymatgen.core.structure import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator


example_text = 'python generate_xrd.py --wavelength'

parser = argparse.ArgumentParser(epilog=example_text)
parser.add_argument('-w', '--wavelength', help='wavelength of the inciding radiation', type=float)
parser.add_argument('-m', '--minangle', help='minimum angle of the generated XRD spectra', type=float)
parser.add_argument('-M', '--maxangle', help='maximum angle of the generated XRD spectra', type=float)
parser.add_argument('-s', '--sigma', help='sigma value for gaussian broadening of peaks', type=float)

args = parser.parse_args()
cif_files = [f for f in os.listdir('.') if f.endswith('.cif')]

for f1 in tqdm(cif_files):
    structure = Structure.from_file(f1)
    calculator = XRDCalculator(wavelength=args.wavelength)
    pattern = calculator.get_pattern(structure, two_theta_range=(args.minangle, args.maxangle))

    # write spectrum file
    with open(f1[:-4] + '.txt', 'w') as f2:
        for angle, intensity in zip(pattern.x, pattern.y):
            f2.write(f'{angle:5.2f}    {intensity:6.2f}\n')

    # write png
    plt.figure(figsize=(16, 9))
    x = np.arange(args.minangle, args.maxangle, 0.01)
    y = np.zeros_like(x)

    for angle, intensity in zip(pattern.x, pattern.y):
        y += intensity * np.exp(-((x - angle) ** 2) / (2 * args.sigma ** 2))

    plt.plot(x, y)

    plt.savefig(f1[:-4] + '.png')
    plt.close()
