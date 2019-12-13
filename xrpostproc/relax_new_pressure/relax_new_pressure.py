import os
from pymatgen.core.structure import Structure
from ase.calculators.vasp import Vasp
from ase.io import read

NEW_PRESSURE = 50

cif_files = [f for f in os.listdir('.') if '.cif' in f]

calc = Vasp(
    xc='PBE',
    setups='recommended',
    istart=0,
    prec='Accurate',
    symprec=1e-6,
    encut=700,
    ediff=1e-5,
    ibrion=2,
    isif=3,
    nsw=200,
    ismear=-5,
    sigma=0.005,
    potim=0.05,
    lcharg=False,
    lwave=False,
    ediffg=1e-3,
    kspacing=0.25,
    pstress=NEW_PRESSURE * 10
)

for f in cif_files:
    atoms = read(f, format='cif')
    atoms.set_calculator(calc)

    # the next line triggers a VASP relaxation
    energy = atoms.get_potential_energy()

    filename = f.strip('.cif') + '_to' + str(NEW_PRESSURE) + 'GPa' + '.cif'
    structure = Structure(atoms.cell.array, atoms.numbers, atoms.positions, coords_are_cartesian=True)
    try:
        structure.to(filename=filename, symprec=0.2)
    except TypeError:
        print('WARNING: {} failed to be created.'.format(filename))
