from ase_psi4.ase_psi4 import Psi4
from ase.build import molecule
import numpy as np

atoms = molecule('H2O')
atoms[0].position += np.array([0,0,0.1])

calc = Psi4(atoms = atoms,method = 'b3lyp')

atoms.set_calculator(calc)

print(atoms.get_potential_energy())
