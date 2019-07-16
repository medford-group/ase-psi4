from ase_psi4.ase_psi4 import Psi4
from ase.optimize import BFGSLineSearch
from ase.build import molecule
import numpy as np

atoms = molecule('H2O')
atoms[0].position += np.array([0,0,0.1])
"""
calc = Psi4(atoms = atoms,xc = 'ccsd(t)', basis = '6-31g**')


atoms.set_calculator(calc)

print(calc.parameters)

# test energy call
print(atoms.get_potential_energy())

calc = Psi4('psi4-calc')
print(calc.directory)
print(calc.parameters)
# test interface with ASE optimizers
relax = BFGSLineSearch(atoms)
relax.run(fmax = 0.05)

# test direct interaction with psi4
energy = calc.psi4.energy('b3lyp/cc-pvdz',
                   molecule = calc.molecule,
                   return_wfn = False)
"""
# test more obscure settings
calc = Psi4(atoms = atoms,
            method = 'b3lyp',
            basis = '6-31g**',
            num_threads = 'max',
            memory = '500MB',
            reference = 'uhf')

#print(calc.atoms)
calc.calculate()


calc.psi4.frequency('scf/cc-pvdz', molecule=calc.molecule, 
                    return_wfn=True, dertype=1)
