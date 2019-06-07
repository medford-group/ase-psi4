# ase-psi4

A simple ase calculator for Psi4. You can load it and instantiate it the way you would any other ase calculator. The implemented methods are `get_forces` and `get_potential_energy`. You can also access the in built psi4 module with `calc.psi4`. This was thrown together rather quickly, so I make no guarantees about it's functionality

## Installation

You can install it simply with pip:

```
pip install git+https://github.com/medford-group/ase-psi4
```

or alternatively you can clone it and use `setup.py`

```
git clone https://github.com/medford-group/ase-psi4
python setup.py install
```

## Example

```
from ase_psi4.ase_psi4 import Psi4
from ase.build import molecule
import numpy as np

atoms = molecule('H2O')

calc = Psi4(atoms = atoms,
            method = 'b3lyp',
            basis = '6-311g_d_p_')

atoms.set_calculator(calc)
atoms.get_potential_energy()
```

have fun!!
