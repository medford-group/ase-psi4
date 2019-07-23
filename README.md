# ase-psi4

[![Build Status](https://travis-ci.org/medford-group/ase-psi4.svg?branch=master)](https://travis-ci.org/medford-group/ase-psi4)

A simple ase calculator for [psi4](http://www.psicode.org/psi4manual/master/index.html). You can load it and instantiate it the way you would any other ase calculator. The implemented methods are `get_forces` and `get_potential_energy`. You can also access the in built psi4 module with `calc.psi4`.

## Installation
You'll need psi4 to use this package, the easiest way to get that is with conda:

```
conda install -c psi4 psi4 
```

You can install this package quite simply with pip:

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
            xc = 'b3lyp',
            basis = '6-311g_d_p_')

atoms.set_calculator(calc)
print(atoms.get_potential_energy())
print(atoms.get_forces())
```

## Implemented Arguments

Here is a list of the things you can pass into the `Psi4` ase calculator with the corresponding default value:

basis: "aug-cc-pvtz"

xc: "hf"

charge: 0

multiplicity: 1

symmetry:'c1'

num\_threads: None


have fun!!
