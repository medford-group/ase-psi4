"""
authors: Ben Comer (Georgia Tech), Xiangyun (Ray) Lei (Georgia Tech)

"""
from ase.calculators.calculator import Calculator, all_properties, all_changes
from ase.calculators.calculator import InputError, CalculationFailed, SCFError 
import numpy as np
from ase.units import Bohr, Hartree
import multiprocessing
import psi4


class Psi4(Calculator):
    """
    An ase calculator for the popular open source Q-chem code
    psi4. This is really rudimentary

    you can always use the in-built psi4 module through:
    calc.psi4
    """
    implemented_properties = ['energy', 'forces']
    
    default_parameters = {
                  "basis": "aug-cc-pvtz",
                  "num_threads": None,
                  "method": "hf",
                  "memory": None,
                  "D_CONVERGENCE":1e-12,
                  "E_CONVERGENCE":1e-12,
                  'DFT_BLOCK_MAX_POINTS': 500000,
                  'DFT_BLOCK_MIN_POINTS': 100000,
                  'charge': 0,
                  'multiplicity': 1,
                  'reference': None,
                  'symmetry':'c1',
#                  'DFT_SPHERICAL_POINTS': 302,
#                  'DFT_RADIAL_POINTS':    75,
                  "SAVE_JK": True, }
    def __init__(self, restart=None, ignore_bad_restart=False,
                 label='psi4-calc', atoms=None, command=None,
                 **kwargs):
        self.results = {}  # calculated properties (energy, forces, ...)

        # Use default parameters if they were not read from file:
        self.parameters = self.get_default_parameters()

        self.name = 'psi4'
        self.atoms = atoms
        self.set_label(label)
        self.set(**kwargs)
        self.psi4 = psi4
        # perform initial setup of psi4 python API
        self.set_psi4(atoms = atoms)

    def set_psi4(self, atoms = None):
        """
        This function sets the imported psi4 module to the settings the user 
        defines
        """
        # Input spin settings
        if self.parameters['reference'] is not None:
            self.psi4.set_options({'reference':
                                   self.parameters['reference']})
        # memory
        if self.parameters['memory'] is not None:
            self.psi4.set_memory(self.parameters['memory'])

        if self.parameters['num_threads'] == 'max':
            self.psi4.set_num_threads(multiprocessing.cpu_count())
        elif  type(self.parameters['num_threads']) == int:
            self.psi4.set_num_threads(kwargs['num_threads'])

        if atoms is None:
            if self.atoms is None:
                return None
            else:
                atoms = self.atoms 
        # Generate the atomic input
        result = ''
        for atom in atoms:
            temp = '{}\t{}\t{}\t{}\n'.format(atom.symbol, \
            atom.position[0],atom.position[1], \
            atom.position[2])
            result += temp
        result += '\t symmetry {}\n'.format(self.parameters['symmetry'])
        result += 'units angstrom\n'
        result += '{} {}'.format(self.parameters['charge'],
                                 self.parameters['multiplicity'])
        self.psi4.core.set_output_file(self.label + '.out', False)
        self.molecule = psi4.geometry(result)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes, symmetry = 'c1'):
        """Do the calculation.

        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these six: 'positions', 'numbers', 'cell',
            'pbc', 'initial_charges' and 'initial_magmoms'.

        Subclasses need to implement this, but can ignore properties
        and system_changes if they want.  Calculated properties should
        be inserted into results dictionary like shown in this dummy
        example::

            self.results = {'energy': 0.0,
                            'forces': np.zeros((len(atoms), 3)),
                            'stress': np.zeros(6),
                            'dipole': np.zeros(3),
                            'charges': np.zeros(len(atoms)),
                            'magmom': 0.0,
                            'magmoms': np.zeros(len(atoms))}

        """
        Calculator.calculate(self, atoms = atoms)
        if atoms == None:
            if self.atoms is None:
                raise InputError('An atoms object must be provided to perform a calculation')
            else:
                atoms = self.atoms 
        # this inputs all the settings into psi4
        self.set_psi4(atoms = atoms)

        # Set up the method
        method = self.parameters['method']
        basis = self.parameters['basis']

        # Do the calculations
        for item in properties:
            if item == 'energy':
                energy = self.psi4.energy('{}/{}'.format(method,basis), 
                                      molecule = self.molecule,)
                # convert to eV
                self.results['energy'] = energy * Hartree
            if item == 'forces':
                energy = self.psi4.energy('{}/{}'.format(method,basis),
                                      molecule = self.molecule,)
                self.results['energy'] = energy * Hartree
                grad = self.psi4.driver.gradient('{}/{}'.format(method,basis),
                                                 return_wfn=False,)
                # convert to eV/A
                # also note that the gradient is -1 * forces
                self.results['forces'] = -1 * np.array(grad) * Hartree * Bohr

