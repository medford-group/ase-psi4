"""
author: Ben Comer (Georgia Tech)
"""
from ase.calculators.calculator import Calculator, all_properties, all_changes
#from ase.utils.timing import Timer
#from .utilities import psi4_to_atoms
import numpy as np
from ase.units import Bohr, Hartree
import multiprocessing
#import threading
import psi4


#all_properties = ['energy', 'forces', 'stress', 'dipole',
#                  'charges', 'magmom', 'magmoms', 'free_energy']


#all_changes = ['positions', 'numbers', 'cell', 'pbc',
#               'initial_charges', 'initial_magmoms']

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
                  "method": "hf",
                  "D_CONVERGENCE":1e-12,
                  "E_CONVERGENCE":1e-12,
                  'DFT_BLOCK_MAX_POINTS': 500000,
                  'DFT_BLOCK_MIN_POINTS': 100000,
                  'MAXITER': 500,
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
        self.parameters = None  # calculational parameters

        if self.parameters is None:
            # Use default parameters if they were not read from file:
            self.parameters = self.get_default_parameters()

        self.name = 'psi4'
        self.atoms = atoms
        self.set_label(label)
        self.set(**kwargs)
        self.psi4 = psi4
        #self.psi4.set_num_threads(threading.active_count())
        self.psi4.set_num_threads(multiprocessing.cpu_count()-1)

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
        Calculator.calculate(self,atoms=atoms)
        if atoms == None:
            atoms = self.atoms 
        # Generate the atomic input
        result = ''
        for atom in atoms:
            temp = '{}\t{}\t{}\t{}\n'.format(atom.symbol, \
            atom.position[0],atom.position[1], \
            atom.position[2])
            result += temp
        result += '\t symmetry {}\n'.format(symmetry)
        result += 'units angstrom\n'
        result += '{} {}'.format(self.parameters['charge'],
                                 self.parameters['multiplicity'])
        self.psi4.core.set_output_file(self.label + '.out', False)
        molecule = psi4.geometry(result)
        #molecule.charge(1)

        # Set up the method
        method = self.parameters['method']
        basis = self.parameters['basis']
        if self.parameters['reference'] is not None:
            psi4.set_options({'reference': 
                              self.parameters['reference']}) 
        # Do the calculations
        for item in properties:
            if item == 'energy':
                energy = self.psi4.energy('{}/{}'.format(method,basis), 
                                      molecule = molecule)
                # convert to eV
                self.results['energy'] = energy * Hartree
            if item == 'forces':
                energy = self.psi4.energy('{}/{}'.format(method,basis),
                                      molecule = molecule, MAXITER = 1000)
                self.results['energy'] = energy * Hartree
                grad = self.psi4.driver.gradient('{}/{}'.format(method,basis),
                                                 return_wfn=False)
                # convert to eV/A
                # also note that the gradient is -1 * forces
                self.results['forces'] = -1 * np.array(grad) * Hartree * Bohr

