import psi4

def atoms_to_psi4(atoms, symmetry = 'c1'):
    """
    a function that takes in an atoms object and returns
    the psi4 input representation of that atoms object

    inputs:
        atoms (ASE atoms object):
            an ASE atoms object of the system you'd like
            to convert to psi4 format

    returns:
        psi4 (dict):
            the psi4 representation of the atoms object
    """
    result = ''
    for atom in atoms:
        temp = '{}\t{}\t{}\t{}\n'.format(atom.symbol, \
        atom.position[0],atom.position[1], \
        atom.position[2])
        result += temp
    result += '\t symmetry {}\n'.format(symmetry)
    result += 'units angstrom\n'
    return psi4.geometry(result)


def psi4_to_atoms(psi4):
    """
    a function to convert a psi4 representation of atoms
    to an ase atoms object. note that symmetry is not
    kept

    inputs:
         psi4 (dict)::
            the psi4 representation of the atoms object

    returns:
        atoms (ASE atoms object):
            the corresponding ase atoms object

    """
    from ase.atoms import Atoms
    
    atoms = Atoms(psi4['atoms'],
                  positions = psi4['coordinates'])
    
    
    
    

