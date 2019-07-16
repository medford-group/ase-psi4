from ase_psi4.ase_psi4 import Psi4
#import Psi4
from ase.build import molecule

def h2dft(name):
    Calculator = Psi4
    par = {}
    calc = Calculator(label=name, xc='LDA', **par)
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    e2 = h2.get_potential_energy()
    calc.set(xc='PBE')
    e2pbe = h2.get_potential_energy()
    h1 = h2.copy()
    del h1[1]
    h1.set_initial_magnetic_moments([1])
    h1.calc = calc
    e1pbe = h1.get_potential_energy()
    calc.set(xc='LDA')
    e1 = h1.get_potential_energy()
    try:
        m1 = h1.get_magnetic_moment()
    except NotImplementedError:
        pass
    else:
        print(m1)
    print(2 * e1 - e2)
    print(2 * e1pbe - e2pbe)
    print(e1, e2, e1pbe, e2pbe)
    calc = Calculator(name)
    print(calc.parameters, calc.results, calc.atoms)
    assert not calc.calculation_required(h1, ['energy'])
    h1 = calc.get_atoms()
    print(h1.get_potential_energy())
    label = 'dir/' + name + '-h1'
    calc = Calculator(label=label, atoms=h1, xc='LDA', **par)
    print(h1.get_potential_energy())
    print(Calculator.read_atoms(label).get_potential_energy())

h2dft('psi4-calc')
