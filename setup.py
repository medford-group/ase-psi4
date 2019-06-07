#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

setup(name='ase_psi4',
      version='0.1',
      description='Python Wrapper for the Psi4 Q-Chem code',
      author='Ben Comer',
      author_email='ben.comer@gatech.edu',
      url='https://github.com/medford-group/ase-psi4',
      packages=find_packages(),
      install_requires=['spglib', 'numpy','ase','scipy'],
     )
