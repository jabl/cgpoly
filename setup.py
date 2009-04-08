#!/usr/bin/python

from distutils.core import setup
from glob import glob

# Scripts whose names end in a-z or 1-9 (avoids emacs backup files)
scripts = glob('scripts/*[a-z,1-9]')

setup(name='cgpoly',
      version='0.1',
      description='Package for doing coarse grained polymer models with GROMACS',
      author='Janne Blomqvist',
      author_email='Janne.Blomqvist@tkk.fi',
      url='http://tfy.tkk.fi/~job/',
      packages=['cgpoly', 'cgpoly.test'],
      scripts = scripts)

