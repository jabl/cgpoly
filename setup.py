#!/usr/bin/python

from distutils.core import setup
from glob import glob

version = '0.2'

# Scripts whose names end in a-z or 1-9 (avoids emacs backup files)
scripts = glob('scripts/*[a-z,0-9]')
data = glob('data/*[a-z,0-9]')
#datadir = 'share/cgpoly-' + version
datadir = 'share/cgpoly'

setup(name='cgpoly',
      version=version,
      description='Package for doing coarse grained polymer models with GROMACS',
      author='Janne Blomqvist',
      author_email='Janne.Blomqvist@tkk.fi',
      url='http://tfy.tkk.fi/~job/',
      #packages=['cgpoly', 'cgpoly.test'],
      packages=['cgpoly'],
      data_files = [(datadir, data)],
      scripts = scripts)

