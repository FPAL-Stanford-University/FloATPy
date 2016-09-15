#!/usr/bin/env python

from distutils.core import setup
import sys

# Python version
if sys.version_info[:2] < (2, 7):
    print('UFloAT requires Python 2.7 or newer')
    sys.exit(-1)

version = '1.0'

# Modules.
modules = ['ufloat.derivatives',
           'ufloat.readers']

setup(name='UFloAT',
      version=version,
      description='Postprocessing utilities for codes in UFPA lab',
      author='Unsteady Flow Physics and Aeroacoustics Laboratory of Stanford University',
      author_email='wongml@stanford.edu',
      packages=['ufloat'] + modules,
      )
