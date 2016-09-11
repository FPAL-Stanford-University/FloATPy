#!/usr/bin/env python

from distutils.core import setup
import sys

# Python version
if sys.version_info[:2] < (2, 7):
    print('PyFR requires Python 2.7 or newer')
    sys.exit(-1)

version = '1.0'

# Modules.
modules = ['fpalpost.derivatives',
           'fpalpost.readers']

setup(name='FPAL-Post',
      version=version,
      description='Postprocessing utilities for FPAL\'s code',
      author='Flow Physics and Aeroacoustics Laboratory of Stanford University',
      author_email='wongml@stanford.edu',
      packages=['fpalpost'] + modules,
      )
