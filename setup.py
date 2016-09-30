#!/usr/bin/env python

from distutils.core import setup
import sys

# Python version
if sys.version_info[:2] < (2, 7):
    print('FloATPy requires Python 2.7 or newer')
    sys.exit(-1)

version = '1.0'

# Modules.
modules = ['floatpy.derivatives',
           'floatpy.readers']

setup(name='FloATPy',
      version=version,
      description='Postprocessing utilities for codes in FPAL',
      author='Flow Physics and Aeroacoustics Laboratory of Stanford University',
      author_email='wongml@stanford.edu',
      packages=['floatpy'] + modules,
      )
