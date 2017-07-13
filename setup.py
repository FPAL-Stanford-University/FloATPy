#!/usr/bin/env python

from numpy.distutils.core import setup, Extension

from setup_fortran import BuildFortranObjects

import sys

# Python version
if sys.version_info[:2] < (2, 7):
    print('FloATPy requires Python 2.7 or newer')
    sys.exit(-1)

version = '1.0'

# Build the fortran extension modules.
objects = BuildFortranObjects(['floatpy/derivatives/compact/kind_parameters.F90',
                               'floatpy/derivatives/compact/constants.F90',
                               'floatpy/derivatives/compact/cd06.F90'])

ext_compact_6th_order = Extension('_compact_sixth_order',
                                    sources = ['floatpy/derivatives/compact/f90wrap_kind_parameters.f90',
                                               'floatpy/derivatives/compact/f90wrap_constants.f90',
                                               'floatpy/derivatives/compact/f90wrap_cd06.f90'],
                                    extra_objects = objects,
                                    f2py_options = [])

# Modules.
modules = ['floatpy.derivatives',
           'floatpy.derivatives.compact',
           'floatpy.upsampling',
           'floatpy.readers']

setup(name='FloATPy',
      version=version,
      description='Postprocessing utilities for codes in FPAL',
      author='Flow Physics and Aeroacoustics Laboratory of Stanford University',
      author_email='wongml@stanford.edu',
      ext_modules = [ext_compact_6th_order],
      packages=['floatpy'] + modules,
      )
