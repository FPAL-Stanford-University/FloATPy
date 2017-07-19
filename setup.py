#!/usr/bin/env python

from numpy.distutils.core import setup, Extension

from setup_fortran import BuildFortranObjects

import sys
import os
import subprocess

try:
    MPIFC = os.environ['MPIFC']
except:
    MPIFC = subprocess.check_output(['which','mpif90']).split()[0]
print "Using MPI Fortran compiler ", MPIFC

MPI_INCLUDE_DIRS = subprocess.check_output([MPIFC,'--showme:incdirs']).split()
MPI_LINK_DIRS = subprocess.check_output([MPIFC,'--showme:libdirs']).split()

f2py_mpi_include = '--include-paths '
for incdir in MPI_INCLUDE_DIRS:
    f2py_mpi_include += incdir + ' '
print "f2py_mpi_include = ", f2py_mpi_include

f2py_mpi_link = subprocess.check_output([MPIFC,'--showme:link']).strip()
print "f2py_mpi_link = ", f2py_mpi_link


# Python version
if sys.version_info[:2] < (2, 7):
    print('FloATPy requires Python 2.7 or newer')
    sys.exit(-1)

version = '1.0'

# Build the fortran extension modules.
obj_compact_6th_order = BuildFortranObjects(['floatpy/derivatives/compact/kind_parameters.F90',
                                             'floatpy/derivatives/compact/constants.F90',
                                             'floatpy/derivatives/compact/cd06.F90'])

ext_compact_6th_order = Extension('_compact_sixth_order',
                                    sources = ['floatpy/derivatives/compact/f90wrap_kind_parameters.f90',
                                               'floatpy/derivatives/compact/f90wrap_constants.f90',
                                               'floatpy/derivatives/compact/f90wrap_cd06.f90'],
                                    extra_objects = obj_compact_6th_order,
                                    f2py_options = [])

obj_pyt3d = BuildFortranObjects(['floatpy/parallel/pyt3d/kind_parameters.F90',
                                 'floatpy/parallel/pyt3d/constants.F90',
                                 'floatpy/parallel/pyt3d/exits.F90',
                                 'floatpy/parallel/pyt3d/reductions.F90',
                                 'floatpy/parallel/pyt3d/t3dMod.F90'], compiler=MPIFC)

ext_pyt3d = Extension('_pyt3d',
                      sources = ['floatpy/parallel/pyt3d/f90wrap_t3dMod.f90'],
                      extra_objects = obj_compact_6th_order,
                      f2py_options = [f2py_mpi_include, f2py_mpi_link])

# Modules.
modules = ['floatpy.parallel',
           'floatpy.parallel.pyt3d',
           'floatpy.derivatives',
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
