#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
from common import *
import statistics as stats
import get_namelist as nml
from SettingLib import NumSetting
from decorr_lscale_y import transpose2y

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_start(default=0)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    if len(sys.argv) > 2:
        tID_list = map(int, sys.argv[2].strip('[]').split(','))
    else:
        start_index = 0

    
    periodic_dimensions = (True,False,True)
    x_bc = (0,0)
    y_bc = (0,0)
    z_bc = (0,0)

    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    procs = comm.Get_size()

    # Set up the serial Miranda reader
    # Set up the parallel reader
    # Set up the reduction object
    serial_reader = por.PadeopsReader(filename_prefix, 
            periodic_dimensions=periodic_dimensions)
    reader = pdr.ParallelDataReader(comm, serial_reader)
    avg = red.Reduction(reader.grid_partition, periodic_dimensions)
    steps = sorted(reader.steps)
    if tID_list is None: tID_list = steps

    # Set up the derivative object
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)

    # setup the inputs object, get grid info
    dirname = os.path.dirname(filename_prefix)
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    Ny = int(Ny)
    inp = nml.inputs(dirname,verbose=(rank==0))
    du = inp.du
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)
   
    # Compute stats at each step:
    for tID in tID_list: 
        reader.step =tID
        
        # density and streamwise vel, means
        T,p = reader.readData( ('T','p') )
        if procs>1: 
            T = transpose2y(settings,T)
            p = transpose2y(settings,p)

        T0 = 1.0
        p0 = 1.0
        ds_cp = np.log(T/T0) - (inp.gam-1)/inp.gam * np.log(p/p0)
        dat = stats.reynolds_average(avg,ds_cp)
        # p/rho^gam = constant? check rms error
        """
        const = p/r**inp.gam
        mean = stats.reynolds_average(avg,const)
        fluct = const-mean
        err = np.sqrt(fluct**2)
        dat = stats.reynolds_average(avg,err)/mean
        """

        if rank==0: 
            dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/mira/')[-1]
            dir_out = '/home/kmatsuno/ShearLayerData/production/' + dir_out + '/'
            outputfile = dir_out+"isentropic_check_%04d.dat"%tID
            print("Writing to {}".format(outputfile))
            np.savetxt(outputfile,np.squeeze(dat),delimiter=' ')
            print('Done')
