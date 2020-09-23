from mpi4py import MPI
import numpy as np
import numpy.linalg as la
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml
from SettingLib import NumSetting
from PoissonSol import * 
from common import *

xdir = 0
zdir = 2

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_start(default=0)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    start_index = 0
    if len(sys.argv) > 2:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
    else: tID_list = None

    dirname = os.path.dirname(filename_prefix)
    #dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/mira/')[-1]
    #dir_out = '/home/kmatsuno/ShearLayerData/production/' + dir_out + '/'
    dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/temporal/')[-1]
    dir_out = '/home/kmatsuno/ShearLayerData/temporal/' + dir_out + '/'

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

    x, y, z = reader.readCoordinates()
    Nx,Ny,Nz = reader.domain_size
    dx,dy,dz = grid_res(x,y,z)
    
    # setup the inputs object, get grid info
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    Ny = int(Ny)
    if rank==0: verbose=True
    else: verbose=False
    inp = nml.inputs(dirname,verbose)
    du = inp.du
    if rank==0: print("du = {}".format(inp.du))
    
    # Get the grid partition information
    nx,ny,nz = reader._full_chunk_size
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)

    # Compute stats at each step:
    for tID in tID_list: 
        reader.step = tID

        q = reader.readData( ('T') )
        T = transpose2y(settings,q[0])

        # Speed of sound, xz average
        c3D = np.sqrt(inp.gam*T)
        cbar = stats.reynolds_average(avg, c3D)
        cbar = np.squeeze(cbar) 
        
        if rank==0: 
            outputfile = dir_out+"/shearlayer_cbar_%04d.dat"%tID
            np.savetxt(outputfile,cbar,delimiter=' ')
            print("Done writing to {}".format(outputfile))
