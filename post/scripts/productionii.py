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
        print "Usage:" 
        print "Computes utilde" 
        print "  python {} <prefix> [tID_list (csv)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    if len(sys.argv) > 2:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
    else: tID_list = None
    
    dirname = os.path.dirname(filename_prefix)
    if 'Mc04' in dirname:
        dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/temporal/')[-1]
        dir_out = '/home/kmatsuno/ShearLayerData/temporal/' + dir_out + '/'
    else:
        dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/mira/')[-1]
        dir_out = '/home/kmatsuno/ShearLayerData/production/' + dir_out + '/'
    
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

    # Set up compact derivative object w/ 10th order schemes
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)
    
    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=(rank==0))
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    du = inp.du
    if rank==0: print("\tdu = {}".format(inp.du))
    nx,ny,nz = reader._full_chunk_size
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)
    
    for tID in tID_list:
        reader.step = tID
        r,u,v = reader.readData( ('rho','u','v') )
        if np.shape(r)[1] < Ny:
            r = transpose2y(settings,r)
            u = transpose2y(settings,u)
            v = transpose2y(settings,v)
        utilde = stats.favre_average(avg,r,u)
        vtilde = stats.favre_average(avg,r,v)
        upp = u - utilde
        vpp = v - vtilde
        
        R12 = stats.reynolds_average(avg,r*upp*vpp)
        upp,vpp = None,None

        qx,dudy,qz = der.gradient(u, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        dudy_tilde = stats.favre_average(avg,r,dudy)
        qx,qz = None,None
        
        P = np.squeeze(R12*dudy_tilde)

        if rank==0: 
            outputfile = dir_out+"/shearlayer_production_%04d.dat"%tID
            np.savetxt(outputfile,P,delimiter=' ')
            print("Done writing to {}".format(outputfile))
