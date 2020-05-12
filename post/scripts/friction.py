#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml
from SettingLib import NumSetting

debug = False
def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz

def mpi_max(comm,value):
    local_max = np.amax(value)
    global_max = comm.reduce(local_max, op=MPI.MAX, root=0)
    return global_max
    
# Transpose from 3d to y
def transpose2y(numSet,q3d): 
    bufy = np.empty(numSet.chunk_y_size, 
            dtype=np.dtype(np.float64, align=True), order='F')
    numSet.grid_partition.transpose_3d_to_y(q3d,bufy)
    q3d = None
    return bufy

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_start(default=0)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    start_index = 0
    if len(sys.argv) > 2:
        start_index = int(sys.argv[2])
    
    periodic_dimensions = (True,False,True)
    x_bc = (0,0)
    y_bc = (1,1)
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

    # Set up the derivative object
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)

    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,rank==0)
    if rank==0: print("du = {}".format(inp.du))
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    yplot = np.linspace(-Ly/2,Ly/2,int(Ny))
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)

    for i,tID in enumerate(steps): 
        reader.step = tID
        time = reader.time
        
        # friction vel 
        if rank==0: print('Computing gradient')
        q = reader.readData(('u'))
        dudx,dudy,dudz = der.gradient(q[0], x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        dudy_bar = stats.reynolds_average(avg,dudy)
        q,dudx,dudy,dudz = None,None,None,None


        if rank==0: print('Reading data')
        r,mu = reader.readData(('rho','mu'))
        if rank==0: print('Computing friction vel')
        tmp = transpose2y(settings,np.sqrt(mu*abs(dudy_bar)/r))
        utau = stats.reynolds_average(avg,tmp)
        utau = np.squeeze(utau) 
        # Write to file 
        if rank==0:
            outputfile = filename_prefix + "utau_%04d.dat"%tID
            np.savetxt(outputfile,utau,delimiter=' ')
            print("Done writing to {}".format(outputfile))
