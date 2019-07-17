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

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    step = int(sys.argv[2])

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
    #steps = sorted(reader.steps)

    x, y, z = reader.readCoordinates()
    Nx,Ny,Nz = reader.domain_size
    dx,dy,dz = grid_res(x,y,z)
    
    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    if rank==0: verbosity=True
    else: verbosity=False
    inp = nml.inputs(dirname,verbose=verbosity)
    du = inp.du
    if rank==0: print("\tdu = {}".format(inp.du))

    # Speed of sound, xz average

    reader.step = step
    T = reader.readData( 'T' )
    T = np.array(T)
    c3D = np.sqrt(inp.gam*T)
    cbar = stats.reynolds_average(avg, c3D)
    cmin = np.amin(c3D)
    cmax = np.amax(c3D)
    norm = np.amax( (cmax-cmin)/cbar )

    print("cmin: {} \t cmax:{} \t norm max:{}".format(cmin,cmax,norm))
