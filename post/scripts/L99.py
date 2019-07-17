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

xdir = 0
zdir = 2

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
   
 
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:" 
        print "Computes the momentum thickness and decorrelation lengths" 
        print "  python {} <prefix> [tID_start(default=0)] ".format(sys.argv[0])
        sys.exit()
    start_index = 0;
    if len(sys.argv) > 2:
        start_index = int(sys.argv[2])
    filename_prefix = sys.argv[1]
    outputfile = filename_prefix + "L99.dat"
    
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
    tid_list = steps[start_index:]

    # Set up compact derivative object w/ 10th order schemes
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)
    
    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=True)
    du = inp.du
    if rank==0: print("\tdu = {}".format(inp.du))
    utop = inp.du/2.*0.99# this is wrong, ubot is neg
    ubot = inp.du/2.*0.01

    # Preallocate for means and derivatives
    Nsteps = np.size(tid_list)
    Nx,Ny,Nz = reader.domain_size
    time = np.zeros([Nsteps])
    L99 = np.zeros([Nsteps])
    
    # Compute stats at each step:
    i = 0
    if rank==0: print("Time \t 99%")
    for tid in tid_list:
        reader.step = tid
        time[i] = reader.time

        # Mean velocity 
        u = reader.readData( ('u') )
        u = np.squeeze(np.array(u))
        ubar = stats.reynolds_average(avg, u)
         
        # Compute momentum thickness
        i1 = np.argmin(abs(ubar[0,:Ny/2,0]-ubot))
        i2 = np.argmin(abs(ubar[0,Ny/2:,0]-utop)) + Ny/2
        L99[i] = (y[0,i2,0] - y[0,i1,0])/2.
        if rank==0: print("{} \t {}".format(time[i],L99[i]))
        i+=1

    # Write to file 
    if rank==0:
        array2save = np.empty((Nsteps,2))
        array2save[:,0] = time
        array2save[:,1] = L99 
        np.savetxt(outputfile,array2save,delimiter=' ')
        print("Done writing to {}".format(outputfile))
