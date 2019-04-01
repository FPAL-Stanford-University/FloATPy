from mpi4py import MPI
import numpy as np
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats

if __name__ == '__main__':
   
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_start(default=0)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    start_index = 0
    if len(sys.argv) > 2:
        start_index = int(sys.argv[2])
    outputfile  = filename_prefix + "TKE_growth.dat"

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

    # Set up the derivative object
    x, y, z = reader.readCoordinates()
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]

    # Compute stats at each step:
    Nsteps = np.size(steps[start_index:])-1
    time = np.empty(Nsteps)
    tke  = np.empty(Nsteps)
    i = 0
    if rank == 0:
        print("Time \t TKE integrated") 
    for step in steps[start_index:-1]:
        reader.step = step
        time[i] = reader.time
        rho, u, v, w = reader.readData( ('rho', 'u', 'v', 'w') )
        tke[i] = stats.TKE(rho, u, v, w, reader.grid_partition, avg, dx, dy, dz)
        
        if rank == 0:
            print("{} \t {}".format(reader.time,tke[i]))
        i = i+1

    # Write to file 
    if rank==0:
        array2save = np.empty((Nsteps,2))
        array2save[:,0] = time
        array2save[:,1] = tke 
        np.savetxt(outputfile,array2save,delimiter=' ')
        print("Done writing to {}".format(outputfile))
