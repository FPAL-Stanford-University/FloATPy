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

debug = False
def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_start(default=0)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    start_index = 0
    if len(sys.argv) > 2:
        start_index = int(sys.argv[2])
    outputfile  = filename_prefix + "anisotropy.dat"
    
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
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)

    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=True)
    if rank==0: print("\tdu = {}".format(inp.du))

    # Compute stats at each step:
    Nsteps = np.size(steps[start_index:])-1
    Ny = np.size(y,1)
    time   = np.empty(Nsteps)
    b12 = np.empty(Nsteps)
    i = 0
    if rank == 0:
        print("Time \t b12_center\t b12_max")
    for step in steps[start_index:-1]:
        reader.step = step
        time[i] = reader.time
        
        rho, u, v, w = reader.readData( ('rho', 'u', 'v', 'w') )
        rho_bar = stats.reynolds_average(avg, rho)
        u_tilde = stats.favre_average(avg, rho, u, rho_bar=rho_bar)
        v_tilde = stats.favre_average(avg, rho, v, rho_bar=rho_bar)
        w_tilde = stats.favre_average(avg, rho, w, rho_bar=rho_bar)
        upp = u - u_tilde
        vpp = v - v_tilde
        wpp = w - w_tilde

        u2 = stats.favre_average(avg, rho, upp*upp + vpp*vpp + wpp*wpp, rho_bar=rho_bar)
        uv = stats.favre_average(avg, rho, upp*vpp, rho_bar=rho_bar)
        b12_vec = uv/u2;

        if rank == 0:
            b12[i] = b12_vec[0,Ny/2,0] 
            print("{} \t {} \t {}".format(reader.time,b12[i],max(b12_vec[0,:,0])))
            i += 1 
            # Save the vector
            #np.savetxt(filename_prefix+"%04d_b12.dat"%step,np.squeeze(b12_vec),delimiter=' ')
    
    # Write to file 
    if rank==0:
        array2save = np.empty((Nsteps,2))
        array2save[:,0] = time
        array2save[:,1] = b12
        outputfile = filename_prefix + 'anisotropy.dat'
        np.savetxt(outputfile,array2save,delimiter=' ')
        print("Done writing to {}".format(outputfile))

