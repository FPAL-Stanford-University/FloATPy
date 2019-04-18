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
    outputfile  = filename_prefix + "growth.dat"
    
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
    du = inp.du
    if rank==0: print("du = {}".format(inp.du))

    # Compute stats at each step:
    Nsteps = np.size(steps[start_index:])-1
    time   = np.empty(Nsteps)
    dtheta = np.empty(Nsteps)
    domega = np.empty(Nsteps)
    i = 0
    if rank == 0:
        print("Time \t dtheta \t domega")
    for step in steps[start_index:-1]:
        reader.step = step
        time[i] = reader.time
        
        # density and streamwise vel, means
        rho, u = reader.readData( ('rho', 'u') )
        rho_bar = stats.reynolds_average(avg,rho)
        u_bar = stats.reynolds_average(avg,rho)
        utilde = stats.favre_average(avg,rho,u,rho_bar)
        
        # Compute momentum thickness
        I = rho_bar*(0.5*du-utilde)*(0.5*du+utilde)
        dtheta[i] = stats.integrate_y(I, dy, reader.grid_partition)/(inp.r_ref*inp.du**2)

        # vorticity thickness
        dudx,dudy,dudz = der.gradient(u, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        dudy_bar = stats.reynolds_average(avg,dudy)
        domega[i] = du/np.amax(np.abs(dudy_bar))

        if rank == 0:
            print("{} \t {} \t {}".format(reader.time,dtheta[i],domega[i]))
        i = i+1; 
    
    # Write to file 
    if rank==0:
        if 0:#Nsteps <= 5:
            print("Tested less than 5 steps, not writing to file.")
        else:
            array2save = np.empty((Nsteps,3))
            array2save[:,0] = time
            array2save[:,1] = dtheta
            array2save[:,2] = domega
            np.savetxt(outputfile,array2save,delimiter=' ')
            print("Done writing to {}".format(outputfile))
