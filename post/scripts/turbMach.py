from mpi4py import MPI
import numpy as np
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats

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
    outputfile  = filename_prefix + "Mt_growth.dat"

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

    x, y, z = reader.readCoordinates()
    Nx,Ny,Nz = reader.domain_size
    dx,dy,dz = grid_res(x,y,z)
    
    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=True)
    du = inp.du
    if rank==0: print("du = {}".format(inp.du))


    # Compute stats at each step:
    Nsteps = np.size(steps[start_index:])-1
    time = np.empty(Nsteps)
    Mt_max = np.empty(Nsteps)   # Max Mt across all y planes
    Mt_center = np.empty(Nsteps)# Mt at center plane Ny/2
    if rank == 0: print("Time \t Max Mt \t Center Mt") 
    
    i = 0
    for step in steps[start_index:-1]:
        reader.step = step
        time[i] = reader.time

        # TKE
        rho, u, v, w, T = reader.readData( ('rho', 'u', 'v', 'w','T') )
        tke = stats.TKE(rho, u, v, w, reader.grid_partition, 
                avg, dx, dy, dz, volume_integrated=False)
        
        # Speed of sound, xz average
        c3D = np.sqrt(inp.gam*T)
        c = stats.reynolds_average(avg, c3D)

        # Turbulent Mach number
        Mt = np.squeeze(tke/c)
        Mt_max[i] = np.amax(Mt);
        Mt_center[i] = Mt[Ny/2]

        if rank == 0:
            print("{} \t {} \t {}".format(reader.time,Mt_max[i],Mt_center[i]))
        i = i+1

    # Write to file 
    if rank==0:
        array2save = np.empty((Nsteps,3))
        array2save[:,0] = time
        array2save[:,1] = Mt_center
        array2save[:,2] = Mt_max
        np.savetxt(outputfile,array2save,delimiter=' ')
        print("Done writing to {}".format(outputfile))
