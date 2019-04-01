from mpi4py import MPI
import numpy as np
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
#import floatpy.derivatives.explicit.explicit_derivative as ed
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats

debug = False

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz

class inputs:
    def __init__(self,Mc,Re,rr):
        gam = 1.4
        p_ref = 1.
        T_ref = 1.
        r_ref = 1.
        mu_ref = 1./Re
        Rgas1 = (1.+rr)/2.			
        Rgas2 = (1.+rr)/(2.*rr) 
        c1 = np.sqrt(gam*p_ref/(r_ref/Rgas1))
        c2 = np.sqrt(gam*p_ref/(r_ref/Rgas2))
        
        self.du = Mc*(c1+c2)
        self.gam = gam  
        self.p_ref = p_ref 
        self.T_ref = T_ref 
        self.r_ref = r_ref 
        self.mu_ref = mu_ref
        
    
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage:"
        print "This script writes means normalized by self-similar coords"
        print "tID_list: csv, no spaces, eg (0,1,5)" 
        print "  python {} <prefix> [tID_list]".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    tid_list = map(int, sys.argv[2].strip('[]').split(',')) 
    outputfile  = filename_prefix + "means.dat"
    
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

    # Set up compact derivative object w/ 10th order schemes
    # Set up explicit derivative object w/ 6th order schemes
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)
    
    # setup the inputs object
    Mc = 0.7; Re = 640; rr = 8 # should have been Re160
    inp = inputs(Mc,Re,rr)
    du = inp.du
    if rank==0:
        if debug: print("du = {}".format(inp.du))

    # Preallocate for means and derivatives
    Nsteps = np.size(tid_list)
    Nx,Ny,Nz = reader.domain_size
    time   = np.empty(Nsteps)
    eta    = np.empty(Nsteps)
    dtheta = np.empty(Nsteps)

    # Compute stats at each step:
    i = 0
    if rank == 0:
        print("Time")
    for tid in tid_list:
        reader.step = steps[tid]
        time[i] = reader.time
       
        # density and streamwise vel, means
        rho, u = reader.readData( ('rho', 'u') )
        rho_bar = stats.reynolds_average(avg,rho)
        u_bar = stats.reynolds_average(avg,rho)
        utilde= stats.favre_average(avg,rho,u,rho_bar)

        if rank == 0:
            print("{}".format(reader.time))
        i = i+1; 

    # Write to file 
    if rank==0:
        if Nsteps <= 5:
            print("Tested less than 5 steps, not writing to file.")
        else:
            array2save = np.empty((Ny,2*Nq+1))
            array2save[:,0:Nq] = np.squeeze(q_bar) 
            array2save[:,Nq+1:] = np.squeeze(q_tilde) 
            np.savetxt(outputfile,array2save,delimiter=' ')
            print("Done writing to {}".format(outputfile))
