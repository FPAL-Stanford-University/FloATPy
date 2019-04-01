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

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
    
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage:" 
        print "This script writes the full field of q, lambda2 criteria" 
        print "tID_list: csv, no spaces, eg [0,1,5]" 
        print "  python {} <prefix> [tID_list]".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    tid_list = map(int, sys.argv[2].strip('[]').split(',')) 
    
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
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)
    
    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=True)
    du = inp.du
    if rank==0: print("\tdu = {}".format(inp.du))

    # Preallocate for means and derivatives
    Nsteps = np.size(tid_list)
    Nx,Ny,Nz = reader.domain_size

    # Compute stats at each step:
    i = 0
    for tid in tid_list:
        reader.step = tid
       
        # vel gradients 
        u, v, w = reader.readData( ('u','v','w') )
        ux,uy,uz = der.gradient(u) #
        vx,vy,vz = der.gradient(v) #
        wx,wy,wz = der.gradient(w) #

        # Strain rate 
        S = np.zeros([Nx,Ny,Nz,9])
        S[:,:,:,0] = 2*ux;  S[:,:,:,1] = uy+vx; S[:,:,:,2] = uz+wx 
        S[:,:,:,3] = uy+vx; S[:,:,:,4] = 2*vy;  S[:,:,:,5] = vz+wy  
        S[:,:,:,6] = uz+wx; S[:,:,:,7] = vz+wy; S[:,:,:,8] = 2*wz
        S = 0.5*S
        S2 = la.norm(S,axis=3)# frob norm

        # Rotation tensor 
        R = np.zeros([Nx,Ny,Nz,9])
        zero = np.zeros([Nx,Ny,Nz])
        R[:,:,:,0] = zero;  R[:,:,:,1] = uy-vx; R[:,:,:,2] = uz-wx 
        R[:,:,:,3] =-uy+vx; R[:,:,:,4] = zero;  R[:,:,:,5] = vz-wy  
        R[:,:,:,6] =-uz+wx; R[:,:,:,7] =-vz+wy; R[:,:,:,8] = zero 
        R = 0.5*R
        R2 = la.norm(R,axis=3)

        # Lambda2 crit
        lambda2crit = np.zeros([Nx,Ny,Nz])
        for i in range(Nx-1):
            for j in range(Ny/4,3*Ny/4,1):
                for k in range(Nz-1):
                    r = np.reshape(R[i,j,k,:], [3,3])
                    s = np.reshape(S[i,j,k,:], [3,3])
                    M = 0.5*(la.matrix_power(r,2)+la.matrix_power(s,2))

                    eigvals,eigvecs = la.eig(M)
                    eigvals = sorted(eigvals)
                    lambda2crit[i,j,k] = eigvals[1]
        
        # Write to file 
        if rank==0:
            outputfile  = filename_prefix+"%04d_criterion_lambda2.npy"%tid
            np.save(outputfile,lambda2crit)
            
            print("Time: {}\n\tDone writing to {}".format(reader.time, outputfile))
            i = i+1; 
