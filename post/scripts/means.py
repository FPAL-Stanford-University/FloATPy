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
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
    
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

    # setup the inputs object, get grid info
    dirname = os.path.dirname(filename_prefix)
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    Ny = int(Ny)
    if rank==0: verbose=True
    else: verbose=False
    inp = nml.inputs(dirname,verbose)
    du = inp.du
    if rank==0: print("du = {}".format(inp.du))
    
    # Get the grid partition information
    nx,ny,nz = reader._full_chunk_size
    nblk_x = int(np.round(Nx/nx))
    nblk_y = int(np.round(Ny/ny))    
    nblk_z = int(np.round(Nz/nz))
    if rank==0: print("Processor decomp: {}x{}x{}".format(nblk_x,nblk_y,nblk_z))
    #if nblk_y==1: sys.exit()

    # Compute stats at each step:
    recvbuf = None
    for step in tID_list: 
        reader.step = step
        
        # density and streamwise vel, means
        r, u, v, w, p, T = reader.readData( ('rho','u', 'v','w','p','T') )
        u = (u-inp.du/2.)/inp.du
        means = np.zeros([np.shape(u)[1],6],dtype='f')
        means[:,0] = np.squeeze(stats.reynolds_average(avg,r))
        means[:,1] = np.squeeze(stats.reynolds_average(avg,u))
        means[:,2] = np.squeeze(stats.reynolds_average(avg,v))
        means[:,3] = np.squeeze(stats.reynolds_average(avg,w))
        means[:,4] = np.squeeze(stats.reynolds_average(avg,p))
        means[:,5] = np.squeeze(stats.reynolds_average(avg,T))
       
        # now gather from all processes for each Rij into another array
        recvbuf = {}
        root=0
        for i in range(6):
            recvbuf[i]=comm.gather(means[:,i], root=0)
            comm.Barrier()
            recvbuf[i] = np.array(recvbuf[i])

        if rank==0:
            total_array = np.zeros([Ny,6],dtype='f')
            for j in range(6):
                vec_array = np.reshape(recvbuf[j],[nblk_x,nblk_y,nblk_z,ny],order='F')
            
                # Now concat one column
                vec = vec_array[0,0,0,:];
                if (nblk_y>1):
                    for jj in range(1,nblk_y):
                        mat = np.squeeze(vec_array[0,jj,0,:])
                        vec = np.hstack([vec,mat])

                total_array[:,j] = vec
                print(vec)
                print("i={}\t{}".format(j,np.sum(vec)))
           
            print(total_array)
            outputfile = filename_prefix + 'means_%04d.dat'%step
            np.savetxt(outputfile,total_array,delimiter=' ')
            print("{}".format(outputfile))
