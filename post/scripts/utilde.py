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
from SettingLib import NumSetting
from PoissonSol import * 

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
        print "Computes utilde" 
        print "  python {} <prefix> [tID_list (csv)] ".format(sys.argv[0])
        sys.exit()
    start_index = 0;
    if len(sys.argv) > 2:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
    filename_prefix = sys.argv[1]
    
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
    inp = nml.inputs(dirname,verbose=(rank==0))
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    du = inp.du
    if rank==0: print("\tdu = {}".format(inp.du))

    # Get the grid partition information
    nx,ny,nz = reader._full_chunk_size
    nblk_x = int(np.round(Nx/nx))
    nblk_y = int(np.round(Ny/ny))
    nblk_z = int(np.round(Nz/nz))
    if rank==0: print("Processor decomp: {}x{}x{}".format(nblk_x,nblk_y,nblk_z))

    # Compute stats at each step:
    i = 0
    thresh = 0.15
    for tid in tID_list:
        reader.step = tid
        r,u = reader.readData( ('rho','v') )
        utilde = stats.favre_average(avg,r,u)
        utilde = np.squeeze(utilde)
        
        offset = 50
        print(utilde[offset])
        print(utilde[-offset])
        offset = 20
        print(utilde[offset])
        print(utilde[-offset])
        sys.exit()
        # now gather from all processes
        root=0
        comm.Barrier()
        recvbuf=comm.gather(utilde, root)
        comm.Barrier()
        recvbuf = np.array(recvbuf)

        if rank==0:
            try: np.shape(recvbuf)[1] #should be a 2d array
            except:
                print("ERROR: Shape mismatch, recvbuf {}".format(np.shape(recvbuf)))
                sys.exit()
            
            # Stack the vectors into the correct order
            R22_array = np.reshape(recvbuf,[nblk_x,nblk_y,nblk_z,ny],order='F')
            
            # Now concat one column
            R22_mean = R22_array[0,0,0,:];
            if (nblk_y>1):
                for i in range(1,nblk_y):
                    mat = np.squeeze(R22_array[0,i,0,:])
                    R22_mean = np.hstack([R22_mean,mat])
            utilde = R22_mean

            if (np.shape(utilde)[0] < Ny):
                print("ERROR: Shape mismatch, utilde {}".format(np.shape(utilde)))
                sys.exit()
            else:
                outputfile = filename_prefix + "utilde_%04d.dat"%tid
                np.savetxt(outputfile,utilde,delimiter=' ')
                print("Done writing to {}".format(outputfile))
