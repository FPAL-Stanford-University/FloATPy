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

    # Setup the fft object
    if rank==0: print('Setting up fft...')
    Nx,Ny,Nz = reader.domain_size
    settings = NumSetting(comm, reader.grid_partition,
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=(Nx)*dx,
             YMIN=-Ny*dy/2.,YMAX=(Ny)*dy/2.,
             ZMIN=0,        ZMAX=(Nz)*dz,
             order=10)
    ffto = PoissonSol(settings)
    nx,ny,nz = ffto.size_3d
    nblk_x = int(np.ceil(float(Nx)/nx))
    nblk_y = int(np.ceil(float(Ny)/ny))
    nblk_z = int(np.ceil(float(Nz)/nz))
    if rank==0: 
        print("Domain size: {}x{}x{} ({}x{}x{})".format(
        ffto.LX,ffto.LY,ffto.LZ,
        ffto.NX,ffto.NY,ffto.NZ))
        print('Processor decomposition: {}x{}x{}'.format(
        nblk_x,nblk_y,nblk_z))

    # Centerline:
    ic = int(Ny)/2

    # Compute stats at each step:
    for tid in tID_list:
        reader.step = tid
        rho, u, v, w = reader.readData( ('rho', 'u', 'v', 'w') )
        vorticity = der.curl(u, v, w, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        enstrophy = stats.get_enstrophy(vorticity, reader.grid_partition, dx, dy, dz)
        enstrophy = enstrophy[:,ic,:]


       
        # print grid partition
        print(rank,reader.grid_partition)
        sys.exit()

        # now gather from all processes
        root=0
        comm.Barrier()
        recvbuf=comm.gather(enstrophy, root)
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
                outputfile = filename_prefix + "enstrophy_%04d.dat"%tid
                np.savetxt(outputfile,utilde,delimiter=' ')
                print("Done writing to {}".format(outputfile))
