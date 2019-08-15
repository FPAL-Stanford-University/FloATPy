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
    if len(sys.argv) < 4:
        print "Usage:" 
        print "Computes the integral correlation lengths of a primitive var" 
        print "  python {} <prefix> [tID_list (csv)] varname ".format(sys.argv[0])
        sys.exit()
    if len(sys.argv) > 3:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
        varname  = sys.argv[3] 
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

    # Preallocate for means and derivatives
    Nsteps = np.size(tID_list)
    Nx,Ny,Nz = reader.domain_size
    time = np.zeros([Nsteps])
    Lx = np.zeros([Nsteps])
    Ly = np.zeros([Nsteps])
    Lz = np.zeros([Nsteps])
    
    # Compute stats at each step:
    i = 0
    thresh = 0.15
    for tid in tID_list:
        reader.step = tid
        r,q = reader.readData( ('r',varname) )
        qtilde = stats.favre_average(avg,r,q)
        qpp = q - qtilde

        # y integral lengthscale
        vpp_hat = ffto._fftY(qpp)
        R22hat = np.abs(vpp_hat)
        R22 = np.abs(ffto._ifftY(R22hat))
        R22_mean = stats.reynolds_average(avg,R22)
        R22_mean = np.squeeze(R22_mean)

        # now gather from all processes
        root=0
        comm.Barrier()
        recvbuf=comm.gather(R22_mean,root)
        recvbuf = np.array(recvbuf)
        comm.Barrier()
        if rank==0:
            try: np.shape(recvbuf)[1] #should be a 2d array
            except:
                print("ERROR: Shape mismatch, recvbuf {}".format(np.shape(recvbuf)))
                sys.exit()
        
            # Stack the vectors into the correct order
            R22_array = np.reshape(recvbuf,[nblk_x,nblk_y,nblk_z,ny],order='F')
            #for zi in range(nblk_z):
            #    for yi in range(nblk_y):
            #        for xi in range(nblk_x):
            #            idx = xi + yi*(nblk_x) + zi*(nblk_x*nblk_y)
            #            R22_array[xi,yi,zi,:] = recvbuf[idx,:]
            
            # Now concat one column
            R22_mean = R22_array[0,0,0,:];
            if (nblk_y>1):
                for i in range(1,nblk_y):
                    mat = np.squeeze(R22_array[0,i,0,:])
                    R22_mean = np.vstack([R22_mean,mat])
            
            print("Shape of R22_mean: {}".format(np.size(R22_mean)))
            R22_mean /= np.amax(R22_mean)

            outputfile = filename_prefix + "lscale_"+varname+varname+"_%04d.dat"%tid
            np.savetxt(outputfile,R22_mean,delimiter=' ')
            print("Done writing to {}".format(outputfile))
