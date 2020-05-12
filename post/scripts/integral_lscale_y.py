from mpi4py import MPI
import numpy as np
import numpy.linalg as la
import os
import sys
import matplotlib.pyplot as plt

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
  
def get_L99(y,utilde):
    utilde = np.squeeze(utilde)
    du = abs(utilde[-1]-utilde[0])
    utop = 0.99*du/2.
    ubot = -0.99*du/2.
    Ny = np.size(y)
    ibot = np.argmin(abs(utilde[Ny/2:]-utop)[::-1])
    itop = np.argmin(abs(utilde[:Ny/2]-ubot)[::-1])+Ny/2
    #itop = np.argmin(abs(utilde[Ny/2:]-utop))+Ny/2
    #ibot = np.argmin(abs(utilde[:Ny/2]-ubot))
    L99 = y[itop]-y[ibot]
    if L99<0: print('utilde or y misoriented. exiting'); sys.exit()
    return L99, itop, ibot

def window_tukey(n,N):
    alpha = 0.7
    tmp = alpha*N/2
    pi = np.pi
    
    if n<tmp:
        return 0.5*(1.+np.cos(pi*(n/tmp-1.)))
    elif n>=tmp and n<N*(1.-alpha/2.):
        return 1.
    elif n>=N*(1.-alpha/2.) and n<=N:
        return 0.5*(1.+np.cos(pi*(n/tmp-2./alpha+1.)))
    else:
        return 0

# Window if supersonic
def window_field(filename_prefix,ychunk,yplot,du,tID,q,utilde=None,pad=0):
    if utilde is None:
        fname = filename_prefix+'utilde_%04d.dat'%tID
        try:
            utilde = np.fromfile(fname,count=-1,sep=' ')
        except:
            if rank==0: 
                print('Write {}'.format(fname))
                sys.exit()
    L99, itop, ibot = get_L99(yplot,utilde)
    #print(L99, itop,ibot,y[itop],y[ibot])
    Ny = np.size(utilde)
    
    # Window based off 99% thickness
    ibot = max(ibot-pad,0)
    itop = min(itop+pad,Ny)
    N = abs(itop-ibot)      # total window size
    ny = np.shape(q)[1]     # size of chunk
    for iy_,ylocal in zip(range(ny),ychunk[0,:,0]):
        if ylocal>yplot[ibot] and ylocal<yplot[itop]:
            iy = np.argmin(abs(ylocal-yplot)) #global y index
            idx = iy-ibot # distance from bottom bound
            assert idx>=0
            assert idx<=N
            q[:,iy_,:] *= window_tukey(idx,N)
        else:
            q[:,iy_,:] = 0
    return q

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

    # Get the y coordinates
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    yplot = np.linspace(-Ly/2,Ly/2,int(Ny))
   
    for tid in tID_list:
        reader.step = tid
        if varname=='rho':
            q = reader.readData( ('rho'))
            rbar = stats.reynolds_average(avg,q[0])
            qpp = q[0] - rbar
        else:
            r,q = reader.readData( ('rho',varname) )
            qtilde = stats.favre_average(avg,r,q)
            qpp = q - qtilde

        # Window if supersonic
        if inp.Mc>=0.4:
            if inp.Mc==0.8: pad=0
            else: pad=Ny/10
            qpp=window_field(filename_prefix,y,yplot,inp.du,tid,
                    qpp,utilde=None,pad=pad)

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
            
            # Now concat one column
            R22_mean = R22_array[0,0,0,:];
            if (nblk_y>1):
                for i in range(1,nblk_y):
                    mat = np.squeeze(R22_array[0,i,0,:])
                    R22_mean = np.vstack([R22_mean,mat])
            
            R22_mean /= np.amax(R22_mean)
            if varname=='rho':
                outputfile = filename_prefix + "lscale_rr_%04d.dat"%tid
            else:
                outputfile = filename_prefix + "lscale_"+varname+varname+"_%04d.dat"%tid
            np.savetxt(outputfile,R22_mean,delimiter=' ')
            print("Done writing to {}".format(outputfile))
