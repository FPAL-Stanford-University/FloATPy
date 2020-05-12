from mpi4py import MPI
import numpy as np
import numpy.linalg as la
import os
import sys
sys.path.insert(0,'/home/kmatsuno/h5py/build/lib.linux-x86_64-2.7/')
import h5py

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml
from SettingLib import NumSetting
from PoissonSol import * 
from common import *
from decorr_lscale_y import get_qp,get_qpp
from mom_decomp import gather_y_plane

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
   
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:" 
        print "Computes autocorr at centerline" 
        print "python {} <prefix> [tID_list (csv)] varname".format(sys.argv[0])
        sys.exit()
    if len(sys.argv) > 2:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
    filename_prefix = sys.argv[1]
    varname  = sys.argv[3] 
    
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
    xmin = 0.;      xmax = Lx
    ymin = -Ly/2.;  ymax = Ly/2.
    zmin = 0.;      zmax = Lz
    yplot = np.linspace(-Ly/2,Ly/2,int(Ny))

    # Get the grid partition information
    nx,ny,nz = reader._full_chunk_size
    szx,szy,szz = np.shape(x)
    nblk_x = int(np.round(Nx/(szx-1)))
    nblk_y = int(np.round(Ny/(szy-1)))    
    nblk_z = int(np.round(Nz/(szz-1)))
    
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

    # Compute stats at each step:
    for tID in tID_list:

        if rank==0: print('Reading data...')
        reader.step = tID
        if varname=='r' or varname=='p': qpp = get_qp(reader,avg,varname)
        else:  qpp = get_qpp(reader,avg,varname)
        
        # Get utilde
        fname = filename_prefix+'utilde_%04d.dat'%tID 
        utilde = np.fromfile(fname,count=-1,sep=' ')

        # get thicknesses
        L99,itop,ibot = get_L99(yplot,utilde)
        dtheta = get_dtheta(dirname,reader.time)

        # Collect a snapshots along centerline and save
        if inp.rr==1:
            yc = 0
        else:
            ic = np.argmin(abs(utilde))
            yc = yplot[ic]
        offset = L99/4.#dtheta/2.
        y0_list = [yc,yc+offset,yc-offset]

        corrx = np.zeros([Nx,3])
        corrz = np.zeros([Nz,3])
        xdir = 0
        zdir = 1
        for j,y0 in enumerate(y0_list): 
            if rank==0: print('Gather planes...')
            plane = gather_y_plane(comm,reader,np.squeeze(qpp),y,yslice=y0)
            
            if rank==0: 
                print('Calc fftx...')
                qhat = np.fft.fft(plane,axis=xdir)
                autocorr2D = np.real(np.fft.ifft(qhat*np.conj(qhat),axis=xdir))
                autocorr2D = np.abs(np.fft.fftshift(autocorr2D))
                for k in range(Nz): autocorr2D /= np.amax(autocorr2D[:,k])
                corrx[:,j] = np.squeeze(np.mean(autocorr2D,axis=zdir))
                
                print('Calc fftz...')
                qhat = np.fft.fft(plane,axis=zdir)
                autocorr2D = np.real(np.fft.ifft(qhat*np.conj(qhat),axis=zdir))
                autocorr2D = np.abs(np.fft.fftshift(autocorr2D))
                for k in range(Nz): autocorr2D /= np.amax(autocorr2D[:,k])
                corrz[:,j] = np.squeeze(np.mean(autocorr2D,axis=xdir))
                #qhat = np.fft.fft(plane,axis=1)
                #planes['qz'] = np.real(np.fft.ifft(qhat*np.conj(qhat),axis=1))

        if rank==0: 
            print('Writing out...')
            dir_out = dirname.split('/projects/ShockInducedMix/')[-1]
            dir_out = '/home/kmatsuno/' + dir_out + '/'
            if varname=='rho':
                outputfile = dir_out+"autocorr_rr_%04d.h5"%tID
            else:
                outputfile = dir_out+"autocorr_"+varname+varname+"_%04d.h5"%tID
            hf = h5py.File(outputfile, 'w')
            hf.create_dataset('corrx', data=corrx)
            hf.create_dataset('corrz', data=corrz)
            hf.close()
            print('Saved to %s'%outputfile)
