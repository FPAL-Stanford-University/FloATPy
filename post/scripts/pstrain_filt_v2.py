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
from h5test import h5_writer
from mom_decomp import gather_y_plane
from decorr_lscale_y import transpose2y

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
   
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:" 
        print "Computes pressure strain fields from filtered pressure and velocities" 
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
    xmin = 0.;      xmax = Lx
    ymin = -Ly/2.;  ymax = Ly/2.
    zmin = 0.;      zmax = Lz
    yplot = np.linspace(-Ly/2,Ly/2,int(Ny))
    kx = np.array([2*np.pi/Lx*j for j in range(int(Nx/2))])
    kz = np.array([2*np.pi/Lz*j for j in range(int(Nz/2))])
    
    # Get the grid partition information
    nx,ny,nz = reader._full_chunk_size
    szx,szy,szz = np.shape(x)
    
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
    
    # set up the writer
    if rank==0: print('Setting up writer...')
    writer = h5_writer(settings)

    def filter_q_all(p,kc):
        tmpx = ffto._fftX(p)
        p = None
        tmpy = ffto._fftY(tmpx)
        tmpx = None
        phat = ffto._fftZ(tmpy)
        tmpy = None
        comm.Barrier()
        for ix,kx in enumerate(ffto.kx):
            if abs(kx)>kc: phat[ix,:,:] = 0. 
        for iy,ky in enumerate(ffto.ky):
            if abs(ky)>kc: phat[:,iy,:] = 0.
        for iz,kz in enumerate(ffto.kz):
            if abs(kz)>kc: phat[:,:,iz] = 0.
        pfilt = ffto._ifftX(ffto._ifftY(ffto._ifftZ(phat))) 
        return np.real(pfilt)
    def filter_q_ydir(p,kc): # for a 3D field p
        phat = ffto._fftY(p)
        print('rank {} done with fft'.format(rank))
        p = None
        for iy,ky in enumerate(ffto.ky):
            if abs(ky)>kc: phat[:,iy,:] = 0.
        print('rank {} done with filtering'.format(rank))
        pfilt = ffto._ifftY(phat) 
        print('rank {} done with ifft'.format(rank))
        return np.real(pfilt)
    def filter_q_xzdir(q,kc): #q must be a 2D xz plane
        qhat = np.fft.fft(np.fft.fft(q,axis=0),axis=1)
        for ix,kx in enumerate(kx):
            if abs(kx)>kc: qhat[ix,:] = 0. 
        for iz,kz in enumerate(kz):
            if abs(kz)>kc: qhat[:,iz] = 0.
        qfilt = np.fft.fft(np.fft.fft(q,axis=0),axis=1)
        return np.real(qfilt)

    def get_cutoff(q,thresh=0.8):# 80% of energy for an xz plane
        spec = np.squeeze(abs(np.fft.fft(q,axis=0)))
        spec = np.mean(spec[:Nx/2,:Nz/2],axis=1)
        spec *= kx

        # integrate the spectrum
        I = 0.
        cumulative = np.zeros(Nx/2)
        for i in range(Nx/2-2):
            I += spec[i]*abs(kx[i]-kx[i+1])
            cumulative[i] = I

        # Normalize the cumulative spectrum andfind 80% point
        cumulative /= np.amax(cumulative)
        ic = np.argmin(abs(cumulative-thresh))
        return kx[ic]
    
    # Compute stats at each step:
    for tID in tID_list:
        if rank==0: print('Reading data...')
        reader.step = tID
        #r,u,v,p = reader.readData( ('rho','u','v','p') )
        q = reader.readData( ('p') ); p=q[0]
        #rbar = stats.reynolds_average(avg,r)
        pbar = stats.reynolds_average(avg,p)
        #utilde = stats.favre_average(avg,r,u,rho_bar=rbar)
        #vtilde = stats.favre_average(avg,r,v,rho_bar=rbar)
        #upp = np.array(u-utilde)
        #vpp = np.array(v-vtilde)
        p = np.array(p-pbar)
        #r = None
        
        # Cutoff wavenumber
        kc = 11.46
        # Filter in y first, then gather, then xz planes
        if rank==0: print('Filtering y') 
        tmp = filter_q_ydir(p,kc)
        print('rank {} done'.format(rank))
        sys.exit() 

        # Unfiltered pressure and strain
        s,planes={},{}
        if rank==0: print('Compute unfiltered fluct strain...')
        qx,qy,qz = der.gradient(upp, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        s['11'] = qx
        s['12'] = qy/2.
        qx,qy,qz = der.gradient(vpp, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        s['22'] = qy
        s['12'] += qx/2.
        qx,qy,qz = None,None,None
        upp,vpp = None,None
        planes['s11'] = gather_y_plane(comm,reader,s['11'],y,yslice=0)
        planes['s22'] = gather_y_plane(comm,reader,s['22'],y,yslice=0)
        planes['s12'] = gather_y_plane(comm,reader,s['12'],y,yslice=0)
        planes['p']   = gather_y_plane(comm,reader,p,y,yslice=0)
        
        # Cutoff wavenumber
        if rank==0: kc = get_cutoff(planes['p'])
        else: kc = None
        kc = comm.bcast(kc,root=0)
        if rank==0: print('Filtering with cutoff kc={}'.format(kc))

        # Filter in y first, then gather, then xz planes
        if rank==0: print('Filtering y') 
        tmp = filter_q_ydir(s['11'],kc)
        if rank==0: print('Filtering xz') 
        qplane = gather_y_plane(comm,reader,tmp,y,yslice=0)
        if rank==0: qfilt = filter_q_ydir(qplane,kc)
        print('Done')
        sys.exit()
        p,s = None,None

        if rank==0: 
            print('Writing to file...')
            dir_out = dirname.split('/projects/ShockInducedMix/')[-1]
            dir_out = '/home/kmatsuno/' + dir_out + '/'
            outputfile = dir_out+ '/xz_pstrain_%04d.h5'%tID
            hf = h5py.File(outputfile, 'w')
            for key in planes.keys():
                hf.create_dataset(key, data=planes[key])
            hf.close()
            print('Saved to %s'%outputfile)
