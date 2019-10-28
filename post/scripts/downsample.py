#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys
sys.path.insert(0,'/home/kmatsuno/h5py/build/lib.linux-x86_64-2.7/h5py/')
import h5py
from shutil import copyfile
from scipy.interpolate import interp1d

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml
from common import *
from h5test import h5_writer

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
    
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print "Usage: "
        print "  python {} <dir_read> [tID_read] <dir_write> ".format(sys.argv[0])
        sys.exit()
    dir_in  = sys.argv[1]
    tID     = int(sys.argv[2])
    dir_out = sys.argv[3]
    fname_in = dir_in + '/restart_%04d'%tID + '.h5' 
    fname_out = dir_out + '/restart_0000.h5' 
    
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
    serial_reader = por.PadeopsReader(dir_in+'/restart_',
            periodic_dimensions=periodic_dimensions)
    reader = serial_reader#pdr.ParallelDataReader(comm, serial_reader)
    #avg = red.Reduction(reader.grid_partition, periodic_dimensions)
    steps = sorted(reader.steps)

    # setup the inputs object
    inp = nml.inputs(dir_in,verbose=True)
    Nx0,Ny0,Nz0,Lx,Ly,Lz = nml.read_grid_params(dir_in,verbose=True)
    Nx,Ny,Nz,Lx,Ly,Lz, = nml.read_grid_params(dir_out,verbose=True)
    if rank==0: print('Original size: {}'.format([Nx0,Ny0,Nz0])) 
    if rank==0: print('New size: {}'.format([Nx,Ny,Nz])) 
    y = np.linspace(-Ly/2,Ly/2,Ny) 
    y0 = np.linspace(-Ly/2,Ly/2,Ny0) 
    
    # set up the settings object and writer
    #if rank==0: print('Setting up writer...')
    #settings = NumSetting( comm, reader.grid_partition, 
    #         NX=Nx, NY=Ny, NZ=Nz,
    #         XMIN=0,        XMAX=Lx,
    #         YMIN=-Ly/2.,   YMAX=Ly/2.,
    #         ZMIN=0,        ZMAX=Lz,
    #         order=10)
    #writer = h5_writer(settings)
    
    # Write attributes time and step
    seedfile = h5py.File(fname_out, 'w' )
    reader.step = tID
    time = np.array([reader.time])
    step = np.array([tID])
    seedfile.attrs.create('Time',time.astype('>d'))
    seedfile.attrs.create('step',step.astype('>i4'))


    qlist = ('rhou','rhov', 'rhow', 'TE', 'rhoY_0001','rhoY_0002')
    for i,qname in enumerate(qlist): 
        # Read input restart file
        if rank==0: print('Reading ' + qname)
        q0 = reader.readData( qname )
        q0 = q0[0]

        # Downsample in x
        if rank==0: print('Downsampling x,z')
        qhat0 = np.fft.fft(np.fft.fft(q0,axis=xdir),axis=zdir)
        qhat = np.zeros([Nx,Ny0,Nz],dtype=np.complex64)
        qhat[0:Nx/2,:,0:Nz/2] = qhat0[0:Nx/2,:,0:Nz/2]
        qhat[Nx/2:,:,Nz/2:]   = qhat0[Nx0-Nx/2:,:,Nz0-Nz/2:]
        qtall = np.fft.ifft(np.fft.ifft(qhat,axis=xdir),axis=zdir)
        qhat0 = None
        qhat = None

        # Interpolate in y
        q = np.zeros([Nx,Ny,Nz])
        if rank==0: print('Downsampling y')
        for i in range(Nx):
            for k in range(Nz):
                interpolater = interp1d(y0,np.real(qtall[i,:,k]),kind='nearest')
                q[i,:,k] = interpolater(y)
        qtall = None

        # write the new field out size 808094096
        if rank==0: print("Writing out")
        q = np.transpose(np.array(np.real(q)))
        shape = np.shape(q)
        sz = shape#(2,shape[1],2) 
        tmp=seedfile.create_dataset(qname, shape, data=q, 
                dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
        q = None
        tmp=None
    seedfile.close()
    print('Wrote to {}'.format(fname_out))
