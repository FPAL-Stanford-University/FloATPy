#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys

sys.path.insert(0,'/home/kmatsuno/h5py/build/lib.linux-x86_64-2.7/')
import h5py

from shutil import copyfile

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml

from h5test import h5_writer
from SettingLib import NumSetting
from PoissonSol import * 

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
    
if __name__ == '__main__':
    if len(sys.argv) < 3:
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

    # Set up the seed reader and new reader
    seed_reader = por.PadeopsReader(dir_in+'/restart_',
        periodic_dimensions=periodic_dimensions)
    serial_reader = por.PadeopsReader(dir_out+'/restart_',
        periodic_dimensions=periodic_dimensions)
    seeder = seed_reader#pdr.ParallelDataReader(comm, seed_reader)
    #avg = red.Reduction(seeder.grid_partition, periodic_dimensions)
    reader = serial_reader#pdr.ParallelDataReader(comm, serial_reader)

    # Get parameters
    inp_s = nml.inputs(dir_in, verbose=(rank==0))
    inp   = nml.inputs(dir_out,verbose=(rank==0))
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dir_in,verbose=(rank==0))
   
    """
    if rank==0: print('Setting up writer...')
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)
    szx, szy, szz = settings.chunk_3d_size
    writer = h5_writer(settings)
    """
    
    # Write attributes time and step
    fname_out = dir_out+'seeded_%04d.h5'%tID
    seedfile = h5py.File(fname_out, 'w' )
    
    seeder.step = tID
    t = np.array([seeder.time])
    step = np.array([tID])
    seedfile.attrs.create('Time',t.astype('>d'))
    seedfile.attrs.create('step',step.astype('>i4'))

    # Read input restart file
    qlist = ('rhou', 'rhov', 'rhow', 'TE', 'rhoY_0001','rhoY_0002')
    for qname in qlist:
        # Get perturbations perturbations
        if rank==0: print('Reading seed %s'%qname)
        q0 = seeder.readData( qname )
        qbar = np.mean(np.mean(q0[0],axis=2),axis=0)#stats.reynolds_average(avg,q0[0])
        q = q0[0]-qbar[np.newaxis,:,np.newaxis]
        q0 = None

        # Get data to add perturbations to
        if rank==0: print('Reading new %s'%qname)
        qnew = reader.readData( qname )
        qbar = np.mean(np.mean(qnew[0],axis=2),axis=0)#stats.reynolds_average(avg,q0[0])
        q += qbar[np.newaxis,:,np.newaxis]
        qnew = None

        # transpose and write out
        print("Writing out")
        q = np.transpose(q)
        shape = np.shape(q)
        sz = (2,shape[1],2) 
        rudat = seedfile.create_dataset(qname, shape,data=q, 
                dtype=np.dtype('>d'),chunks=sz, fillvalue=0)

   
    """
    if rank==0: print('Writing to h5')
    q,q0 = None,None
    writer.update(tID,seeder.time)
    writer.saveFields(qDict,filePrefix=fname_out) 
    """
    seedfile.close()
