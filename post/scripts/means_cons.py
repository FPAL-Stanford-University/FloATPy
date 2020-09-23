#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys
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

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_list(optional)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    start_index = 0
    tID_list = None
    if len(sys.argv) > 2:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
    
    dirname = os.path.dirname(filename_prefix)
    if 'Mc04' in dirname:
        dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/temporal/')[-1]
        dir_out = '/home/kmatsuno/ShearLayerData/temporal/' + dir_out + '/'
    else:
        dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/mira/')[-1]
        dir_out = '/home/kmatsuno/ShearLayerData/production/' + dir_out + '/'
    
    periodic_dimensions = (True,False,True)
    x_bc = (0,0)
    y_bc = (1,1)
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
    if tID_list is None: tID_list = sorted(reader.steps)

    # Set up the derivative object
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)

    # setup inputs and writer
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=(rank==0))
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    du = inp.du
    if rank==0: print("\tdu = {}".format(inp.du))
    nx,ny,nz = reader._full_chunk_size
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)
    
    # Store means:
    means = np.zeros([6,Ny])
    qlist = {'rhoY_0001':0, 'rhoY_0002':1, 'rhou':2,'rhov':3,'rhow':4,'TE':5}
    for tID in tID_list: 
        reader.step = tID
        
        # Density first
        #print("Reading density")
        r1, r2  = reader.readData( ('rhoY_0001','rhoY_0002') )
        if (procs>1) and (np.shape(r1)[1]<Ny): 
            r1 = transpose2y(settings,r1)
            r2 = transpose2y(settings,r2)
        r = r1+r2; 
        rbar = stats.reynolds_average(avg,r)

        for key in qlist.keys():
            tmp = reader.readData(key)
            q = tmp[0]
            if (procs>1) and (np.shape(q)[1]<Ny): 
                q = transpose2y(settings,q)
            qtilde = stats.favre_average(avg,r,q/r,rho_bar=rbar)
            idx = qlist[key]
            means[idx,:] = np.squeeze(qtilde)
        
        if rank==0: 
            outputfile = dir_out+"restart_means_%04d.h5"%tID
            hf = h5py.File(outputfile, 'w')
            for key in qlist.keys():
                hf.create_dataset(key, data=means[qlist[key],:])
            time = reader.time
            hf.attrs.create('Time',time.astype('>d'))
            hf.close()
            print('Saved to %s'%outputfile)
    print('Done')
