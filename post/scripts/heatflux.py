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
from decorr_lscale_y import transpose2y
from common import *
   
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:" 
        print "Computes ubar,vbar,utilde,vtilde, r'u' and r'v'" 
        print "  python {} <prefix> [tID_list (csv)] ".format(sys.argv[0])
        sys.exit()
    start_index = 0;
    if len(sys.argv) > 2:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
    else: tID_list = None
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
    if tID_list is None: tID_list = steps

    # Set up compact derivative object w/ 10th order schemes
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)
    
    # setup the inputs object, get grid info
    dirname = os.path.dirname(filename_prefix)
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    Ny = int(Ny)
    inp = nml.inputs(dirname,verbose=(rank==0))
    du = inp.du
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)
    
    for tID in tID_list:
        reader.step = tID
        r,T,v = reader.readData( ('rho','T','v') )
        if (procs>1) and (np.shape(r)[1]<Ny): 
            if rank==0: print('Transposing')
            r = transpose2y(settings,r)
            T = transpose2y(settings,T)
            if rank==0: print('Done')
        rbar = stats.reynolds_average(avg,r)
        Tbar = stats.reynolds_average(avg,T)
        vbar = stats.reynolds_average(avg,v)
        Ttilde = stats.favre_average(avg,r,T,rho_bar=rbar)
        vtilde = stats.favre_average(avg,r,v,rho_bar=rbar)
        Tpp = T - Ttilde
        vpp = v - vtilde
        
        dat = np.zeros([np.shape(r)[1],6],dtype='f')
        dat[:,0] = np.squeeze(Tbar)
        dat[:,1] = np.squeeze(vbar)
        dat[:,2] = np.squeeze(Ttilde) 
        dat[:,3] = np.squeeze(vtilde)
        dat[:,4] = np.squeeze(stats.favre_average(avg,r,vpp*Tpp))
        #dat[:,5] = np.squeeze(stats.reynolds_average(avg,))
        
        if rank==0: 
            dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/mira/')[-1]
            dir_out = '/home/kmatsuno/ShearLayerData/production/' + dir_out + '/'
            outputfile = dir_out+"heatflux_%04d.dat"%tID
            print("Writing to {}".format(outputfile))
            np.savetxt(outputfile,np.squeeze(dat),delimiter=' ')
            print('Done')
