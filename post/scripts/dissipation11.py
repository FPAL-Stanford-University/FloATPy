#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml
from SettingLib import NumSetting
from decorr_lscale_y import transpose2y
from common import *
    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_start(default=0)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    start_index = 0
    if len(sys.argv) > 2:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
    else: tID_list = None
    
    dirname = os.path.dirname(filename_prefix)
    if 'Mc04' in dirname:
        dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/temporal/')[-1]
        dir_out = '/home/kmatsuno/ShearLayerData/temporal/' + dir_out + '/'
    else:
        dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/mira/')[-1]
        dir_out = '/home/kmatsuno/ShearLayerData/production/' + dir_out + '/'

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

    # Set up the derivative object
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

    # Compute stats at each step:
    for tID in tID_list: 
        reader.step = tID
        
        u, v, w = reader.readData( ('u','v','w') )

        bulk = 0.0
        mu =1./ inp.Re
        bmbda = (4./3.)*mu + bulk
        lmbda = bulk - (2./3.)*mu
        if rank==0: print('Computing grad u') 
        dudx,dudy,dudz = der.gradient(u, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        if rank==0: print('Computing grad v') 
        dvdx,dvdy,dvdz = der.gradient(v, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        if rank==0: print('Computing grad w') 
        dwdx,dwdy,dwdz = der.gradient(w, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        if rank==0: print('Computing tau_ij') 
        tau12 = mu*(dudy + dvdx); dudy,dvdx=None,None    
        tau13 = mu*(dudz + dwdx); dudz,dwdx=None,None 
        dvdz,dwdy=None,None
        tau11 = bmbda*dudx + lmbda*(dvdy + dwdz) 
        dudx,dvdy,dwdz=None,None,None


        if rank==0: print('Computing upp') 
        q = reader.readData( ('rho') )
        r = q[0]; q=None
        utilde = stats.favre_average(avg,r,u)
        upp = u - utilde
        r,u,v,w = None,None,None,None
        if rank==0: print('Computing grad upp')
        dudx,dudy,dudz = der.gradient(upp, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        tmp =2*( tau11*dudx + tau12*dudy + tau13*dudz)
       
        tau12,tau13,tau11=None,None,None

        if (procs>1) and (np.shape(tmp)[1]<Ny): 
            if rank==0: print('Transposing')
            tmp = transpose2y(settings,tmp)
        D11 = np.squeeze(stats.reynolds_average(avg,tmp))
        tmp=None

        if rank==0: 
            outputfile = dir_out+"/dissipation11_%04d.dat"%tID
            np.savetxt(outputfile,D11,delimiter=' ')
            print("Done writing to {}".format(outputfile))
