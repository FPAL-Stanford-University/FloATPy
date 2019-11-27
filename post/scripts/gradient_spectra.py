#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys
sys.path.insert(0,'/home/kmatsuno/h5py/build/lib.linux-x86_64-2.7/h5py/')
import h5py

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red

from SettingLib import NumSetting
from PoissonSol import * 
import statistics as stats
import get_namelist as nml
from common import *

"""

        reader.step = tid
        rho, u, v, w = reader.readData( ('rho', 'u', 'v', 'w') )
        vorticity = der.curl(u, v, w, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        enstrophy = stats.get_enstrophy(vorticity, reader.grid_partition, dx, dy, dz)
        enstrophy = enstrophy[:,ic,:]
"""


debug = False
def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_start(default=0)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    start_index = 0
    if len(sys.argv) > 2:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
    
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

    # Set up the derivative object
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)

    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=(rank==0))
    du = inp.du
    if rank==0: print("du = {}".format(inp.du))
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    
    yplot = np.linspace(Ly/2.,-Ly/2., int(Ny)) 
    szx,szy,szz = np.shape(x)
    xmin = 0.;      xmax = Lx
    ymin = -Ly/2.;  ymax = Ly/2.
    zmin = 0.;      zmax = Lz

    chunkx = abs(x[0,0,0]-x[-1,0,0])
    chunky = abs(y[0,0,0]-y[0,-1,0])
    chunkz = abs(z[0,0,0]-z[0,0,-1])
    blkID = [int((np.amax(x)-xmin)/chunkx)-1,
            int((np.amax(y)-ymin)/chunky)-1,
            int((np.amax(z)-zmin)/chunkz)-1]
    if debug:
        comm.Barrier()
        print('rank {} is block ID {}'.format(rank,blkID))

    # Setup the fft object
    if rank==0: print('Setting up fft...')
    Nx,Ny,Nz = reader.domain_size
    settings = NumSetting(comm, reader.grid_partition,
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)
    ffto = PoissonSol(settings)
    nx,ny,nz = ffto.size_3d
    nblk_x = int(np.ceil(float(Nx)/nx))
    nblk_y = int(np.ceil(float(Ny)/ny))
    nblk_z = int(np.ceil(float(Nz)/nz))
    
    #if nblk_y == 1: sys.exit()
    # set up the writer
    #settings_writer = NumSetting(comm, reader.grid_partition,
    #         NX=Nx, NY=1, NZ=Nz,
    #         XMIN=0,        XMAX=Lx,
    #         YMIN=-Ly/2.,   YMAX=Ly/2.,
    #         ZMIN=0,        ZMAX=Lz,
    #         order=10)
    #writer = h5_writer(settings_writer)
    
    for tID in tID_list:
        reader.step = tID 
        
        # density gradient 
        q = reader.readData( ('rho') )
        drdx,drdy,drdz = der.gradient(q[0], x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)

        # Spectra in x,z:
        drdx_hat = ffto._fftX(drdx) 
        drdz_hat = ffto._fftZ(drdz) 
        drdx,drdy,drdz = None,None,None

        # Get centerline slice
        ic,yc = get_centerline(dirname,yplot,tID)
        yc = float(yc) # global centerline y
        this_ymin = y[0,0,0]-dy
        this_ymax = y[0,-1,0]
        rank_yc = ((yc>=this_ymin) and (yc<=this_ymax))
        print('rank {} with yrange ({},{}) contains yc={}? {}'.format(
            rank,y[0,0,0],y[0,-1,0],yc,rank_yc))
        if rank_yc:
            specx = abs(np.power(drdx_hat,2))
            specz = abs(np.power(drdz_hat,2))
        else:
            specx = np.zeros([szx,szy,szz])
            specz = np.zeros([szx,szy,szz])

        # Collect a snapshot along centerline
        blankx = np.zeros([Nx,Nz])
        blankz = np.zeros([Nx,Nz])
        if rank_yc:
            idy =  np.argmin(abs(y[0,:,0]-yc))
            idx0 = blkID[0]*szx
            idz0 = blkID[2]*szz
            idx = min(Nx,(blkID[0]+1)*szx)
            idz = min(Nz,(blkID[2]+1)*szz)
            blankx[idx0:idx,idz0:idz] = specx[:,idy,:]
            blankz[idx0:idx,idz0:idz] = specz[:,idy,:]
        planex = comm.reduce(blankx,root=0)
        planez = comm.reduce(blankz,root=0)
        if rank==0:
            data2save = np.zeros([2,Nx,Nz])
            data2save[0,:,:] = planex
            data2save[1,:,:] = planez
            outputfile = filename_prefix + 'grad_spec_%04d.dat'%tID
            np.save(outputfile, data2save)
            print('Saved to %s'%outputfile)
