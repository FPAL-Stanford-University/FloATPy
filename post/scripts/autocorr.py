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

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
   
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:" 
        print "Computes autocorr at centerline" 
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

    # Get the grid partition information
    nx,ny,nz = reader._full_chunk_size
    szx,szy,szz = np.shape(x)
    nblk_x = int(np.round(Nx/(szx-1)))
    nblk_y = int(np.round(Ny/(szy-1)))    
    nblk_z = int(np.round(Nz/(szz-1)))
    
    chunkx = abs(x[0,0,0]-x[-1,0,0])
    chunky = abs(y[0,0,0]-y[0,-1,0])
    chunkz = abs(z[0,0,0]-z[0,0,-1])
    blkID = [int((np.amax(x)-xmin)/chunkx)-1,
            int((np.amax(y)-ymin)/chunky)-1,
            int((np.amax(z)-zmin)/chunkz)-1]

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
    kxmax = ffto.NX * np.pi/ffto.LX
    kymax = ffto.NY * np.pi/ffto.LY
    kzmax = ffto.NZ * np.pi/ffto.LZ
    if rank==0: 
        print("Domain size: {}x{}x{} ({}x{}x{})".format(
        ffto.LX,ffto.LY,ffto.LZ,
        ffto.NX,ffto.NY,ffto.NZ))
        print('Processor decomposition: {}x{}x{}'.format(
        nblk_x,nblk_y,nblk_z))
    rank_0x = (0 in ffto.kx)
    rank_0y = (0 in ffto.ky)
    rank_0z = (0 in ffto.kz)
    thresh = 1e-10
    rank_oddx = (abs(kxmax -np.amax(abs(ffto.kx)))<thresh)
    rank_oddy = (abs(kymax -np.amax(abs(ffto.ky)))<thresh)
    rank_oddz = (abs(kzmax -np.amax(abs(ffto.kz)))<thresh)

    def gather_centerline(comm,dat2save,ic=None):
        blank = np.zeros([Nx,Nz])
        if rank_oddy:
            if ic is None: ic = np.argmin(abs(y[0,:,0]))
            idx0 = blkID[0]*szx
            idz0 = blkID[2]*szz
            idx = min(Nx,(blkID[0]+1)*szx)
            idz = min(Nz,(blkID[2]+1)*szz)
            blank[idx0:idx,idz0:idz] = dat2save[:,ic,:]
        plane = comm.reduce(blank,root=0)
        return plane
    
    # Compute stats at each step:
    for tID in tID_list:

        if rank==0: print('Reading data...')
        reader.step = tID
        r,v = reader.readData( ('rho','v') )
        vtilde = stats.favre_average(avg,r,v)
        vpp = np.array(v-vtilde)
        r,v = None,None

        qDict = {}
        qDict['qx'] = ffto._ifftX(abs(ffto._fftX(vpp)))
        qDict['qz'] = ffto._ifftZ(abs(ffto._fftZ(vpp)))

        tmp,yc =  get_centerline(dirname,yplot,tID)
        ic = np.argmin(abs(y[0,:,0]-yc))

        # Collect a snapshots along centerline and save
        planes = {}
        for key in qDict.keys():
            planes[key] = gather_centerline(comm,qDict[key],ic=ic)

        if rank==0: 
            print('Writing to file...')
            outputfile = dirname + '/autocorr_%04d.h5'%tID
            hf = h5py.File(outputfile, 'w')
            for key in planes.keys():
                hf.create_dataset(key, data=planes[key])
            hf.close()
            print('Saved to %s'%outputfile)
