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
    
    # set up the writer
    if rank==0: print('Setting up writer...')
    writer = h5_writer(settings)

    def filter_q(p,kc):
        phat = ffto._fftX(ffto._fftY(ffto._fftZ(p)))
        for ix,kx in enumerate(ffto.kx):
            if abs(kx)>kc: phat[ix,:,:] = 0. 
        for iy,ky in enumerate(ffto.ky):
            if abs(ky)>kc: phat[:,iy,:] = 0.
        for iz,kz in enumerate(ffto.kz):
            if abs(kz)>kc: phat[:,:,iz] = 0.
        pfilt = ffto._ifftX(ffto._ifftY(ffto._ifftZ(phat))) 
        return np.real(pfilt)

    def gather_centerline(comm,dat2save):
        blank = np.zeros([Nx,Nz])
        if rank_oddy:
            idy = np.argmin(abs(y[0,:,0]))
            idx0 = blkID[0]*szx
            idz0 = blkID[2]*szz
            idx = min(Nx,(blkID[0]+1)*szx)
            idz = min(Nz,(blkID[2]+1)*szz)
            blank[idx0:idx,idz0:idz] = dat2save[:,idy,:]
        plane = comm.reduce(blank,root=0)
        return plane
    
    # Compute stats at each step:
    i = 0
    thresh = 0.15
    for tID in tID_list:

        if rank==0: print('Reading data...')
        reader.step = tID
        r,u,v,p = reader.readData( ('rho','u','v','p') )
        rbar = stats.reynolds_average(avg,r)
        pbar = stats.reynolds_average(avg,p)
        utilde = stats.favre_average(avg,r,u,rho_bar=rbar)
        vtilde = stats.favre_average(avg,r,v,rho_bar=rbar)
        upp = np.array(u-utilde)
        vpp = np.array(v-vtilde)
        p = np.array(p-pbar)
        r = None

        # Unfiltered pressure and strain
        qDict,s={},{}
        if rank==0: print('Compute unfiltered fluct strain...')
        qx,qy,qz = der.gradient(upp, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        s['11'] = qx
        s['12'] = qy/2.
        qx,qy,qz = der.gradient(vpp, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        s['22'] = qy
        s['12'] += qx/2.
        qx,qy,qz = None,None,None
        qDict['s11'] = s['11']
        qDict['s22'] = s['22']
        qDict['s12'] = s['12']
        qDict['p']   = p
        upp,vpp = None,None

        # filter p and pstrain
        lscale = Lx/5
        #kc = 2.*np.pi/lscale
        kc = ffto.kx[50]
        if rank==0: print('Filtering with cutoff kc={}'.format(kc))
        qDict['s11f'] = filter_q(s['11'],kc)
        qDict['s22f'] = filter_q(s['22'],kc)  
        qDict['s12f'] = filter_q(s['12'],kc)
        qDict['pf']   = filter_q(p,kc)

        # Collect a snapshots along centerline and save
        planes = {}
        for key in qDict.keys():
            planes[key] = gather_centerline(comm,qDict[key])

        if rank==0: 
            print('Writing to file...')
            outputfile = dirname + '/xz_pstrain_%04d.h5'%tID
            hf = h5py.File(outputfile, 'w')
            for key in planes.keys():
                hf.create_dataset(key, data=planes[key])
            hf.close()
            print('Saved to %s'%outputfile)
        
        # write 3d file
        #writer.update(tid,reader.time)
        #writer.saveFields(qDict,filePrefix=directory+'/pstrainFilt_') 
        #qDict = None
