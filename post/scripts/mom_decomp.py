#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import pickle

sys.path.insert(0,'/home/kmatsuno/h5py/build/lib.linux-x86_64-2.7/')
import h5py

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
reload(pdr)

import statistics as stats
import get_namelist as nml
from SettingLib import NumSetting
from PoissonSol import * 
from h5test import h5_writer
from common import *

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz

def mpi_max(comm,value):
    local_max = np.amax(value)
    global_max = comm.reduce(local_max, op=MPI.MAX, root=0)
    return global_max

# 99% thickness (might be diff from one in common)
def get_L99(y,utilde,du):
    utilde = np.squeeze(utilde)
    utop = 0.99*du/2.
    ubot = -0.99*du/2.
    Ny = np.size(y)
    itop = np.argmin(abs(utilde[Ny/2:]-utop))+Ny/2
    ibot = np.argmin(abs(utilde[:Ny/2]-ubot))
    L99 = abs(y[itop]-y[ibot])
    return L99, itop, ibot

def window_tukey(n,N):
    alpha = 0.7
    tmp = alpha*N/2
    pi = np.pi
    
    if n<tmp:
        return 0.5*(1.+np.cos(pi*(n/tmp-1.)))
    elif n>=tmp and n<N*(1.-alpha/2.):
        return 1.
    elif n>=N*(1.-alpha/2.) and n<=N:
        return 0.5*(1.+np.cos(pi*(n/tmp-2./alpha+1.)))
    else:
        return 0

# Window if supersonic
def window_field(filename_prefix,ychunk,yplot,du,tID,q,utilde=None):
    if utilde is None: 
        fname = filename_prefix+'utilde_%04d.dat'%tID 
        try:
            utilde = np.fromfile(fname,count=-1,sep=' ')
        except:
            if rank==0:
                print('Write {}'.format(fname))
            sys.exit()
    L99, itop, ibot = get_L99(yplot,utilde,du)
    Ny = np.size(utilde)

    # Window based off 99% thickness
    pad = 75
    ibot = max(ibot-pad,0)
    itop = min(itop+pad,Ny)
    N = abs(itop-ibot)      # total window size
    ny = np.shape(q)[1]     # size of chunk
    for iy_,ylocal in zip(range(ny),ychunk[0,:,0]):
        if ylocal>yplot[ibot] and ylocal<yplot[itop]:
            iy = np.argmin(abs(ylocal-yplot)) #global y index
            idx = iy-ibot # distance from bottom bound
            assert idx>=0
            assert idx<=N
            q[:,iy_,:] *= window_tukey(idx,N)
        else:
            q[:,iy_,:] = 0
    return q

def get_fluct(reader,avg):
    r, p = reader.readData( ('rho','p') )
    rbar = stats.reynolds_average(avg, r)
    pbar = stats.reynolds_average(avg, p)
    c2 = inp.gam*p/r
    r -= rbar
    p -= pbar
    return r,p,c2

def get_mom(reader,avg):
    r, u, v, w = reader.readData( ('rho','u','v','w') )
    ru = r*u - stats.reynolds_average(avg,r*u)
    rv = r*v - stats.reynolds_average(avg,r*v)
    rw = r*w - stats.reynolds_average(avg,r*w)
    return ru,rv,rw

# flist, tlist are lists of np arrays
def ddt(flist,tlist,method='fwd'):
    if method=='fwd':
        dfdt = (flist[1]-flist[0])/(tlist[1]-tlist[0])
    return dfdt

def solve(fhathat):
    # 2. Solve for particular soln and derivative, phi_L, dphi_L, zero means
    if rank==0: print('Computing scalar phi fft...')

    k2 = (kx**2 + ky**2 + kz**2)
    if rank_0x and rank_0y and rank_0z: k2[0,0,0] = 1 #suppress warning
    phihathat  = -fhathat/k2
    dphihathat = 1j*ky*phihathat
    if rank_0x and rank_0y and rank_0z: 
        phihathat[0,0,0]  = 0
        dphihathat[0,0,0] = 0
    if rank_0x: dphihathat[0,:,:] = 0
    if rank_0z: dphihathat[:,:,0] = 0
    phihat_L  = ffto._ifftY(phihathat);
    dphihat_L = ffto._ifftY(dphihathat);
    fhathat = None
    k2 = None

    # 3. Coefficients for homogenous exponetial soln
    # Zero coefficients at zero wavenumber kx=kz=0
    if rank==0: print('Computing scalar coeffs...')
    aplus = np.zeros(np.shape(kx))
    aminus = np.zeros(np.shape(kx))
    comm.Barrier()
    if rank_0y: # Do on rank 0y (boundaries)
        kmat = (kx[:,0,:]**2 + kz[:,0,:]**2)**0.5;
        if rank_0x and rank_0z: kmat[0,0] = 1 #to suppress warning
        tmp = 0.5/kmat
        if rank_0x and rank_0z: tmp[0,0] = 0;
        aplus = -tmp*(kmat*phihat_L[:,0,:] + dphihat_L[:,0,:]);
        aminus = tmp*(kmat*phihat_L[:,0,:] - dphihat_L[:,0,:]);
    
    # broadcast a+,a- coeffs so each proc has the entire bottom array:
    if rank_0y:
        aplus_  = [np.array([blkID[0],blkID[2]]), aplus]
        aminus_ = [np.array([blkID[0],blkID[2]]), aminus]
    else:
        aplus_  = [np.array([None]), np.array([None])]
        aminus_ = [np.array([None]), np.array([None])]
    # bcast 
    ap0 = comm.allreduce(aplus_)
    am0 = comm.allreduce(aminus_)
    ap_ = [val for val in ap0 if val.all() is not None]
    am_ = [val for val in am0 if val.all() is not None]
    ap,am = {},{} #make into dict for easier access
    for i in range(len(ap_)/2):
        ind = ap_[i*2]
        ap[ind[0],ind[1]] = ap_[i*2+1]
    for i in range(len(am_)/2):
        ind = am_[i*2]
        am[ind[0],ind[1]] = am_[i*2+1]
    #comm.Barrier()
    aplus=None
    aminus=None
    aplus_=None
    aminus_=None
    ap_ = None
    am_ = None
    ap0 = None
    am0 = None

    # 4. Make total solution for phi: 
    if rank==0: print('Computing scalar phi ...')
    kmat   = (kx[:,0,:]**2 + kz[:,0,:]**2)**0.5
    aplus  = ap[blkID[0],blkID[2]]
    aminus = am[blkID[0],blkID[2]]
    phihat = phihat_L
    for j in range(szy):
        phihat[:,j,:] += aplus*np.exp(kmat*(y[0,j,0]-ffto.LY)) 
        phihat[:,j,:] += aminus*np.exp(-kmat*(y[0,j,0]+ffto.LY))
    kmat = None
    aplus = None
    aminus = None
    phihat_L = None
    dphihat_L = None
    phi = ffto._ifftZ(ffto._ifftX(phihat))
    phihat=None
    return phi

def gather_y_plane(comm,reader,dat2save,y,yslice=0):
    Nx,Ny,Nz = reader.domain_size
    blank = np.zeros([Nx,Nz])
    flag=((yslice<np.amax(y)) and (yslice>np.amin(y)))
    if flag:
        idy = np.argmin(abs(y[0,:,0]-yslice))
        lo = reader._full_chunk_lo
        hi = reader._full_chunk_hi
        blank[lo[0]:(hi[0]+1),lo[2]:(hi[2]+1)] = dat2save[:,idy,:]
    plane = comm.reduce(blank,root=0)
    return plane

def gather_z_plane(comm,reader,dat2save,zslice=1):
    blank = np.zeros([Nx,Ny])
    flag = ((zslice<np.amax(z)) and (zslice>np.amin(z)))
    #print('range {}-{} has {}? {}'.format(np.amin(z),np.amax(z),zslice,flag))
    if flag:
        idz = np.argmin(abs(z[0,0,:]-zslice))
        lo = reader._full_chunk_lo
        hi = reader._full_chunk_hi
        blank[lo[0]:(hi[0]+1),lo[1]:(hi[1]+1)] = dat2save[:,:,idz]
        #idx0 = blkID[0]*szx
        #idy0 = blkID[1]*szy
        #idx = min(Nx,(blkID[0]+1)*szx)
        #idy = min(Ny,(blkID[1]+1)*szy)
        #blank[idx0:idx,idy0:idy] = dat2save[:,:,idz]
    plane = comm.reduce(blank,root=0)
    return plane

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage:" 
        print "This script writes the full field of dilation for some times."
        print "tID_list: csv, no spaces, eg [0,1,5]" 
        print "  python {} <prefix> [tID_list]".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    tID0 = int(sys.argv[2])
    if len(sys.argv)>3: 
        tIDmax = int(sys.argv[3])
    else: 
        tIDmax = tID0+1
    
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
    reader = pdr.ParallelDataReader(comm, serial_reader)#,pdecomp=[2,2,1])
    avg = red.Reduction(reader.grid_partition, periodic_dimensions)
    steps = sorted(reader.steps)
    tIDmax = max(tIDmax,max(steps))

    lo = reader._full_chunk_lo
    hi = reader._full_chunk_hi
    print('{}: {},{}'.format(rank,lo,hi))

    
    # Set up compact derivative object w/ 10th order schemes
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)
    
    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=(rank==0))
    du = inp.du
    if rank==0: print("\tdu = {}".format(inp.du))
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0)) 
    xmin = 0.;      xmax = Lx
    ymin = -Ly/2.;  ymax = Ly/2.
    zmin = 0.;      zmax = Lz
    yplot = np.linspace(-Ly/2,Ly/2,int(Ny))

    chunkx = abs(x[0,0,0]-x[-1,0,0])
    chunky = abs(y[0,0,0]-y[0,-1,0])
    chunkz = abs(z[0,0,0]-z[0,0,-1])
    blkID = [int((np.amax(x)-xmin)/chunkx)-1,
            int((np.amax(y)-ymin)/chunky)-1,
            int((np.amax(z)-zmin)/chunkz)-1]

    # Get the grid partition information
    szx,szy,szz = np.shape(x)
    nblk_x = int(np.round(Nx/(szx-1)))
    nblk_y = int(np.round(Ny/(szy-1)))    
    nblk_z = int(np.round(Nz/(szz-1)))
    if rank==0: print("Processor decomp: {}x{}x{}".format(nblk_x,nblk_y,nblk_z))
    #if nblk_y>1: sys.exit()

    # Setup the fft object
    if rank==0: print('Setting up fft...')
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)
    ffto = PoissonSol(settings)
    szx, szy, szz = ffto.size_3d 
    kxmax = ffto.NX * np.pi/ffto.LX
    kymax = ffto.NY * np.pi/ffto.LY
    kzmax = ffto.NZ * np.pi/ffto.LZ
   
    # set up the writer
    writer = h5_writer(settings)

    # What are the max wavenumbers? Need to check if ranks have these
    # when we need to zero oddballs (half of max)
    rank_0x = (0 in ffto.kx)
    rank_0y = (0 in ffto.ky)
    rank_0z = (0 in ffto.kz)
    thresh = 1e-10
    rank_oddx = (abs(kxmax -np.amax(abs(ffto.kx)))<thresh)
    rank_oddy = (abs(kymax -np.amax(abs(ffto.ky)))<thresh)
    rank_oddz = (abs(kzmax -np.amax(abs(ffto.kz)))<thresh)

    if rank==0: print('Computing fluctuations at t0...')
    reader.step = tID0 
    t0 = reader.time
    r0,p0,c2_0 = get_fluct(reader,avg)
    
    phiDict = {}
    yplanesDict = {}
    zplanesDict = {}
    ylist = range(0,int(Ly/2.),10)
     
    #for tID in steps[2:]:#tid_list[1:]:
    tID = tID0
    while tID < tIDmax:
        tID += 1

        if rank==0: print('Computing fluctuations...')
        reader.step = tID 
        t = reader.time
        r,p,c2 = get_fluct(reader,avg)
        
        kx = np.tile(ffto.kx[:,np.newaxis,np.newaxis],[1,szy,szz])
        ky = np.tile(ffto.ky[np.newaxis,:,np.newaxis],[szx,1,szz])
        kz = np.tile(ffto.kz[np.newaxis,np.newaxis,:],[szx,szy,1])

        for name in ['total','acoustic']:
            
            if rank==0: print('Getting ddt terms')
            if  name=='total': 
                rhs = -ddt([r,r0],[t,t0])
                r0 = None
            elif name=='acoustic': 
                rhs = -ddt([p,p0],[t,t0])
                rhs /= c2_0
                p0 = None

            # 1. Get RHS
            if rank==0: print('Computing RHS...') 
            fhat = ffto._fftZ(ffto._fftX(rhs))
            fhathat = ffto._fftY(fhat);
            fhat,rhs = None,None
            # need to zero any oddballs or means here?

            # Solve for scalar field phi
            phi = np.real(solve(fhathat))
            phiDict[name] = np.real(phi)
            fhathat,phi = None,None

        # Save planes
        for name in ['acoustic','thermal','b']:

            if rank==0: print('Saving {} planes'.format(name))
            if name=='acoustic':
                fx,fy,fz = der.gradient(phiDict[name])
            elif name=='thermal':
                fx,fy,fz = der.gradient(phiDict['total']-phiDict['acoustic'], 
                        y_bc=y_bc)
            elif name=='b':
                # Solve for solenoidal part
                if rank==0: print('Computing solenoidal components...')
                reader.step = tID0
                fx,fy,fz = get_mom(reader,avg)
                phix,phiy,phiz = der.gradient(phiDict['total'])
                fx -= phix
                fy -= phiy
                fz -= phiz
                phix,phiy,phiz = None,None,None
                #phiDict = {}
            
            # Save z=0 plane for viz
            #if rank==0: print('Saving planes at z=0')
            zplane_fx = np.zeros([Nx,Ny])  
            zplane_fy = np.zeros([Nx,Ny])  
            zplane_fz = np.zeros([Nx,Ny])  
            comm.Barrier()
            zplane_fx = gather_z_plane(comm,reader,fx,zslice=10)
            zplane_fy = gather_z_plane(comm,reader,fy,zslice=10)
            zplane_fz = gather_z_plane(comm,reader,fz,zslice=10)
            comm.Barrier()
            if rank==0:
                zplanesDict[name] = {}
                zplanesDict[name]['x'] = zplane_fx
                zplanesDict[name]['y'] = zplane_fy
                zplanesDict[name]['z'] = zplane_fz

            # Save horizontal planes
            planes_fx = np.zeros([Nx,Nz,len(ylist)])  
            planes_fy = np.zeros([Nx,Nz,len(ylist)])  
            planes_fz = np.zeros([Nx,Nz,len(ylist)])  
            for idx,yslice in enumerate(ylist):
                #if rank==0: print('Saving planes at y={}'.format(yslice))
                planes_fx[:,:,idx] = gather_y_plane(comm,reader,fx,y,yslice=yslice)
                planes_fy[:,:,idx] = gather_y_plane(comm,reader,fy,y,yslice=yslice)
                planes_fz[:,:,idx] = gather_y_plane(comm,reader,fz,y,yslice=yslice)
            if rank==0:
                yplanesDict[name] = {}
                yplanesDict[name]['x'] = planes_fx 
                yplanesDict[name]['y'] = planes_fy 
                yplanesDict[name]['z'] = planes_fz 


        # write 3d file
        #writer.update(tID0,t0)
        #writer.saveFields(phiDict,filePrefix=dirname+'/momdecomp_') 
        if rank==0:
            print('Saving...')
            savename = dirname+'/momdecomp_%04d.pkl'%tID0
            output = open(savename, 'wb')
            dat = {}
            dat['yplanes'] = yplanesDict
            dat['zplanes'] = zplanesDict
            pickle.dump(dat, output)
            output.close()
            print('Wrote to {}'.format(savename))
         
        # reset fields
        tID0,t0,r0,p0,c2_0 = tID,t,r,p,c2
