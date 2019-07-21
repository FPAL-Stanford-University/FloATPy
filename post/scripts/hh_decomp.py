from mpi4py import MPI
import numpy as np
import os
import sys
import scipy.io
import matplotlib.pyplot as plt

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml
from SettingLib import NumSetting
from PoissonSol import * 

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage:" 
        print "This script writes the full field of dilation for some times."
        print "tID_list: csv, no spaces, eg [0,1,5]" 
        print "  python {} <prefix> [tID_list]".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    tid_list = map(int, sys.argv[2].strip('[]').split(',')) 
    outputfile  = filename_prefix + "enstrophy_spectrum.dat"
    
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
    verbose = False
    if rank==0: verbose=True
    inp = nml.inputs(dirname,verbose=verbose)
    du = inp.du
    if rank==0: print("\tdu = {}".format(inp.du))

    # Setup the fft object
    if rank==0: print('Setting up fft...')
    Nx,Ny,Nz = reader.domain_size
    settings = NumSetting( NX=Nx, NY=Ny, NZ=Nz,
                 XMIN=x[0,0,0], XMAX=x[-1,0,0]+dx,
                 YMIN=y[0,0,0], YMAX=y[0,-1,0]+dy,
                 ZMIN=z[0,0,0], ZMAX=z[0,0,-1]+dz,
                 order=10)
    ffto = PoissonSol(settings)
    szx, szy, szz = ffto.size_3d 
    kx = np.tile(ffto.kx[:,np.newaxis,np.newaxis],[1,szy,szz])
    ky = np.tile(ffto.ky[np.newaxis,:,np.newaxis],[szx,1,szz])
    kz = np.tile(ffto.kz[np.newaxis,np.newaxis,:],[szx,szy,1])

    if rank==1:
        print(settings.chunk_3d_size)
        print(settings.chunk_x_size)
        print(settings.chunk_y_size)
        print(settings.chunk_z_size)
        print(settings.chunk_z_lo)
        print(settings.chunk_z_hi)

    # Which ranks have idx=0, N/2?
    rank_0x = 0 
    rank_0y = 0 
    rank_0z = 0 
    rank_nx = 0 
    rank_ny = 0 
    rank_nz = 0 
    if 0 in kx: rank_0x = rank
    if 0 in ky: rank_0y = rank
    if 0 in kz: rank_0z = rank
    if 2*np.pi/ffto.LX in kx: rank_nx = rank
    if 2*np.pi/ffto.LY in ky: rank_ny = rank
    if 2*np.pi/ffto.LZ in kz: rank_nz = rank

    # vel fluct
    if rank==0: print('Computing fluctuations...')
    reader.step = tid_list[0]
    r, u, v, w = reader.readData( ('rho','u','v','w') )
    rbar = stats.reynolds_average(avg, r)
    utilde = stats.favre_average(avg, r, u, rho_bar=rbar)
    u -= utilde

    # Get perturbations fft in x,z, zero the means
    if rank==1: print('Computing fluctuations fft...')
    uhat = ffto._fftZ(ffto._fftX(u))
    vhat = ffto._fftZ(ffto._fftX(v))
    what = ffto._fftZ(ffto._fftX(w))
    if rank_0x:
        uhat[0,:,:] = 0
        vhat[0,:,:] = 0
        what[0,:,:] = 0
    if rank_0z:
        uhat[:,:,0] = 0
        vhat[:,:,0] = 0
        what[:,:,0] = 0

    # Get derivative ddy{vhat}, zero the oddball
    if rank==1: print('Computing ddy vhat ...')
    ddy_vhat_hat = 1j*ky*ffto._fftY(vhat)
    if rank_ny:
        ddy_vhat_hat[:,Ny/2,:] = 0
    ddy_vhat = ffto._ifftY(ddy_vhat_hat)

    # 1. Get divergence RHS, zero oddballs
    if rank==1: 
        print('Computing divergence fft...')
        print(np.shape(kx))
        print(np.shape(kz))
        print(np.shape(uhat))
        print(np.shape(what))
        print(np.shape(ddy_vhat))
    fhat = 1j*(kx*uhat + kz*what) + ddy_vhat;
    if rank_nx: fhat[Nx/2,:,:] = 0
    if rank_nz: fhat[:,:,Nz/2] = 0
    fhathat = ffto._fftY(fhat);

    # 2. Solve for particular soln and derivative, phi_L, dphi_L, zero means
    k2 = (kx**2 + ky**2 + kz**2)
    if rank_0x and rank_0y and rank_0z: k2[0,0,0] = 1 #suppress warning
    phihathat  = -fhathat/k2
    dphihathat = 1j*ky*phihathat
    if rank_0x and rank_0y and rank_0z: 
        phihathat[0,0,0]  = 0
        dphihathat[0,0,0] = 0
    if rank_0x: dphihathat[0,:,:] = 0
    if rank_0z: dphihathat[:,:,0] = 0
    if rank==1: print('Computing scalar phi fft...')
    phihat_L  = ffto._ifftY(phihathat);
    dphihat_L = ffto._ifftY(dphihathat);
    
    # 3. Coefficients for homogenous exponetial soln
    # Zero coefficients at zero wavenumber kx=kz=0
    #if rank_0y:
    phiL  = dphihat_L[:,0,:];
    dphiL = dphihat_L[:,0,:];
    kmat = (kx[:,0,:]**2 + kz[:,0,:]**2)**0.5;
    if rank_0x and rank_0z: kmat[0,0] = 1 #to suppress warning
    tmp = 0.5/kmat
    if rank_0x and rank_0z: tmp[0,0] = 0;
    a_plus = -tmp*(kmat*phiL + dphiL);
    a_minus = tmp*(kmat*phiL - dphiL);
    aplus  = np.tile(a_plus[:,np.newaxis,:], [1,szy,1])
    aminus = np.tile(a_minus[:,np.newaxis,:],[1,szy,1])
    
    # 4. Make total solution for phi
    if rank==1: print('Computing scalar phi ...')
    kmat   = (kx**2 + kz**2)**0.5 
    phihat = phihat_L 
    phihat += aplus*np.exp(kmat*(y-ffto.LY)) 
    phihat += aminus*np.exp(-kmat*(y+ffto.LY))
    tmp = ffto._ifftX(phihat)
    phi = ffto._ifftZ(tmp)

    # 5. Get velocities
    if rank==0: print('Computing ud, us...')
    dphihat = 1j*kx*ffto._fftX(phi)
    if rank_nx: dphihat[Nx/2,:,:] = 0
    ud = np.real(ffto._ifftX(dphihat))
    
    dphihat = 1j*ky*ffto._fftY(phi)
    if rank_ny: dphihat[:,Ny/2,:] = 0
    vd = np.real(ffto._ifftY(dphihat))

    dphihat = 1j*kz*ffto._fftZ(phi)
    if rank_nz: dphihat[:,:,Nz/2] = 0
    wd = np.real(ffto._ifftZ(dphihat))

    # Normalize
    u /= du; ud /= du
    v /= du; vd /= du
    w /= du; wd /= du
    us = u - ud
    vs = v - vd
    ws = w - wd
    
    scale = 1e4
    uu   = stats.reynolds_average(avg,u**2)
    vv   = stats.reynolds_average(avg,v**2)
    uv   = stats.reynolds_average(avg,u*v)
    udud = stats.reynolds_average(avg,ud*ud) * scale
    vdvd = stats.reynolds_average(avg,vd*vd) * scale
    udvd = stats.reynolds_average(avg,ud*vd) * scale
    usud = stats.reynolds_average(avg,us*ud)*2 * scale
    vsvd = stats.reynolds_average(avg,vs*vd)*2 * scale
    cross = stats.reynolds_average(avg,us*vd+ud*vs) * scale

    # Print stats
    if rank==0:
        print("\n%0.3f & %0.3f & %0.3f & %0.3f & %0.3f & %0.3f & %0.3f & %0.3f & %0.3f"
        %(uu.max(),udud.max(),usud.max(),
        vv.max(),vdvd.max(),vsvd.min(),
        uv.min(),udvd.min(),cross.max()))
        #0.011 & 0.057 & 0.281 & 0.006 & 0.040 & -1.024 -0.004 & -0.025 & 0.928 # matlab
        #0.011 & 0.061 & 0.297 & 0.006 & 0.043 & -1.077 & -0.004 & -0.025 & 0.959 #serial
