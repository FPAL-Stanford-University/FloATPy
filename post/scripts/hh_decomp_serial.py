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

# debug printing?
debug = True

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
    outputfile  = filename_prefix + "hhdecomp.dat"
    
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
    inp = nml.inputs(dirname,verbose=True)
    du = inp.du
    print("\tdu = {}".format(inp.du))

    # Setup the fft object
    print('Setting up fft...')
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0)) 
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)
    ffto = PoissonSol(settings)
    szx, szy, szz = ffto.size_3d 
    kx = np.tile(ffto.kx[:,np.newaxis,np.newaxis],[1,szy,szz])
    ky = np.tile(ffto.ky[np.newaxis,:,np.newaxis],[szx,1,szz])
    kz = np.tile(ffto.kz[np.newaxis,np.newaxis,:],[szx,szy,1])
    #kxmax = ffto.NX * np.pi/ffto.LX
    #kymax = ffto.NY * np.pi/ffto.LY
    #kzmax = ffto.NZ * np.pi/ffto.LZ
    #print(ffto.kx[Nx/2-1:Nx/2+1],kxmax); 
    #print(ffto.ky[Ny/2-1:Ny/2+1],kymax); 
    #print(ffto.kz[Nz/2-1:Nz/2+1],kzmax); 
    #sys.exit()

    # vel fluct
    print('Computing fluctuations...')
    reader.step = tid_list[0]
    r, u, v, w = reader.readData( ('rho','u','v','w') )
    rbar = stats.reynolds_average(avg, r)
    utilde = stats.favre_average(avg, r, u, rho_bar=rbar)
    u -= utilde

    # Get perturbations fft in x,z, zero the means
    print('Computing fluctuations fft...')
    uhat = ffto._fftZ(ffto._fftX(u))
    vhat = ffto._fftZ(ffto._fftX(v))
    what = ffto._fftZ(ffto._fftX(w))
    uhat[0,:,:] = 0
    vhat[0,:,:] = 0
    what[0,:,:] = 0
    uhat[:,:,0] = 0
    vhat[:,:,0] = 0
    what[:,:,0] = 0
    if debug:
        tmp = stats.reynolds_average(avg,uhat)
        print(np.amax(tmp))
        tmp = stats.reynolds_average(avg,vhat)
        print(np.amax(tmp))
        tmp = stats.reynolds_average(avg,what)
        print(np.amax(tmp))
        #(0.5077095987465812-6.938893903907228e-17j)
        #(0.23412663621636276+8.326672684688674e-17j)
        #(0.32394727601018203+0j)

    # Get derivative ddy{vhat}, zero the oddball
    print('Computing ddy vhat ...')
    ddy_vhat_hat = 1j*ky*ffto._fftY(vhat)
    ddy_vhat_hat[:,Ny/2,:] = 0
    ddy_vhat = ffto._ifftY(ddy_vhat_hat)
    if debug:
        tmp = stats.reynolds_average(avg,ddy_vhat_hat)
        print(np.amax(tmp))
        tmp = stats.reynolds_average(avg,ddy_vhat)
        print(np.amax(tmp))
        # 2.212896723832909-0.4901971997534521j
        # 0.23175621471702385-6.938893903907228e-17j

    # 1. Get divergence RHS, zero oddballs
    print('Computing divergence fft...')
    fhat = 1j*(kx*uhat + kz*what) + ddy_vhat;
    fhat[Nx/2,:,:] = 0
    fhat[:,:,Nz/2] = 0
    fhathat = ffto._fftY(fhat);
    if debug:
        tmp = stats.reynolds_average(avg,fhathat)
        print(np.amax(tmp))
        # 0.24985901906407545-0.09510710620261113j

    # 2. Solve for particular soln and derivative, phi_L, dphi_L, zero means
    print('Computing scalar phi fft...')
    k2 = (kx**2 + ky**2 + kz**2)
    k2[0,0,0] = 1 #suppress warning
    phihathat  = -fhathat/k2
    dphihathat = 1j*ky*phihathat
    phihathat[0,0,0]  = 0
    dphihathat[0,0,0] = 0
    dphihathat[0,:,:] = 0
    dphihathat[:,:,0] = 0
    phihat_L  = ffto._ifftY(phihathat);
    dphihat_L = ffto._ifftY(dphihathat);
    if debug:
        tmp = stats.reynolds_average(avg,phihat_L)
        print(np.amax(tmp))
        tmp = stats.reynolds_average(avg,dphihat_L)
        print(np.amax(tmp))
        # 0.12184347132611192+2.6020852139652106e-17j
        # 0.019613656482489538-8.778941053593642e-06j

    # 3. Coefficients for homogenous exponetial soln
    # Zero coefficients at zero wavenumber kx=kz=0
    if rank==0: print('Computing scalar coeffs...')
    phiL  = phihat_L[:,0,:];
    dphiL = dphihat_L[:,0,:];
    kmat = (kx[:,0,:]**2 + kz[:,0,:]**2)**0.5;
    kmat[0,0] = 1 #to suppress warning
    tmp = 0.5/kmat
    tmp[0,0] = 0;
    a_plus = -tmp*(kmat*phiL + dphiL);
    a_minus = tmp*(kmat*phiL - dphiL);
    aplus  = np.tile(a_plus[:,np.newaxis,:], [1,1,1])
    aminus = np.tile(a_minus[:,np.newaxis,:],[1,1,1])
    if debug:
        tmp = stats.reynolds_average(avg,np.tile(abs(phiL)[:,np.newaxis,:], [1,1,1]))
        print(np.amax(tmp))
        tmp = stats.reynolds_average(avg,abs(aplus))
        print(np.amax(tmp))
        tmp = stats.reynolds_average(avg,abs(aminus))
        print(np.amax(tmp))
        #0.06407680256581094
        #0.04322990287105939
        #0.036996887750879046

    
    # 4. Make total solution for phi
    print('Computing scalar phi ...')
    kmat   = (kx**2 + kz**2)**0.5 
    phihat = phihat_L 
    phihat += aplus*np.exp(kmat*(y-ffto.LY)) 
    phihat += aminus*np.exp(-kmat*(y+ffto.LY))
    phi = ffto._ifftZ(ffto._ifftX(phihat))
    if debug:
        tmp = stats.reynolds_average(avg,abs(phi))
        print(np.amax(tmp))
        # 0.06352658697119544

    # 5. Get velocities
    print('Computing ud, us...')
    dphihat = 1j*kx*ffto._fftX(phi)
    dphihat[Nx/2,:,:] = 0
    ud = np.real(ffto._ifftX(dphihat))
    
    dphihat = 1j*ky*ffto._fftY(phi)
    dphihat[:,Ny/2,:] = 0
    vd = np.real(ffto._ifftY(dphihat))

    dphihat = 1j*kz*ffto._fftZ(phi)
    dphihat[:,:,Nz/2] = 0
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
        #0.011 & 0.061 & 0.297 & 0.006 & 0.043 & -1.077 & -0.004 & -0.025 & 0.959 #serial, 
        # Mc12 256x384x128 with LAD shearlayer_0025
