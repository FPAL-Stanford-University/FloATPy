#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

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
debug = False 

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

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage:" 
        print "This script writes the full field of dilation for some times."
        print "tID_list: csv, no spaces, eg [0,1,5]" 
        print "  python {} <prefix> [tID_list]".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    tid_list = map(int, sys.argv[2].strip('[]').split(',')) 
    outputfile  = filename_prefix + "hhdecomp"
    
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
    if debug:
        comm.Barrier()
        print('rank {} is block ID {}'.format(rank,blkID))

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
    #print('kxmax {},{}'.format(ffto.kx[int(Nx/2)],np.pi/ffto.LX*ffto.NX))
    # 0...kmax, -kmax, -kmin, the max (oddball) is negative at index N/2
    # kxmax = 4.63932869 -4.67585883
    # kymax = 4.6515054  -4.67585883
   
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
    if debug:
        comm.Barrier()
        print('rank {}: k0 [{} {} {}], kmax [{} {} {}]'.format(
            rank,rank_0x,rank_0y,rank_0z,rank_oddx,rank_oddy,rank_oddz))

    for tID in tid_list:
        
        kx = np.tile(ffto.kx[:,np.newaxis,np.newaxis],[1,szy,szz])
        ky = np.tile(ffto.ky[np.newaxis,:,np.newaxis],[szx,1,szz])
        kz = np.tile(ffto.kz[np.newaxis,np.newaxis,:],[szx,szy,1])

        # vel fluct
        if rank==0: print('Computing fluctuations...')
        reader.step = tID
        r, u, v, w = reader.readData( ('rho','u','v','w') )
        rbar = stats.reynolds_average(avg, r)
        utilde = stats.favre_average(avg, r, u, rho_bar=rbar)
        u -= utilde
        u /= du
        v /= du
        w /= du
        r = None
        if inp.Mc>0.6:
            u = window_field(filename_prefix,yplot,inp.du,tID,u)
            v = window_field(filename_prefix,yplot,inp.du,tID,v)
            w = window_field(filename_prefix,yplot,inp.du,tID,w)
        
        # Get perturbations fft in x,z, zero the means
        if rank==0: print('Computing fluctuations fft...')
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
        if debug:
            tmp = mpi_max(comm,stats.reynolds_average(avg,uhat))
            if rank==0: print(tmp)
            tmp = mpi_max(comm,stats.reynolds_average(avg,vhat))
            if rank==0: print(tmp)
            tmp = mpi_max(comm,stats.reynolds_average(avg,what))
            if rank==0: print(tmp)
            #(0.5077095987465812-6.938893903907228e-17j)
            #(0.23412663621636276+8.326672684688674e-17j)
            #(0.32394727601018203+0j)

        # Get derivative ddy{vhat}, zero the oddball
        if rank==0: print('Computing ddy vhat ...')
        ddy_vhat_hat = 1j*ky*ffto._fftY(vhat)
        if rank_oddy: 
            oddy = np.argmax(abs(ffto.ky))
            ddy_vhat_hat[:,oddy,:] = 0
        ddy_vhat = ffto._ifftY(ddy_vhat_hat)
        if debug:
            #if rank_oddy: print("rank {} zeroed y oddball at ky[{}]={}".format(
            #    rank,oddy,ffto.ky[oddy]))
            tmp = mpi_max(comm,stats.reynolds_average(avg,ddy_vhat_hat))
            if rank==0: print(tmp)
            tmp = mpi_max(comm,stats.reynolds_average(avg,ddy_vhat))
            if rank==0: print(tmp)
            # 2.212896723832909-0.4901971997534521j
            # 0.23175621471702385-6.938893903907228e-17j
            # 0.23175621471702362-1.9181274109061897e-05j without zero oddball
        ddy_vhat_hat = None

        # 1. Get divergence RHS, zero oddballs
        if rank==0: print('Computing divergence fft...')
        fhat = 1j*(kx*uhat + kz*what) + ddy_vhat;
        if rank_oddx: 
            oddx = np.argmax(abs(ffto.kx)); 
            fhat[oddx,:,:] = 0
        if rank_oddz: 
            oddz = np.argmax(abs(ffto.kz)); 
            fhat[:,:,oddz] = 0
        fhathat = ffto._fftY(fhat);
        if debug:
            tmp = mpi_max(comm,stats.reynolds_average(avg,fhathat))
            if rank==0: print(tmp)
            # 0.2498590190640754-0.0951071062026113j
        fhat = None
        uhat, vhat, what = None, None, None


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
        if debug:
            tmp = mpi_max(comm,stats.reynolds_average(avg,phihat_L))
            if rank==0: print(tmp)
            tmp = mpi_max(comm,stats.reynolds_average(avg,dphihat_L))
            if rank==0: print(tmp)
            # 0.12184347132611192+2.6020852139652106e-17j
            # 0.019613656482489538-8.778941053593642e-06j
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
        if debug: 
            tmp = mpi_max(comm,stats.reynolds_average(avg,abs(aplus)))
            if rank==0: print(tmp)
            tmp = mpi_max(comm,stats.reynolds_average(avg,abs(aminus)))
            if rank==0: print(tmp)
            #0.04322990287105939
            #0.036996887750879046
        
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
        if debug:
            print("rank {} will access coeff array sized {}".format(rank,np.shape(ap[blkID[0],blkID[2]])))
        comm.Barrier()
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
        phi = ffto._ifftZ(ffto._ifftX(phihat))
        if debug:
            tmp = mpi_max(comm,stats.reynolds_average(avg,abs(phi)))
            if rank==0: print(tmp)
            # 0.06352658697119544
        kmat = None
        aplus = None
        aminus = None
        phihat = None
        phihat_L = None
        dphihat_L = None

        # 5. Get velocities
        if rank==0: print('Computing ud, us...')
        dphihat = 1j*kx*ffto._fftX(phi)
        if rank_oddx: 
            oddx = np.argmax(abs(ffto.kx)); 
            dphihat[oddx,:,:] = 0
        ud = np.real(ffto._ifftX(dphihat))
        
        dphihat = 1j*ky*ffto._fftY(phi)
        if rank_oddy: 
            oddy = np.argmax(abs(ffto.ky)); 
            dphihat[:,oddy,:] = 0
        vd = np.real(ffto._ifftY(dphihat))

        dphihat = 1j*kz*ffto._fftZ(phi)
        if rank_oddz: 
            oddz = np.argmax(abs(ffto.kz)); 
            dphihat[:,:,oddz] = 0
        wd = np.real(ffto._ifftZ(dphihat))
        
        kx, ky, kz = None,None,None
        phi = None
        dphihat = None
        
        # Collect a snapshot along centerline
        blank = np.zeros([Nx,Nz])
        if rank_oddy:
            idy = np.argmin(abs(y[0,:,0]))
            idx0 = blkID[0]*szx
            idz0 = blkID[2]*szz
            idx = min(Nx,(blkID[0]+1)*szx)
            idz = min(Nz,(blkID[2]+1)*szz)
            blank[idx0:idx,idz0:idz] = ud[:,idy,:]
        plane = comm.reduce(blank,root=0)
        if rank==0:
            outputfile = filename_prefix + 'hhdecomp_%04d_xz'%tID
            np.save(outputfile, plane)
            print('Saved xz centerline to %s'%outputfile)

        us = u - ud
        vs = v - vd
        ws = w - wd
        uu   = stats.reynolds_average(avg,u**2)
        vv   = stats.reynolds_average(avg,v**2)
        uv   = stats.reynolds_average(avg,u*v)
        u,v,w = None,None,None
        
        # write 3d file
        """
        if rank==0: print('Writing to h5 file...')
        qDict = {}
        qDict['us'] = us
        qDict['ud'] = ud
        qDict['vs'] = vs
        qDict['vd'] = vd
        qDict['ws'] = ws
        qDict['wd'] = wd
        writer.update(tID,reader.time)
        writer.saveFields(qDict) 
        qDict = None
        """

        # Write 1d profiles
        udud = stats.reynolds_average(avg,ud*ud)
        vdvd = stats.reynolds_average(avg,vd*vd)
        udvd = stats.reynolds_average(avg,ud*vd)
        usud = stats.reynolds_average(avg,us*ud)*2
        vsvd = stats.reynolds_average(avg,vs*vd)*2
        cross = stats.reynolds_average(avg,us*vd+ud*vs)
        ud,vd,wd = None,None,None
        us,vs,ws = None,None,None

        # Normalize
        Rij = np.zeros([szy,9],dtype='f')
        Rij[:,0] = np.squeeze(uu)  
        Rij[:,1] = np.squeeze(vv)   
        Rij[:,2] = np.squeeze(uv)  
        Rij[:,3] = np.squeeze(udud) 
        Rij[:,4] = np.squeeze(vdvd)
        Rij[:,5] = np.squeeze(udvd) 
        Rij[:,6] = np.squeeze(usud) 
        Rij[:,7] = np.squeeze(vsvd) 
        Rij[:,8] = np.squeeze(cross)
        

        # Print stats
        if debug:
            if rank==0:
                scale = 1e4
                print("\n%0.3f & %0.3f & %0.3f & %0.3f & %0.3f & %0.3f & %0.3f & %0.3f & %0.3f"
                %(uu.max(),udud.max()*scale,usud.max()*scale,
                vv.max(),vdvd.max()*scale,vsvd.min()*scale,
                uv.min(),udvd.min()*scale,cross.max()*scale))
                #0.011 & 0.057 & 0.281 & 0.006 & 0.040 & -1.024 -0.004 & -0.025 & 0.928 # matlab
                #0.011 & 0.061 & 0.297 & 0.006 & 0.043 & -1.077 & -0.004 & -0.025 & 0.959 #serial

        # now gather from all processes into another array

        recvbuf = {}
        for i in range(9):
            recvbuf[i]=comm.gather(Rij[:,i], root=0)
            recvbuf[i] = np.array(recvbuf[i])

        if rank==0:
            total_array = np.zeros([Ny,9],dtype='f')
            for j in range(9):
                vec_array = np.reshape(recvbuf[j],[nblk_x,nblk_y,nblk_z,szy],order='F')
            
                # Now concat one column
                vec = vec_array[0,0,0,:];
                if (nblk_y>1):
                    for jj in range(1,nblk_y):
                        mat = np.squeeze(vec_array[0,jj,0,:])
                        vec = np.hstack([vec,mat])

                total_array[:,j] = vec
                print("i={}\t{}".format(j,np.sum(vec)))
           
            outputfile = filename_prefix + 'hhdecomp_%04d.dat'%tID
            np.savetxt(outputfile,total_array,delimiter=' ')
            print("{}".format(outputfile))
