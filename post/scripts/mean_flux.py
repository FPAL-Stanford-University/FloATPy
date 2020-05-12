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

xdir = 0
zdir = 2

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:" 
        print "Computes ubar,vbar,utilde,vtilde, r'u' and r'v'" 
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
    gp = reader.grid_partition

    # Set up compact derivative object w/ 10th order schemes
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(gp, (dx, dy, dz), (10, 10, 10), periodic_dimensions)
    
    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=(rank==0))
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    du = inp.du
    if rank==0: print("\tdu = {}".format(inp.du))
    yplot = np.linspace(Ly/2,-Ly/2,int(Ny))
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)
    
    # grid indicies accessed by this proc
    lo = reader._full_chunk_lo
    hi = reader._full_chunk_hi

    #
    tID_list = steps
    for tID in tID_list:
        reader.step = tID
        r,u,v,w,p = reader.readData( ('rho','u','v','w','p') )
        e = p/r/(inp.gam-1) 
       
        # Get utilde and L99
        fname = filename_prefix+'utilde_%04d.dat'%tID 
        utilde = np.fromfile(fname,count=-1,sep=' ')
        L99,ibot,itop = get_L99(yplot,utilde)
    
        # masks
        tmp = np.ones(np.shape(x))
        masktop,maskbot,maskmix = tmp,tmp,tmp
        assert(yplot[itop]>0)
        assert(yplot[ibot]<0)
        masktop[:,itop:,:] = 1.
        maskbot[:,:ibot,:] = 1.
        maskmix[:,ibot:itop,:] = 1.

        # mass, momenta and energy flux
        F_mass = stats.reynolds_average(avg,transpose2y(settings,r*v))
        F_momx = stats.reynolds_average(avg,transpose2y(settings,r*u*v))
        F_momy = stats.reynolds_average(avg,transpose2y(settings,r*v*v))
        F_momz = stats.reynolds_average(avg,transpose2y(settings,r*w*v))
        F_e    = stats.reynolds_average(avg,transpose2y(settings,r*v*e))
            
        # mass, momenta and energy internal generation
        G_mass = stats.integrate_volume(r,  dx,dy,dz,gp,mask=maskmix)
        G_momx = stats.integrate_volume(r*u,dx,dy,dz,gp,mask=maskmix)
        G_momy = stats.integrate_volume(r*v,dx,dy,dz,gp,mask=maskmix)
        G_momz = stats.integrate_volume(r*w,dx,dy,dz,gp,mask=maskmix)
        G_e    = stats.integrate_volume(r*e,dx,dy,dz,gp,mask=maskmix)

        # Output profiles to file
        if rank==0: 
            print('Writing to file...')
            dir_out = dirname.split('/projects/ShockInducedMix/')[-1]
            dir_out = '/home/kmatsuno/' + dir_out + '/'
            outputfile = dir_out+ '/mean_flux_%04d.h5'%tID
            hf = h5py.File(outputfile, 'w')
            hf.create_dataset('itop',   data=itop)
            hf.create_dataset('ibot',   data=ibot)
            hf.create_dataset('F_mass', data=np.squeeze(F_mass))
            hf.create_dataset('F_momx', data=np.squeeze(F_momx))
            hf.create_dataset('F_momy', data=np.squeeze(F_momy))
            hf.create_dataset('F_momz', data=np.squeeze(F_momz))
            hf.create_dataset('F_e',    data=np.squeeze(F_e))
            hf.create_dataset('G_mass', data=G_mass)
            hf.create_dataset('G_momx', data=G_momx)
            hf.create_dataset('G_momy', data=G_momy)
            hf.create_dataset('G_momz', data=G_momz)
            hf.create_dataset('G_e',    data=G_e )
            hf.close()
            print('Saved to %s'%outputfile)
        
