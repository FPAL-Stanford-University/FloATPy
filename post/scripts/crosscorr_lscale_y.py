from mpi4py import MPI
import numpy as np
import numpy.linalg as la
import os
import sys
import matplotlib.pyplot as plt

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml
from SettingLib import NumSetting
from PoissonSol import * 
from common import *

def get_crosscorr(avg,q1,q2,yplot,y0):
    Ny = np.size(yplot)

    # center and bounding indices
    i0 = np.argmin(abs(yplot-y0))
    imax = min(abs(Ny-i0),i0)
    
    # Get numerator
    tmp = np.zeros(np.shape(q1))
    for i in range(imax):
        tmp[:,i0+i,:] = q1[:,i0,:]*q2[:,i0+i,:]
        tmp[:,i0-i,:] = q1[:,i0,:]*q2[:,i0-i,:]

    # Normalize by mean value at y0 then avg
    num = stats.reynolds_average(avg,tmp)
    tmp = stats.reynolds_average(avg,abs(q1)*abs(q2))
    denom = tmp[0,i0,0]
    corr = num[0,:,0]/denom
    return corr

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:" 
        print "Computes the corss correlation lengths of u' and v'' " 
        print "  python {} <prefix> [tID_list (csv)]".format(sys.argv[0])
        sys.exit()
    elif len(sys.argv)==2:
        tID_list = None
    elif len(sys.argv)==3:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
    filename_prefix = sys.argv[1]

    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    procs = comm.Get_size()
    
    # Set up the serial Miranda reader
    # Set up the parallel reader
    # Set up the reduction object
    periodic_dimensions = (True,False,True)
    serial_reader = por.PadeopsReader(filename_prefix, 
            periodic_dimensions=periodic_dimensions)
    reader = pdr.ParallelDataReader(comm, serial_reader)
    avg = red.Reduction(reader.grid_partition, periodic_dimensions)
    steps = sorted(reader.steps)
    if tID_list is None: tID_list = steps

    # Set up compact derivative object w/ 10th order schemes
    x, y, z = reader.readCoordinates()
    Nx,Ny,Nz = reader.domain_size
    dx,dy,dz = grid_res(x,y,z)
    lo = reader._full_chunk_lo
    hi = reader._full_chunk_hi
    print('{}: {},{}'.format(rank,lo,hi))
    
    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=(rank==0))
    du = inp.du
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    yplot = np.linspace(-Ly/2,Ly/2,int(Ny))
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)

    for tID in tID_list:
        reader.step = tID
        up = get_qp(reader,avg,'u')
        vp = get_qp(reader,avg,'v')
        vpp = get_qpp(reader,avg,'v')
        if procs > 1: 
            if rank==0: print('Transposing')
            up = transpose2y(settings,up) 
            vp = transpose2y(settings,vp) 
            vpp = transpose2y(settings,vpp) 

        # Get utilde
        try:
            fname = filename_prefix+'utilde_%04d.dat'%tID 
            utilde = np.fromfile(fname,count=-1,sep=' ')
        except:
            dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/mira/')[-1]
            dir_out = '/home/kmatsuno/ShearLayerData/production/' + dir_out + '/'
            fname = dir_out+'shearlayer_utilde_%04d.dat'%tID 
            utilde = np.fromfile(fname,count=-1,sep=' ')

        # get thicknesses
        L99,itop,ibot = get_L99(-yplot,utilde)
        dtheta = get_dtheta(dirname,reader.time)

        # Start from y0
        ic = np.argmin(abs(utilde))
        yc = yplot[ic]
        offset = L99/4.#dtheta/2.
        y0_list = [yc,yc+offset,yc-offset]
        corr1 = np.zeros([Ny,len(y0_list)])
        corr2 = np.zeros([Ny,len(y0_list)])
        for j,y0 in enumerate(y0_list):
            corr1[:,j] = get_crosscorr(avg,up,vp,yplot,y0)
            corr2[:,j] = get_crosscorr(avg,up,vpp,yplot,y0)

        if rank==0: 
            dir_out = dirname.split('/lus/theta-fs0/projects/HighMachTurbulence/ShearLayerData/mira/')[-1]
            dir_out = '/home/kmatsuno/ShearLayerData/production/' + dir_out + '/'
            dir_out = dirname 
            outputfile = dir_out+"crosscorr_upvp_%04d.dat"%tID
            np.savetxt(outputfile,corr1,delimiter=' ')
            print("Done writing to {}".format(outputfile))
            outputfile = dir_out+"crosscorr_upvpp_%04d.dat"%tID
            np.savetxt(outputfile,corr2,delimiter=' ')
            print("Done writing to {}".format(outputfile))
