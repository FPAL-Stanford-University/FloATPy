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

xdir = 0
zdir = 2
  
# Transpose from 3d to y
def transpose2y(numSet,q3d): 
    bufy = np.empty(numSet.chunk_y_size, 
            dtype=np.dtype(np.float64, align=True), order='F')
    numSet.grid_partition.transpose_3d_to_y(q3d,bufy)
    q3d = None
    return bufy

# Get 99% thickness
def get_L99(y,utilde):
    utilde = np.squeeze(utilde)
    du = abs(utilde[-1]-utilde[0])
    utop = 0.99*du/2.
    ubot = -0.99*du/2.
    Ny = np.size(y)
    ibot = np.argmin(abs(utilde[Ny/2:]-utop)[::-1])
    itop = np.argmin(abs(utilde[:Ny/2]-ubot)[::-1])+Ny/2
    #itop = np.argmin(abs(utilde[Ny/2:]-utop))+Ny/2
    #ibot = np.argmin(abs(utilde[:Ny/2]-ubot))
    L99 = y[itop]-y[ibot]
    if L99<0: print('utilde or y misoriented. exiting'); sys.exit()
    return L99, itop, ibot

def get_qpp(reader,avg,varname):
    r,q = reader.readData(('rho',varname))
    qtilde = stats.favre_average(avg,r,q)
    return q-qtilde[None,:,None]

def get_qp(reader,avg,varname):
    q = reader.readData(varname)
    qbar = stats.reynolds_average(avg,q[0])
    return q-qbar[None,:,None]

def get_corr(avg,qpp,yplot,y0):
    Ny = np.size(yplot)

    # center and bounding indices
    i0 = np.argmin(abs(yplot-y0))
    imax = min(abs(Ny-i0),i0)
    
    # Get numerator
    tmp = np.zeros(np.shape(qpp))
    for i in range(imax):
        tmp[:,i0+i,:] = qpp[:,i0-i,:]*qpp[:,i0+i,:]
        tmp[:,i0-i,:] = tmp[:,i0+i,:]

    # Normalize by mean value at y0 then avg
    num = stats.reynolds_average(avg,tmp)
    tmp = stats.reynolds_average(avg,qpp*qpp)
    denom = tmp[0,i0,0]
    corr = num[0,:,0]/denom
    return corr

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print "Usage:" 
        print "Computes the integral correlation lengths of a primitive var" 
        print "  python {} <prefix> [tID_list (csv)] varname ".format(sys.argv[0])
        sys.exit()
    if len(sys.argv) > 3:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
        varname  = sys.argv[3] 
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
    if rank==0: print("\tdu = {}".format(inp.du))
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
        qpp = get_qp(reader,avg,varname)
        qpp = transpose2y(settings,qpp) 

        # Get utilde
        fname = filename_prefix+'utilde_%04d.dat'%tID 
        utilde = np.fromfile(fname,count=-1,sep=' ')

        # get thicknesses
        L99,itop,ibot = get_L99(yplot,utilde)
        dtheta = get_dtheta(dirname,reader.time)

        # Start from y0
        if inp.rr==1:
            yc = 0
        else:
            ic = np.argmin(abs(utilde))
            yc = yplot[ic]
        offset = L99/4.#dtheta/2.
        y0_list = [yc,yc+offset,yc-offset]
        corr = np.zeros([Ny,len(y0_list)])
        for j,y0 in enumerate(y0_list):
            corr[:,j] = get_corr(avg,qpp,yplot,y0)

        if rank==0: 
            dir_out = dirname 
            if varname=='rho':
                outputfile = dir_out+"lscale_rr_%04d.dat"%tID
            else:
                outputfile = dir_out+"lscale_"+varname+varname+"_%04d.dat"%tID
            np.savetxt(outputfile,corr,delimiter=' ')
            print("Done writing to {}".format(outputfile))
