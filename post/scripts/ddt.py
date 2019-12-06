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

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage:" 
        print "This script writes the full field of dilation for some times."
        print "tID_list: csv, no spaces, eg [0,1,5]" 
        print "  python {} <prefix> [tID_list]".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    tid_list = map(int, sys.argv[2].strip('[]').split(',')) 
    if len(tid_list)<2: 
        print('Need at least 2 steps')
        sys.exit()

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

    # set up the writer
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)
    writer = h5_writer(settings)

    def get_fluct(reader,avg):
        r, p = reader.readData( ('rho','p') )
        rbar = stats.reynolds_average(avg, r)
        pbar = stats.reynolds_average(avg, p)
        r -= rbar
        p -= pbar
        return r,p
    
    # flist, tlist are lists of np arrays
    def ddt(flist,tlist,method='fwd'):
        if method=='fwd':
            dfdt = (flist[1]-flist[0])/(tlist[1]-tlist[0])
        return dfdt

    if rank==0: print('Computing fluctuations at t0...')
    tID0 = tid_list[0]
    reader.step = tID0 
    t0 = reader.time
    r0,p0 = get_fluct(reader,avg)
    c2 = inp.gam*p0/r0

    qDict = {}
    for tID in tid_list[1:]:
        if rank==0: print('Computing fluctuations...')
        reader.step = tID 
        t = reader.time
        r,p = get_fluct(reader,avg)

        qDict['drdt'] = ddt([r,r0],[t,t0])
        qDict['dpdt'] = ddt([p,p0],[t,t0])
        qDict['c2'] = c2 
        
        # write 3d file
        if rank==0: print('Writing to h5 file...')
        writer.update(tID0,t0)
        writer.saveFields(qDict,filePrefix=dirname+'/ddt_') 
    
        # reset fields
        tID0,t0,r0,p0 = tID,t,r,p
        c2 = inp.gam*p0/r0
