#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys
from shutil import copyfile
#sys.path.insert(0,'/home/kmatsuno/h5py/build/lib.linux-x86_64-2.7/')
import h5py

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
    
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage: "
        print "  python {} <dir_read> [tID_read] <dir_write> ".format(sys.argv[0])
        sys.exit()
    dir_in  = sys.argv[1]
    tID     = int(sys.argv[2])
    dir_out = sys.argv[3]
    fname_in = dir_in + '/restart_%04d'%tID + '.h5' 
    fname_out = dir_out + '/seed.h5' 
    
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
    reader = por.PadeopsReader(dir_in+'/restart_',
            periodic_dimensions=periodic_dimensions)

    # setup the inputs object
    inp_s = nml.inputs(dir_in,verbose=False)
    inp   = nml.inputs(dir_out,verbose=False)
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dir_out,verbose=(rank==0)) 
    if rank==0: 
        print("\tSeed Mc,rr,du = {},{},{}".format(inp_s.Mc,inp_s.rr,inp_s.du))
        print("\tNew  Mc,rr,du = {},{},{}".format(inp.Mc,inp.rr,inp.du))
    
    #  setup seed file
    seedfile = h5py.File(fname_out, 'w' )
    zero = np.zeros(1)
    seedfile.attrs.create('Time',zero.astype('>d'))
    seedfile.attrs.create('step',zero.astype('>i4'))
    shape = np.array([Nz,Ny,Nx])
    sz = (4,shape[1],4) 

    # Read input restart file, copy these over
    reader.step = tID
   
    print("Reading density")
    r1, r2  = reader.readData( ('rhoY_0001','rhoY_0002') )
    r = r1+r2; 
    r1 = np.transpose(r1)
    r2 = np.transpose(r2)
    print("Writing out r1,r2")
    r1dat = seedfile.create_dataset('rhoY_0001', shape,data=r1, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    r2dat = seedfile.create_dataset('rhoY_0002', shape,data=r2, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    r1,r2 = None, None
    r1dat,r2dat = None, None
    
    # Rescale u so the freestream values are correct
    print("Updating velocity and total energy")
    ru,TE = reader.readData(('rhou','TE'))
    u = ru/r; ru = None
    TE -= 0.5*r*u*u
    u *= inp.du/inp_s.du
    TE += 0.5*r*u*u 
    ru = r*u 
    r,u = None,None
    print("Writing out u,TE")
    ru = np.transpose(ru)
    TE = np.transpose(TE)
    rudat = seedfile.create_dataset('rhou', shape,data=ru, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    TEdat = seedfile.create_dataset('TE', shape,data=TE, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    ru,TE = None,None
    rudat,TEdat = None,None

    # write the new field out size 808094096
    rv,rw = reader.readData(('rhov','rhow'))
    rv = np.transpose(rv)
    rw = np.transpose(rw)
    print("Writing out rv,rw")
    rvdat = seedfile.create_dataset('rhov', shape,data=rv, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    rwdat = seedfile.create_dataset('rhow', shape,data=rw, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    rv,rw = None,None
    rvdat,rwdat = None,None
    seedfile.close()

    print('Done')
