#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys
from shutil import copyfile
sys.path.insert(0,'/home/kmatsuno/h5py/build/lib.linux-x86_64-2.7/')
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
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <fname> ".format(sys.argv[0])
        sys.exit()
    fname_in = sys.argv[1]
    dir_in = './'
    dir_out = './'
    fname_out = './transposed.h5' 
    
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
    #reader = por.PadeopsReader(dir_in+'/restart_',
    #        periodic_dimensions=periodic_dimensions)

    # setup the inputs object
    inp_s = nml.inputs(dir_in,verbose=False)
    inp   = nml.inputs(dir_out,verbose=False)
    if rank==0: 
        print("\tSeed Mc,rr,du = {},{},{}".format(inp_s.Mc,inp_s.rr,inp_s.du))
        print("\tNew  Mc,rr,du = {},{},{}".format(inp.Mc,inp.rr,inp.du))

    # Read input restart file, copy these over
    print("Reading data")
    hf = h5py.File(fname_in,'r') 
    qlist = ('rhou', 'rhov', 'rhow', 'TE', 'rhoY_0001','rhoY_0002')
    
    # Write attributes time and step
    seedfile = h5py.File(fname_out, 'w' )
    zero = np.zeros(1)
    seedfile.attrs.create('Time',zero.astype('>d'))
    seedfile.attrs.create('step',zero.astype('>i4'))

    # write the new field out size 808094096
    #print("Transposing")
    #ru = np.transpose(ru)
    #rv = np.transpose(rv)
    #rw = np.transpose(rw)
    #TE = np.transpose(TE)
    #r1 = np.transpose(r1)
    #r2 = np.transpose(r2)
    shape = np.shape(hf[qlist[0]]) 
    sz = (shape[0]/2,shape[1],shape[2]/2) 
    print("Writing out")
    for qname in qlist:
        print(qname)
        dat = seedfile.create_dataset(qname, shape,data=hf[qname], 
                dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    seedfile.close()
