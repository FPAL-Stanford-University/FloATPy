#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys
import h5py
from shutil import copyfile

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
    fname_out = dir_out + '/restart_0000.h5' 
    
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
    try:
        serial_reader = por.PadeopsReader(dir_in+'/restart_',
            periodic_dimensions=periodic_dimensions)
    except:
        print('Linking coords...')
        os.symlink(dir_in + '/shearlayer_coords.h5', dir_in+'/restart_coords.h5')
        serial_reader = por.PadeopsReader(dir_in+'/restart_',
            periodic_dimensions=periodic_dimensions)
    reader = pdr.ParallelDataReader(comm, serial_reader)
    avg = red.Reduction(reader.grid_partition, periodic_dimensions)
    steps = sorted(reader.steps)

    # Set up the derivative object
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)

    # setup the inputs object
    verbosity = False
    if rank == 0: verbosity = True
    inp_s = nml.inputs(dir_in,verbose=verbosity)
    inp   = nml.inputs(dir_out,verbose=verbosity)
    if rank==0: 
        print("\tSeed Mc,rr,du = {},{},{}".format(inp_s.Mc,inp_s.rr,inp_s.du))
        print("\tNew  Mc,rr,du = {},{},{}".format(inp.Mc,inp.rr,inp.du))

    # Read input restart file, copy these over
    reader.step = tID
    qlist = ('rhou', 'rhov', 'rhow', 'TE', 'rhoY_0001','rhoY_0002')
    ru, rv, rw, TE, r1, r2  = reader.readData( qlist )
    
    # Subtract the KE of the initial u velocity:
    # TE = r*(e + 0.5*uu)
    print("Getting primitives")
    r = r1+r2
    u = ru/r
    v = rv/r
    w = rw/r
    TE -= 0.5*r*u*u

    # Rescale u so the freestream values are correct
    u *= inp.du/inp_s.du
    
    # Rescale r1 so the freestream values are correct
    #def get_rho_vals(inputs):
    #    lmbda = (inputs.rr-1.)/(inputs.rr+1.)
    #    r1_top = (1.+lmbda)/(1+inputs.rr)
    #    r2_top = (1.+lmbda)/(1.+1./inputs.rr)
    #    r1_bot = (1.-lmbda)/(1+inputs.rr)
    #    r2_bot = (1.-lmbda)/(1.+1./inputs.rr)
    #    return r1_top, r2_top, r1_bot, r2_bot
    #def rescale_rho(r,rmax_old,rmin_old,rmax,rmin):
    #    r -= rmin_old #subtract -infty to zero
    #    r /= rmax_old #scale +infty to 1
    #    r *= (rmax-rmin) #scale to new range
    #    r += rmin    # scale to new -infty
    #    print(r[0,0,0],rmin)
    #    print(r[0,-1,0],rmax)
    #    return r
    #r1_top0, r2_top0, r1_bot0, r2_bot0 = get_rho_vals(inp_s)
    #r1_top, r2_top, r1_bot, r2_bot = get_rho_vals(inp)
    #r1 = rescale_rho(r1,r1_top0,r1_bot0,r1_top,r1_bot); print('r1 rescaled')
    #r2 = rescale_rho(r2,r2_top0,r2_bot0,r2_top,r2_bot); print('r2 rescaled')
    #r = r1+r2;
    #print(r[0,-1,0]/r[0,0,0])
    #sys.exit()

    # Add back to total energy
    print("Remaking convservatives")
    TE += 0.5*r*u*u 
    ru = r*u

    # Write attributes time and step
    seedfile = h5py.File(fname_out, 'w' )
    zero = np.zeros(1)
    seedfile.attrs.create('Time',zero.astype('>d'))
    seedfile.attrs.create('step',zero.astype('>i4'))

    # write the new field out size 808094096
    print("Transposing")
    ru = np.transpose(ru)
    rv = np.transpose(rv)
    rw = np.transpose(rw)
    TE = np.transpose(TE)
    r1 = np.transpose(r1)
    r2 = np.transpose(r2)
    shape = np.shape(ru)
    sz = (2,shape[1],2) 
    print("Writing out")
    rudat = seedfile.create_dataset(qlist[0], shape,data=ru, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    rvdat = seedfile.create_dataset(qlist[1], shape,data=rv, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    rwdat = seedfile.create_dataset(qlist[2], shape,data=rw, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    TEdat = seedfile.create_dataset(qlist[3], shape,data=TE, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    r1dat = seedfile.create_dataset(qlist[4], shape,data=r1, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)
    r2dat = seedfile.create_dataset(qlist[5], shape,data=r2, dtype=np.dtype('>d'),chunks=sz, fillvalue=0)

    seedfile.close()
