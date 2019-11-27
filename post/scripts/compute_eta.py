from mpi4py import MPI
import numpy as np
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml
from common import *

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <dirname> tID_tkeb tID ".format(sys.argv[0])
        sys.exit()
    dirname = sys.argv[1]
    tID_tkeb = int(sys.argv[2]) 
    tID = int(sys.argv[3]) 


    periodic_dimensions = (True,False,True)
    x_bc = (0,0)
    y_bc = (0,0)
    z_bc = (0,0)

    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    procs = comm.Get_size()

    # Set up TKE reader
    TKEreader = por.PadeopsReader(dirname+'TKEBudget_', 
            periodic_dimensions=periodic_dimensions)
    reader = por.PadeopsReader(dirname+'shearlayer_',
            periodic_dimensions=periodic_dimensions)
    
    # Get parameters
    inp   = nml.inputs(dirname,verbose=(rank==0))
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    y = np.linspace(Ly/2.,-Ly/2.,Ny)
    dy = Ly/Ny

    # Centerline slice:
    if inp.rr == 1: ic = Ny/2
    else: 
        ic,yc = get_centerline(dirname,y,tID)
        print('Using yc = {}'.format(yc))
    
    # Get dissipation
    print('Reading dissipation')
    TKEreader.step = tID_tkeb
    diss = TKEreader.readData(('dissipation'))
    eps = np.squeeze(np.array(diss[0]))
    eps = eps[ic]

    # Get nu
    print('Reading nu')
    reader.step = tID
    rho,mu = reader.readData(('rho','mu'))
    nu = np.array(mu)/np.array(rho)
    nu = np.mean(np.mean(nu,axis=2),axis=0) 
    nu = np.squeeze(nu)
    nu = nu[ic]

    # Compute eta
    eta = (nu**0.75)*(eps**-0.25)
    u_eta = (nu*eps)**0.25
    t_eta = (nu/eps)**0.5

    # Write to file
    nstats = 3
    dat = np.zeros([nstats],dtype='f')
    dat[0] = eta
    dat[1] = u_eta
    dat[2] = t_eta
    
    print(dat)
    print('Eta/dy = {}'.format(eta/dy))
    sys.exit()




    # now gather from all processes for each Rij into another array
    recvbuf = {}
    root=0
    recvbuf=comm.gather(dat, root=0)
    comm.Barrier()
    recvbuf = np.array(recvbuf)

    if rank==0:
        total_array = np.zeros([Ny,nstats],dtype='f')
        for j in range(6):
            vec_array = np.reshape(recvbuf[j],[nblk_x,nblk_y,nblk_z,ny],order='F')
        
            # Now concat one column
            vec = vec_array[0,0,0,:];
            if (nblk_y>1):
                for jj in range(1,nblk_y):
                    mat = np.squeeze(vec_array[0,jj,0,:])
                    vec = np.hstack([vec,mat])

            total_array[:,j] = vec
            print(vec)
            print("i={}\t{}".format(j,np.sum(vec)))
       
        print(total_array)
        outputfile = dirname + 'eta_%04d.dat'%step
        np.savetxt(outputfile,total_array,delimiter=' ')
        print("{}".format(outputfile))
        print("Min eta = {}".format(np.amin(eta)))
