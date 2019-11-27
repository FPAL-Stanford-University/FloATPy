from mpi4py import MPI
import numpy as np
import numpy.linalg as la
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml
from SettingLib import NumSetting
from PoissonSol import * 

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
   
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:" 
        print "Computes pressure strain from a restart file" 
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


    # Set up compact derivative object w/ 10th order schemes
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)
    
    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=(rank==0))
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    du = inp.du
    if rank==0: print("\tdu = {}".format(inp.du))

    # Get the grid partition information
    nx,ny,nz = reader._full_chunk_size
    szx,szy,szz = np.shape(x)
    nblk_x = int(np.round(Nx/(szx-1)))
    nblk_y = int(np.round(Ny/(szy-1)))    
    nblk_z = int(np.round(Nz/(szz-1)))
    if rank==0: print("Processor decomp: {}x{}x{}".format(nblk_x,nblk_y,nblk_z))
    
    # Compute stats at each step:
    i = 0
    thresh = 0.15
    for tid in tID_list:
        reader.step = tid
        r,u,v,p = reader.readData( ('rho','u','v','p') )
        rbar = stats.reynolds_average(avg,r)
        utilde = stats.favre_average(avg,r,u,rho_bar=rbar)
        vtilde = stats.favre_average(avg,r,v,rho_bar=rbar)
        upp = np.array(u-utilde)
        vpp = np.array(v-vtilde)
        p = np.array(p)

        # Compute production
        R12 = stats.reynolds_average(avg,r*upp*vpp)
        qx,qy,qz = der.gradient(u, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        dUdy = stats.favre_average(avg,r,qy,rho_bar=rbar)
        prod = R12*dUdy
        r,u,v = None,None,None

        # Copmute sij_pp
        s = {}
        qx,qy,qz = der.gradient(upp, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        s['11'] = qx
        s['12'] = qy/2.
        upp = None
        qx,qy,qz = der.gradient(vpp, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        s['22'] = qy
        s['12'] += qx/2.
        vpp = None
        qx,qy,qz = None,None,None
       
        # now gather from all processes 
        nstats = 4
        mat = np.zeros([np.shape(p)[1],nstats],dtype='f')
        mat[:,0] = np.squeeze(prod)
        mat[:,1] = np.squeeze(stats.reynolds_average(avg,p*s['11']))
        mat[:,2] = np.squeeze(stats.reynolds_average(avg,p*s['22']))
        mat[:,3] = np.squeeze(stats.reynolds_average(avg,p*s['12']))
        recvbuf = {}
        root=0
        for i in range(nstats):
            recvbuf[i]=comm.gather(mat[:,i], root=0)
            comm.Barrier()
            recvbuf[i] = np.array(recvbuf[i])

        if rank==0:
            total_array = np.zeros([Ny,nstats],dtype='f')
            for j in range(nstats):
                vec_array = np.reshape(recvbuf[j],[nblk_x,nblk_y,nblk_z,ny],order='F')
            

                # Now concat one column
                vec = vec_array[0,0,0,:];
                if (nblk_y>1):
                    for jj in range(1,nblk_y):
                        mat = np.squeeze(vec_array[0,jj,0,:])
                        vec = np.hstack([vec,mat])
                total_array[:,j] = vec
           
            outputfile = filename_prefix + 'pstrain_%04d.dat'%tid
            np.savetxt(outputfile,total_array,delimiter=' ')
            print("{}".format(outputfile))
