#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
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

debug = False
def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_start(default=0)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    start_index = 0
    if len(sys.argv) > 2:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 
    
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

    # Set up the derivative object
    Nx,Ny,Nz = reader.domain_size
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    der = cd.CompactDerivative(reader.grid_partition, 
            (dx, dy, dz), (10, 10, 10), periodic_dimensions)

    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    if rank==0: verbose=True
    else: verbose=False
    inp = nml.inputs(dirname,verbose)
    du = inp.du
    if rank==0: print("du = {}".format(inp.du))

    # Compute stats at each step:

    recvbuf = None
    for step in tID_list: 
        reader.step = step
        
        # density and streamwise vel, means
        rho, u, vpp, wpp = reader.readData( ('rho', 'u', 'v','w') )
        rbar = stats.reynolds_average(avg,rho)
        utilde = stats.favre_average(avg,rho,u,rho_bar=rbar)
        upp = u - utilde
       
        Rij = np.zeros([np.shape(u)[1],6],dtype='f')
        Rij[:,0] = np.squeeze(stats.favre_average(avg,rho,upp*upp,rho_bar=rbar))
        Rij[:,1] = np.squeeze(stats.favre_average(avg,rho,upp*vpp,rho_bar=rbar))
        Rij[:,2] = np.squeeze(stats.favre_average(avg,rho,upp*wpp,rho_bar=rbar))
        Rij[:,3] = np.squeeze(stats.favre_average(avg,rho,vpp*vpp,rho_bar=rbar))
        Rij[:,4] = np.squeeze(stats.favre_average(avg,rho,vpp*wpp,rho_bar=rbar))
        Rij[:,5] = np.squeeze(stats.favre_average(avg,rho,wpp*wpp,rho_bar=rbar))

        #recvbuf = comm.gather(Rij, root=0)
        if rank==0:
            #recvbuf = np.empty([Ny,6], dtype='f')
            #comm.Gather(np.float32(Rij), recvbuf, root=0)
            #print(np.shape(recvbuf))
            #recvbuf = np.array(recvbuf)
            #print(recvbuf[0,:,1])
            #print(recvbuf[1,:,1])
            recvbuf = Rij
            outputfile = filename_prefix + 'Rij_%04d.dat'%step
            np.savetxt(outputfile,recvbuf,delimiter=' ')
            print("{}".format(outputfile))
