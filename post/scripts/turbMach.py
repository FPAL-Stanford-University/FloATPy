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

    x, y, z = reader.readCoordinates()
    Nx,Ny,Nz = reader.domain_size
    dx,dy,dz = grid_res(x,y,z)
    
    # setup the inputs object, get grid info
    dirname = os.path.dirname(filename_prefix)
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=(rank==0))
    Ny = int(Ny)
    if rank==0: verbose=True
    else: verbose=False
    inp = nml.inputs(dirname,verbose)
    du = inp.du
    if rank==0: print("du = {}".format(inp.du))
    
    # Get the grid partition information
    nx,ny,nz = reader._full_chunk_size
    nblk_x = int(np.round(Nx/nx))
    nblk_y = int(np.round(Ny/ny))
    nblk_z = int(np.round(Nz/nz))
    if rank==0: print("Processor decomp: {}x{}x{}".format(nblk_x,nblk_y,nblk_z))

    # Compute stats at each step:
    for step in tID_list: 
        reader.step = step

        # TKE
        rho, u, v, w, T = reader.readData( ('rho', 'u', 'v', 'w','T') )
        tke = stats.TKE(rho, u, v, w, reader.grid_partition, 
                avg, dx, dy, dz, volume_integrated=False)
        
        # Speed of sound, xz average
        c3D = np.sqrt(inp.gam*T)
        cbar = stats.reynolds_average(avg, c3D)
        Mt = np.squeeze(tke**0.5/cbar)

        # now gather from all processes
        root=0
        comm.Barrier()
        recvbuf=comm.gather(Mt, root)
        comm.Barrier()
        recvbuf = np.array(recvbuf)

        if rank==0:
            try: np.shape(recvbuf)[1] #should be a 2d array
            except:
                print("ERROR: Shape mismatch, recvbuf {}".format(np.shape(recvbuf)))
                sys.exit()
            
            # Stack the vectors into the correct order
            vec_array = np.reshape(recvbuf,[nblk_x,nblk_y,nblk_z,ny],order='F')
            
            # Now concat one column
            vec = vec_array[0,0,0,:];
            if (nblk_y>1):
                for i in range(1,nblk_y):
                    mat = np.squeeze(vec_array[0,i,0,:])
                    vec = np.hstack([vec,mat])

            if (np.shape(vec)[0] < Ny):
                print("ERROR: Shape mismatch, vec {}".format(np.shape(vec)))
                sys.exit()
            else:
                outputfile = filename_prefix + 'Mt_%04d.dat'%step
                np.savetxt(outputfile,vec,delimiter=' ')
                print("{}".format(outputfile))
