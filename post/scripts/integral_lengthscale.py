from mpi4py import MPI
import numpy as np
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import custom_functions as fn

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz

def autocorr(uprime2D):
    uhat = np.fft.fftshift(np.fft.fft2(uprime2D))
    R_fft = np.square(np.abs(uhat))
    R = np.abs(np.fft.fftshift(np.fft.ifft2(R_fft)))
    return R/R.max()
    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_start(default=0)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    start_index = 0
    if len(sys.argv) > 2:
        start_index = int(sys.argv[2])
    outputfile  = filename_prefix + "integral_lengthscale.dat"
    
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
    x, y, z = reader.readCoordinates()
    dx,dy,dz = grid_res(x,y,z)
    [Nx,Ny,Nz] = reader.domain_size
    [Lx,Ly,Lz] = [x[-1,0,0],y[0,-1,0]-y[0,0,0],z[0,0,-1]]

    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    inp = nml.inputs(dirname,verbose=True)
    du = inp.du
    if rank==0: print("du = {}".format(inp.du))

    # Which slice?
    yi = Ny/2

    # Compute stats at each step:
    Nsteps = np.size(steps[start_index:])-1
    qlist  = ['u','v','p']
    Nq     = np.size(qlist)
    time   = np.empty(Nsteps)
    lx = np.empty([Nsteps,Nq])
    lz = np.empty([Nsteps,Nq])
    i = 0
    if rank == 0:print("Time \t lx \t lz")
    for step in steps[start_index:-1]:
        reader.step = step
        time[i] = reader.time
        
        # get mean density 
        r = reader.readData( 'rho' )
        rho = np.squeeze(r)
        rho_bar= stats.reynolds_average(avg,rho)

        nq = 0
        for qname in qlist:
            # Favre fluctuations q"
            q = reader.readData(qname)
            q = np.squeeze(q)
            qtilde = stats.favre_average(avg,rho,q,rho_bar)
            qpp2D = q[:,yi,:]-qtilde[:,yi,:]
         
            # Autocorrelations
            R = fn.autocorr(qpp2D)

            # Integral lengthscales in x and z: 
            lx[i,nq] = 1./dx * sum(R[Nx/2:, Nz/2]) 
            lz[i,nq] = 1./dz * sum(R[Nx/2, Nz/2:]) 
            nq += 1

        if rank == 0:
            print("{} \t {} \t {}".format(reader.time,lx[i],lz[i]))
            i = i+1; 
    
    # Write to file 
    if rank==0:
        array2save = np.empty((Nsteps,1+2*Nq+1))
        array2save[:,0] = time
        array2save[:,1:Nq+1] = lx/Lx
        array2save[:,Nq+2:] = lz/Lz
        np.savetxt(outputfile,array2save,delimiter=' ')
        print("Done writing to {}".format(outputfile))
