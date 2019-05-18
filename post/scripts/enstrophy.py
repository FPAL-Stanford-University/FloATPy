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
    if len(sys.argv) < 3:
        print "Usage:" 
        print "This script writes the full enstrophy spectrum for one time."
        print "  python {} <prefix> <tID>".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    tid_list = map(int, sys.argv[2].strip('[]').split(',')) 
    
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
    inp = nml.inputs(dirname,verbose=True)
    du = inp.du
    if rank==0: print("du = {}".format(inp.du))

    # Preallocate for means and derivatives
    Nsteps = np.size(tid_list)
    Nx,Ny,Nz = reader.domain_size
    time   = np.empty(Nsteps)
    eta    = np.empty(Nsteps)
    dtheta = np.empty(Nsteps)

    # Compute stats at each step:
    i = 0
    for tid in tid_list:
        reader.step = tid
        time[i] = reader.time
       
        # vorticity and enstrophy 
        u, v, w = reader.readData( ('u', 'v', 'w') )
        vorticity = der.curl(u, v, w, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        enstrophy = stats.get_enstrophy(vorticity, reader.grid_partition, 
                dx, dy, dz, volume_integrated=False)
        print("Max/min enstrophy: {}, {}".format(np.amax(enstrophy),np.amin(enstrophy)) );

        # Spectra in x,z at the centerline:
        spec3D_x = 1./Nx * np.fft.fft(enstrophy[:,Ny/2,:],axis=0)
        spec3D_z = 1./Nx * np.fft.fft(enstrophy[:,Ny/2,:],axis=1)
        spec_x = np.mean(abs(spec3D_x),axis=1) # avg in z for kx spectrum
        spec_z = np.mean(abs(spec3D_z),axis=0) # avg in x for kz spectrum

        # Write to file 
        if rank==0:
            # Save x spectrum
            array2save = np.empty([Nx,1])
            array2save[:,0] = spec_x;
            outputfile  = "%04d_enstrophy_spectrum_x.dat"%tid
            np.savetxt(filename_prefix+outputfile,array2save,delimiter=' ')
            # Save z spectrum
            array2save = np.empty([Nz,1])
            array2save[:,0] = spec_z;
            outputfile  = "%04d_enstrophy_spectrum_z.dat"%tid
            np.savetxt(filename_prefix+outputfile,array2save,delimiter=' ')
            
            print("Time: {}\tDone writing to {}".format(time[i], outputfile))
            i = i+1; 
