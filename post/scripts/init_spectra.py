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
        print "  python {} <prefix> [tID_list]".format(sys.argv[0])
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
    if rank==0: print("\tdu = {}".format(inp.du))

    # Preallocate for means and derivatives
    Nsteps = np.size(tid_list)
    Nx,Ny,Nz = reader.domain_size
    
    # Compute stats at each step:
    i = 0
    for tid in tid_list:
        reader.step = tid
       
        # fluctations 
        u, v, w, r, p = reader.readData( ('u', 'v', 'w','rho','p') )
        rbar = stats.reynolds_average(avg,r)
        utilde = stats.favre_average(avg,r,u,rbar)
        vtilde = stats.favre_average(avg,r,v,rbar)
        wtilde = stats.favre_average(avg,r,w,rbar)
        ptilde = stats.favre_average(avg,r,p,rbar)
        rpp = r-rbar
        upp = u-utilde
        vpp = v-vtilde
        wpp = w-wtilde
        ppp = p-ptilde
   
        # Spectra in x at centerline, avged in z direction
        # take only first half and 3/2 rule, exclude the mean
        idx = range(1,3/2*Nx/2)
        rhat_3D_x = 1./Nx * np.fft.fft(rpp[:,Ny/2,:],axis=0)
        uhat_3D_x = 1./Nx * np.fft.fft(upp[:,Ny/2,:],axis=0)
        vhat_3D_x = 1./Nx * np.fft.fft(vpp[:,Ny/2,:],axis=0)
        what_3D_x = 1./Nx * np.fft.fft(wpp[:,Ny/2,:],axis=0)
        phat_3D_x = 1./Nx * np.fft.fft(ppp[:,Ny/2,:],axis=0)
        rhat_x = np.mean(rhat_3D_x[idx,:],axis=1)
        uhat_x = np.mean(uhat_3D_x[idx,:],axis=1) 
        vhat_x = np.mean(vhat_3D_x[idx,:],axis=1)
        what_x = np.mean(what_3D_x[idx,:],axis=1)
        phat_x = np.mean(phat_3D_x[idx,:],axis=1)

        # Write x spec to file 
        if rank==0:
            print("Time: {}".format(reader.time))
            array2save = np.zeros([idx[-1],5])
            array2save[:,0] = np.real(uhat_x);
            array2save[:,1] = np.real(vhat_x);
            array2save[:,2] = np.real(what_x);
            array2save[:,3] = np.real(rhat_x);
            array2save[:,4] = np.real(phat_x);
            outputfile = filename_prefix + "%04d_kx_real.dat"%tid
            np.savetxt(outputfile,array2save,delimiter=' ')
            array2save[:,0] = np.imag(uhat_x);
            array2save[:,1] = np.imag(vhat_x);
            array2save[:,2] = np.imag(what_x);
            array2save[:,3] = np.imag(rhat_x);
            array2save[:,4] = np.imag(phat_x);
            outputfile = filename_prefix + "%04d_kx_imag.dat"%tid
            np.savetxt(outputfile,array2save,delimiter=' ')
            print("\tDone writing to {}".format(outputfile))
        
        # Spectra in z at centerline, avged in x direction
        # take only first half and 3/2 rule, exclude the mean
        idx = range(1,3/2*Nz/2)
        rhat_3D_z = 1./Nz * np.fft.fft(rpp[:,Ny/2,:],axis=1)
        uhat_3D_z = 1./Nz * np.fft.fft(upp[:,Ny/2,:],axis=1)
        vhat_3D_z = 1./Nz * np.fft.fft(vpp[:,Ny/2,:],axis=1)
        what_3D_z = 1./Nz * np.fft.fft(wpp[:,Ny/2,:],axis=1)
        phat_3D_z = 1./Nz * np.fft.fft(ppp[:,Ny/2,:],axis=1)
        rhat_z = np.mean(rhat_3D_z[:,idx],axis=0)
        uhat_z = np.mean(uhat_3D_z[:,idx],axis=0) 
        vhat_z = np.mean(vhat_3D_z[:,idx],axis=0)
        what_z = np.mean(what_3D_z[:,idx],axis=0)
        phat_z = np.mean(phat_3D_z[:,idx],axis=0)

        # Save z spectrum
        if rank==0:
            print("Time: {}".format(reader.time))
            array2save = np.zeros([idx[-1],5])
            array2save[:,0] = np.real(uhat_z);
            array2save[:,1] = np.real(vhat_z);
            array2save[:,2] = np.real(what_z);
            array2save[:,3] = np.real(rhat_z);
            array2save[:,4] = np.real(phat_z);
            outputfile = filename_prefix + "%04d_kz_real.dat"%tid
            np.savetxt(outputfile,array2save,delimiter=' ')
            array2save[:,0] = np.imag(uhat_z);
            array2save[:,1] = np.imag(vhat_z);
            array2save[:,2] = np.imag(what_z);
            array2save[:,3] = np.imag(rhat_z);
            array2save[:,4] = np.imag(phat_z);
            outputfile = filename_prefix + "%04d_kz_imag.dat"%tid
            np.savetxt(outputfile,array2save,delimiter=' ')
            print("\tDone writing to {}".format(outputfile))
