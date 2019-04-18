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
        print "This script writes the mean momentum budget for one time."
        print "  python {} <prefix> <tID list>".format(sys.argv[0])
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

    # Some temporary variables to write to
    tmp1 = np.empty( der._chunk_3d_size,dtype=np.float64, order='F' )
    tmp2 = tmp1; tmp3 = tmp1


    # Compute stats at each step:
    for tid in tid_list:
        reader.step = tid
       
        # mean density
        q = reader.readData( 'rho' )
        r = np.squeeze(q)
        rbar = stats.reynolds_average(avg,r)

        # fluct
        u, v, w, p, mu = reader.readData( ('u','v','w','p','mu') )
        ubar = stats.reynolds_average(avg,u)
        vbar = stats.reynolds_average(avg,v)
        wbar = stats.reynolds_average(avg,w)
        pbar = stats.reynolds_average(avg,p)
        utilde = stats.favre_average(avg, r, u, rho_bar=rbar)
        vtilde = stats.favre_average(avg, r, v, rho_bar=rbar)
        wtilde = stats.favre_average(avg, r, w, rho_bar=rbar)


        # make the means 3D
        #rbar = np.tile(rbar,(Nx,1,Nz));
        #ubar = np.tile(ubar,(Nx,1,Nz));
        #pbar = np.tile(pbar,(Nx,1,Nz));
        #utilde = np.tile(utilde,(Nx,1,Nz));
        #vtilde = np.tile(vtilde,(Nx,1,Nz));
        #wtilde = np.tile(wtilde,(Nx,1,Nz));
        rbar3D = np.zeros([Nx,Ny,Nz])
        ubar3D = np.zeros([Nx,Ny,Nz])
        pbar3D = np.zeros([Nx,Ny,Nz])
        utilde3D = np.zeros([Nx,Ny,Nz])
        vtilde3D = np.zeros([Nx,Ny,Nz])
        wtilde3D = np.zeros([Nx,Ny,Nz])
        for i in range(Nx-1):
            for j in range(Nz-1):
                rbar3D[i,:,j]   = np.squeeze(rbar)
                ubar3D[i,:,j]   = np.squeeze(ubar)
                pbar3D[i,:,j]   = np.squeeze(pbar)
                utilde3D[i,:,j] = np.squeeze(utilde)
                vtilde3D[i,:,j] = np.squeeze(vtilde)
                wtilde3D[i,:,j] = np.squeeze(wtilde)
        
        # Fluctuations
        upp = u-utilde
        vpp = v-vtilde
        wpp = w-wtilde

        ###### MEAN X MOMENTUM #######

        # Term I with gradients of favre averages
        der.ddx(rbar3D*utilde3D*utilde3D,tmp1,bc=x_bc)
        der.ddy(rbar3D*utilde3D*vtilde3D,tmp2,bc=y_bc)
        der.ddz(rbar3D*utilde3D*wtilde3D,tmp3,bc=z_bc)
        I = stats.reynolds_average(avg,-(tmp1+tmp2+tmp3))

        # Term II with gradients of favre flucts
        der.ddx(rbar3D*upp*upp,tmp1,bc=x_bc)
        der.ddy(rbar3D*upp*vpp,tmp2,bc=y_bc)
        der.ddz(rbar3D*upp*wpp,tmp3,bc=z_bc)
        II = stats.reynolds_average(avg,-(tmp1+tmp2+tmp3))

        # Term III Pressure gradient
        der.ddx(pbar3D,tmp1,bc=x_bc)
        III = stats.reynolds_average(avg,-tmp1)

        # Term IV Shear stresses
        ux,uy,uz = der.gradient(ubar3D,x_bc,y_bc,z_bc)
        uyx,uyy,uzy = der.gradient(mu*uy,x_bc,y_bc,z_bc)
        IV = stats.reynolds_average(avg,uyy)
        dudx,dudy,dudz = der.gradient(u, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        uy = stats.reynolds_average(avg,dudy)

        # ddt term: one sided fd 3fi-4fi+1+fi+2 / 2h
        # reader.step = tid-1


        # Write to file 
        if rank==0:
            time = reader.time
            outputfile = filename_prefix + "%04d_budget_mean_momentum.npz"%tid
            np.savez(outputfile,time=time,I=I,II=II,III=III,IV=IV,
                    ubar3D=ubar3D,uy=uy,mu=mu)
            print("Time: {}\tDone writing to {}".format(reader.time, outputfile))
