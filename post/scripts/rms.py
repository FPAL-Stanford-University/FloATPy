from mpi4py import MPI
import numpy as np
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
from mom_decomp import gather_y_plane

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:" 
        print "Computes pressure strain fields from filtered pressure and velocities" 
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

    # Set up the derivative object
    x, y, z = reader.readCoordinates()
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]

    def get_rms(q):
        q2 = q**2
        qmean = np.mean(np.mean(q2,axis=-1),axis=0)
        return np.sqrt(qmean)

    # Compute stats at each step:
    for tID in tID_list:
        if rank==0: print('Reading data...')
        reader.step = tID
        r,u,v,w,p = reader.readData( ('rho','u','v','w','p') )
        rbar = stats.reynolds_average(avg,r)
        pbar = stats.reynolds_average(avg,p)
        utilde = stats.favre_average(avg,r,u,rho_bar=rbar)
        vtilde = stats.favre_average(avg,r,v,rho_bar=rbar)
        wtilde = stats.favre_average(avg,r,w,rho_bar=rbar)
        upp = np.array(u-utilde)
        vpp = np.array(v-vtilde)
        wpp = np.array(w-wtilde)
        pp = np.array(p-pbar)
        rp = np.array(r-rbar)


        
        planes={}
        planes['upp'] = gather_y_plane(comm,reader,upp,y,yslice=0)
        planes['vpp'] = gather_y_plane(comm,reader,vpp,y,yslice=0)
        planes['wpp'] = gather_y_plane(comm,reader,wpp,y,yslice=0)
        planes['pp']   = gather_y_plane(comm,reader,pp,y,yslice=0)
        planes['rp']   = gather_y_plane(comm,reader,rp,y,yslice=0)

        if rank==0: 
            rms = {}
            for key in planes.keys():
                rms[key] = get_rms(planes[key])
            planes = rms

            print('Writing to file...')
            dir_out = dirname.split('/projects/ShockInducedMix/')[-1]
            dir_out = '/home/kmatsuno/' + dir_out + '/'
            outputfile = dir_out+ '/centerline_rms_%04d.h5'%tID
            hf = h5py.File(outputfile, 'w')
            for key in planes.keys():
                hf.create_dataset(key, data=planes[key])
            hf.close()
            print('Saved to %s'%outputfile)
