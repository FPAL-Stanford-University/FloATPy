#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
from mpi4py import MPI
import numpy as np
import os
import sys
#sys.path.insert(0,'/home/kmatsuno/h5py/')
sys.path.insert(0,'/home/kmatsuno/h5py/build/lib.linux-x86_64-2.7/')
#print(sys.path)
import h5py
from shutil import copyfile

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red
import statistics as stats
import get_namelist as nml
from SettingLib import NumSetting

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
   
class h5_writer:
    def __init__(self,numSet):
        self.numSet = numSet
        self.count = 0
        self.time = 0
    
    def update(self,count,time):
        self.count = count
        self.time = time

    def saveFields(self, fDict, path='./', filePrefix='hhdecomp_', time=None, count=None, numZeroFill=4):
        """
        Save computation ressult to *.h5 file.
        """
        domain_size = (self.numSet.NX, self.numSet.NY, self.numSet.NZ)
        ixl, ixh = self.numSet.chunk_3d_lo[0], self.numSet.chunk_3d_hi[0]+1
        iyl, iyh = self.numSet.chunk_3d_lo[1], self.numSet.chunk_3d_hi[1]+1
        izl, izh = self.numSet.chunk_3d_lo[2], self.numSet.chunk_3d_hi[2]+1
        if time==None: time = self.time
        if count==None: count = self.count

        filename = path + filePrefix + str(self.count).zfill(numZeroFill) + ".h5"
        h5_file = h5py.File(filename, 'w', driver='mpio', comm=self.numSet.comm)

        # Data to be saved
        # Fields
        for fname in fDict:
            field = fDict[fname]
            dset = h5_file.create_dataset(fname, domain_size, dtype=np.float64)
            dset[ixl:ixh, iyl:iyh, izl:izh] = field

        # Domain size
        dset = h5_file.create_dataset('NX', dtype=np.float64, data=domain_size[0])
        dset = h5_file.create_dataset('NY', dtype=np.float64, data=domain_size[1])
        dset = h5_file.create_dataset('NZ', dtype=np.float64, data=domain_size[2])

        # Grid coordinates
        dset = h5_file.create_dataset('x', (domain_size[0],), dtype=np.float64)
        dset = np.linspace(self.numSet.XMIN, self.numSet.XMAX-self.numSet.dx, num=self.numSet.NX)
        dset = h5_file.create_dataset('y', (domain_size[1],), dtype=np.float64)
        dset = np.linspace(self.numSet.YMIN, self.numSet.YMAX-self.numSet.dy, num=self.numSet.NY)
        dset = h5_file.create_dataset('z', (domain_size[2],), dtype=np.float64)
        dset = np.linspace(self.numSet.ZMIN, self.numSet.ZMAX, num=self.numSet.NZ)

        # Other information
        dset = h5_file.create_dataset('time', data=time, dtype=np.float64)
        #dset = time
        dset = h5_file.create_dataset('count', data=count, dtype=np.uint32)
        #dset = count

        h5_file.close()
        self.numSet.cout("Exported data to file: '{}'.".format(filename))

if __name__ == '__main__':
    dir_in  = sys.argv[1]
    tID     = int(sys.argv[2])
    dir_out = dir_in 
    fname_in = dir_in + '/restart_%04d'%tID + '.h5' 
    fname_out = dir_out + '/shearlayer_tmp.h5' 
    
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
    serial_reader = por.PadeopsReader(dir_in+'/restart_',
            periodic_dimensions=periodic_dimensions)
    reader = pdr.ParallelDataReader(comm, serial_reader)
    avg = red.Reduction(reader.grid_partition, periodic_dimensions)

    # setup the inputs object
    inp = nml.inputs(dir_out,verbose=(rank==0))
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dir_in,verbose=(rank==0))
    
    # Get grid partition size
    if rank==0: print('Setting up fft...')
    settings = NumSetting( comm, reader.grid_partition, 
             NX=Nx, NY=Ny, NZ=Nz,
             XMIN=0,        XMAX=Lx,
             YMIN=-Ly/2.,   YMAX=Ly/2.,
             ZMIN=0,        ZMAX=Lz,
             order=10)
    szx, szy, szz = settings.chunk_3d_size
    
    # Set up parallel writer
    writer = h5_writer(settings)

    # Read input restart file
    reader.step = tID
    qlist = ('rhou', 'rhov', 'rhow', 'TE', 'rhoY_0001','rhoY_0002')
    #ru, rv, rw, TE, r1, r2  = reader.readData( qlist )
    q = reader.readData( qlist )
    qDict = {}
    for i in range(len(qlist)): 
        qDict[qlist[i]] = q[i]
    
    writer.saveFields(qDict) 

    comm.Barrier()
    sys.exit()
