import numpy
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=FutureWarning)
    import h5py

def saveFieldsToHdf5(self, fDict, path='./', filePrefix='hhdecomp_', time=None, count=None, numZeroFill=4):
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
        dset = h5_file.create_dataset(fname, domain_size, dtype=numpy.float64)
        dset[ixl:ixh, iyl:iyh, izl:izh] = field

    # Domain size
    dset = h5_file.create_dataset('NX', dtype=numpy.float64, data=domain_size[0])
    dset = h5_file.create_dataset('NY', dtype=numpy.float64, data=domain_size[1])
    dset = h5_file.create_dataset('NZ', dtype=numpy.float64, data=domain_size[2])

    # Grid coordinates
    dset = h5_file.create_dataset('x', (domain_size[0],), dtype=numpy.float64)
    dset = numpy.linspace(self.numSet.XMIN, self.numSet.XMAX-self.numSet.dx, num=self.numSet.NX)
    dset = h5_file.create_dataset('y', (domain_size[1],), dtype=numpy.float64)
    dset = numpy.linspace(self.numSet.YMIN, self.numSet.YMAX-self.numSet.dy, num=self.numSet.NY)
    dset = h5_file.create_dataset('z', (domain_size[2],), dtype=numpy.float64)
    dset = numpy.linspace(self.numSet.ZMIN, self.numSet.ZMAX, num=self.numSet.NZ)

    # Other information
    dset = h5_file.create_dataset('time', data=time, dtype=numpy.float64)
    #dset = time
    dset = h5_file.create_dataset('count', data=count, dtype=numpy.uint32)
    #dset = count

    h5_file.close()
    self.numSet.cout("Exported data to file: '{}'.".format(filename))
