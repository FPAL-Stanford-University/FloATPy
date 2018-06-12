from mpi4py import MPI
import numpy

from floatpy.parallel import t3dmod

class Reduction(object):
    """
    Class to sum or average a 3D parallel decomposed field
    """

    def __init__(self, grid_partition, directions):
        """
        Constructor of the class.
        
        grid_partition : t3d object or the grid_partition property of the parallel data reader class
        directions : tuple of logicals of size 3 defining the summation/averaging directions
        """
        
        if not isinstance(grid_partition, t3dmod.t3d):
            raise TypeError("The given grid partition object is not an instance of the t3d class!")
        self._grid_partition = grid_partition
        
        if len(directions) != 3:
            raise ValueError('Directions is invalid!')
        self._directions = directions
        
        # Get size of chunk from all-direction domain decomposition of this processor and lo and hi of the chunk.
        self._3d_size = numpy.empty(3, dtype=numpy.int32)
        self._3d_lo   = numpy.empty(3, dtype=numpy.int32)
        self._3d_hi   = numpy.empty(3, dtype=numpy.int32)
        
        self._grid_partition.get_sz3d(self._3d_size)
        self._grid_partition.get_st3d(self._3d_lo)
        self._grid_partition.get_en3d(self._3d_hi)
        
        # Convert to 0 based indexing.
        self._3d_lo = self._3d_lo - 1
        self._3d_hi = self._3d_hi - 1

        x_size = numpy.empty(3, dtype=numpy.int32)
        y_size = numpy.empty(3, dtype=numpy.int32)
        z_size = numpy.empty(3, dtype=numpy.int32)

        self._grid_partition.get_szx(x_size)
        self._grid_partition.get_szy(y_size)
        self._grid_partition.get_szz(z_size)

        self._grid_size = numpy.array([ x_size[0], y_size[1], z_size[2] ])

        self._comm  = MPI.Comm.f2py(self._grid_partition.comm3d())
        self._commx = MPI.Comm.f2py(self._grid_partition.commx())
        self._commy = MPI.Comm.f2py(self._grid_partition.commy())
        self._commz = MPI.Comm.f2py(self._grid_partition.commz())

        self._avg_lo = numpy.copy(self._3d_lo)
        self._avg_hi = numpy.copy(self._3d_hi)

        for i in range(3):
            if ( self._directions[i] ):
                self._avg_lo[i] = 0
                self._avg_hi[i] = 0

        if ( self._directions[0] ):
            if (not (self._directions[1] or self._directions[2])):
                self._comm_gather = MPI.Comm.f2py(self._grid_partition.commyz())

            else:
                if ( self._directions[1] and self._directions[2] ):
                    rank = self._comm.Get_rank()
                    self._comm_gather = self._comm.Split(color=rank, key=rank)

                else:
                    if ( self._directions[1] ):
                        self._comm_gather = MPI.Comm.f2py(self._grid_partition.commz())

                    elif ( self._directions[2] ):
                        self._comm_gather = MPI.Comm.f2py(self._grid_partition.commy())

        else:
            if (not (self._directions[1] or self._directions[2])):
                self._comm_gather = MPI.Comm.f2py(self._grid_partition.comm3d())

            else:
                if ( self._directions[1] and self._directions[2] ):
                    self._comm_gather = MPI.Comm.f2py(self._grid_partition.commx())
                else:

                    if ( self._directions[1] ):
                        self._comm_gather = MPI.Comm.f2py(self._grid_partition.commxz())

                    elif ( self._directions[2] ):
                        self._comm_gather = MPI.Comm.f2py(self._grid_partition.commxy())

        self._gather_procs = self._comm_gather.Get_size()
        
        self._all_avg_lo = self._comm_gather.allgather(self._avg_lo)
        self._all_avg_hi = self._comm_gather.allgather(self._avg_hi)


    def sum_x(self, array):
        """
        Method that returns the sum of a 3D decomposed array along X
        
        array : 3D decomposed numpy array
        """

        shape = list(array.shape)
        shape[0] = 1

        mysum = numpy.empty(shape, dtype=array.dtype, order='F')
        array.sum(axis=0, out=mysum, keepdims=True)

        sumx = numpy.empty(shape, dtype=array.dtype, order='F')
        self._commx.Allreduce(mysum, sumx, op=MPI.SUM)

        return sumx


    def average_x(self, array):
        """
        Method that returns the average of a 3D decomposed array along X
        
        array : 3D decomposed numpy array
        """

        return self.sum_x(array) / self._grid_size[0]


    def sum_y(self, array):
        """
        Method that returns the sum of a 3D decomposed array along Y
        
        array : 3D decomposed numpy array
        """

        shape = list(array.shape)
        shape[1] = 1

        mysum = numpy.empty(shape, dtype=array.dtype, order='F')
        array.sum(axis=1, out=mysum, keepdims=True)

        sumy = numpy.empty(mysum.shape, dtype=array.dtype, order='F')
        self._commy.Allreduce(mysum, sumy, op=MPI.SUM)

        return sumy


    def average_y(self, array):
        """
        Method that returns the average of a 3D decomposed array along Y
        
        array : 3D decomposed numpy array
        """

        return self.sum_y(array) / self._grid_size[1]


    def sum_z(self, array):
        """
        Method that returns the sum of a 3D decomposed array along Z
        
        array : 3D decomposed numpy array
        """

        shape = list(array.shape)
        shape[2] = 1

        mysum = numpy.empty(shape, dtype=array.dtype, order='F')
        array.sum(axis=2, out=mysum, keepdims=True)

        sumz = numpy.empty(mysum.shape, dtype=array.dtype, order='F')
        self._commz.Allreduce(mysum, sumz, op=MPI.SUM)

        return sumz


    def average_z(self, array):
        """
        Method that returns the average of a 3D decomposed array along Z
        
        array : 3D decomposed numpy array
        """

        return self.sum_z(array) / self._grid_size[2]


    def sum(self, array):
        """
        Method that returns the sum of a 3D decomposed array along directions
        that the object was initialized with
        
        array : 3D decomposed numpy array
        """

        summation = array

        if self._directions[0]:
            summation = self.sum_x(summation)

        if self._directions[1]:
            summation = self.sum_y(summation)

        if self._directions[2]:
            summation = self.sum_z(summation)

        return summation


    def average(self, array):
        """
        Method that returns the average of a 3D decomposed array along directions
        that the object was initialized with
        
        array : 3D decomposed numpy array
        """

        avg = array

        if self._directions[0]:
            avg = self.average_x(avg)

        if self._directions[1]:
            avg = self.average_y(avg)

        if self._directions[2]:
            avg = self.average_z(avg)

        return avg

    def allgather(self, array):
        """
        Method that returns the full domain data for a reduced field generated by
        this object.

        array : reduced 2D decomposed numpy array
        """

        shape = numpy.copy(self._grid_size)

        for i in range(3):
            if self._directions[i]:
                shape[i] = 1

        array_full = numpy.zeros(shape, dtype=numpy.float64, order='F')

        data = self._comm_gather.allgather(array)

        for i in range(self._gather_procs):
            lo, hi = self._all_avg_lo[i], self._all_avg_hi[i]
            array_full[lo[0]:hi[0]+1, lo[1]:hi[1]+1, lo[2]:hi[2]+1] = data[i]

        return array_full

