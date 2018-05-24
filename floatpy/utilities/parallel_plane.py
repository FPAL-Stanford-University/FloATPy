from mpi4py import MPI
import numpy

from floatpy.parallel import t3dmod

class ParallelPlane(object):
    """
    Class to aggregate a plane from a 3D parallel decomposed field
    """

    def __init__(self, grid_partition, direction, index):
        """
        Constructor of the class.
        
        grid_partition : t3d object or the grid_partition property of the parallel data reader class
        direction : direction of plane normal (0 => x, 1 => y, 2 => z)
        index : plane normal index of points on the plane
        """
        
        if not isinstance(grid_partition, t3dmod.t3d):
            raise TypeError("The given grid partition object is not an instance of the t3d class!")
        self._grid_partition = grid_partition
        
        if direction < 0 or direction > 2:
            raise ValueError('Direction < 0 or > 2 is invalid!')
        self._direction = direction
        
        # Get size of chunk from all-direction domain decomposition of this processor and lo and hi of the chunk.
        self._3d_size = numpy.empty(3, dtype=numpy.int32)
        self._3d_lo   = numpy.empty(3, dtype=numpy.int32)
        self._3d_hi   = numpy.empty(3, dtype=numpy.int32)
        
        self._grid_partition.get_sz3d(self._3d_size)
        self._grid_partition.get_st3d(self._3d_lo)
        self._grid_partition.get_en3d(self._3d_hi)
        
        x_size = numpy.empty(3, dtype=numpy.int32)
        y_size = numpy.empty(3, dtype=numpy.int32)
        z_size = numpy.empty(3, dtype=numpy.int32)

        self._grid_partition.get_szx(x_size)
        self._grid_partition.get_szy(y_size)
        self._grid_partition.get_szz(z_size)

        # Convert to 0 based indexing.
        self._3d_lo = self._3d_lo - 1
        self._3d_hi = self._3d_hi - 1

        self._grid_size = numpy.array([ x_size[0], y_size[1], z_size[2] ])
        if (index < 0) or (index >= self._grid_size[direction]):
            raise ValueError('index has to be within bounds!')
        self._index = index

        self._is_involved = False
        if (self._index >= self._3d_lo[direction]) and (self._index <= self._3d_hi[direction]):
            self._is_involved = True
        
        if direction == 0:
            self._comm = MPI.Comm.f2py(self._grid_partition.commyz())
            self._plane_size = (self._grid_size[1], self._grid_size[2])
            self._lo = (self._3d_lo[1], self._3d_lo[2])
            self._hi = (self._3d_hi[1], self._3d_hi[2])
        
        elif direction == 1:
            self._comm = MPI.Comm.f2py(self._grid_partition.commxz())
            self._plane_size = (self._grid_size[0], self._grid_size[2])
            self._lo = (self._3d_lo[0], self._3d_lo[2])
            self._hi = (self._3d_hi[0], self._3d_hi[2])
        
        else:
            self._comm = MPI.Comm.f2py(self._grid_partition.commxy())
            self._plane_size = (self._grid_size[0], self._grid_size[1])
            self._lo = (self._3d_lo[0], self._3d_lo[1])
            self._hi = (self._3d_hi[0], self._3d_hi[1])

        self._rank  = self._comm.Get_rank()
        self._procs = self._comm.Get_size()

        self._all_lo = self._comm.gather(self._lo, root=0)
        self._all_hi = self._comm.gather(self._hi, root=0)


    def get_plane(self, array):
        """
        Method that returns the aggregated plane given a 3D decomposed array
        
        array : 3D decomposed numpy array
        """

        plane = None
        has_plane = False

        if self._is_involved:

            if self._direction == 0:
                my_plane = array[self._index-self._3d_lo[0],:,:]
            
            elif self._direction == 1:
                my_plane = array[:,self._index-self._3d_lo[1],:]
            
            else:
                my_plane = array[:,:,self._index-self._3d_lo[2]]

            data = self._comm.gather(my_plane, root=0)

            if self._rank == 0:
                plane = numpy.empty(self._plane_size, dtype=numpy.float64, order='F')

                for i in range(self._procs):
                    lo, hi = self._all_lo[i], self._all_hi[i]
                    plane[lo[0]:hi[0]+1, lo[1]:hi[1]+1] = data[i]

                has_plane = True


        return has_plane, plane

