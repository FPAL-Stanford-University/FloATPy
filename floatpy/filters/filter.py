from mpi4py import MPI
import numpy
import pycf90
import pygaussian

from floatpy.parallel import t3dmod
from floatpy.utilities import data_reshaper

class Filter(object):
    """
    Class to perform parallel filter operations
    """

    def __init__(self, grid_partition, filter_type, dimension=3, periodic_dimensions=(False, False, False)):
        """
        Constructor of the class.

        grid_partition : grid_partition property (t3d object) of the parallel data reader class
        filter_type : string iterable of size 3 with the type of filter to use in each direction
                      options are "compact" and "gaussian"
        periodic_dimensions : boolean iterable of size 3
        """

        if not isinstance(grid_partition, t3dmod.t3d):
            raise RuntimeError("The given grid partition object is not an instance of the t3d class!")
        
        if dimension < 1 or dimension > 3:
            raise RuntimeError("Class only works with data with number of dimensions between 1 and 3!")
        
        self._dim = dimension
        
        if len(filter_type) < self._dim:
            raise RuntimeError("Size of 'filter_type' is smaller than problem dimension!")

        if len(periodic_dimensions) < self._dim:
            raise RuntimeError("Size of 'periodic_dimensions' is smaller than problem dimension!")

        for i in range(self._dim):
            if filter_type[i] not in ['compact', 'gaussian']:
                raise RuntimeError("filter_type[%d] has to be one of {'compact', 'gaussian'}" %i)

        self._filter_type = tuple(filter_type)
        self._grid_partition = grid_partition

        self._periodic = tuple(periodic_dimensions)

        self._chunk_3d_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_3d_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_3d_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._grid_partition.get_sz3d(self._chunk_3d_size)
        self._grid_partition.get_st3d(self._chunk_3d_lo)
        self._grid_partition.get_en3d(self._chunk_3d_hi)
        self._chunk_3d_lo = self._chunk_3d_lo - 1 # Convert to 0 based indexing
        self._chunk_3d_hi = self._chunk_3d_hi - 1 # Convert to 0 based indexing
        
        self._chunk_x_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_x_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_x_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._grid_partition.get_szx(self._chunk_x_size)
        self._grid_partition.get_stx(self._chunk_x_lo)
        self._grid_partition.get_enx(self._chunk_x_hi)
        self._chunk_x_lo = self._chunk_x_lo - 1 # Convert to 0 based indexing
        self._chunk_x_hi = self._chunk_x_hi - 1 # Convert to 0 based indexing
        
        self._chunk_y_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_y_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_y_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._grid_partition.get_szy(self._chunk_y_size)
        self._grid_partition.get_sty(self._chunk_y_lo)
        self._grid_partition.get_eny(self._chunk_y_hi)
        self._chunk_y_lo = self._chunk_y_lo - 1 # Convert to 0 based indexing
        self._chunk_y_hi = self._chunk_y_hi - 1 # Convert to 0 based indexing
        
        self._chunk_z_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_z_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_z_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._grid_partition.get_szz(self._chunk_z_size)
        self._grid_partition.get_stz(self._chunk_z_lo)
        self._grid_partition.get_enz(self._chunk_z_hi)
        self._chunk_z_lo = self._chunk_z_lo - 1 # Convert to 0 based indexing
        self._chunk_z_hi = self._chunk_z_hi - 1 # Convert to 0 based indexing

        if self._dim == 1:
            if self._chunk_3d_size[1] != 1 or \
               self._chunk_3d_size[2] != 1:
                raise RuntimeError("Make sure grid_partition is consistent with 1D problem!")
        
        if self._dim == 2:
            if self._chunk_3d_size[2] != 1:
                raise RuntimeError("Make sure grid_partition is consistent with 2D problem!")
        
        self._nx = self._chunk_x_size[0]
        self._ny = self._chunk_y_size[1]
        self._nz = self._chunk_z_size[2]
        
        self._xfil = None
        self._yfil = None
        self._zfil = None
        
        if self._filter_type[0] == 'compact':
            self._xfil = pycf90.cf90stuff.cf90( self._nx, self._periodic[0] )
        elif self._filter_type[0] == 'gaussian':
            self._xfil = pygaussian.gaussianstuff.gaussian( self._nx, self._periodic[0] )
        
        if self._dim > 1:
            if self._filter_type[1] == 'compact':
                self._yfil = pycf90.cf90stuff.cf90( self._ny, self._periodic[1] )
            elif self._filter_type[1] == 'gaussian':
                self._yfil = pygaussian.gaussianstuff.gaussian( self._ny, self._periodic[1] )
        
        if self._dim > 2:
            if self._filter_type[2] == 'compact':
                self._zfil = pycf90.cf90stuff.cf90( self._nz, self._periodic[2] )
            elif self._filter_type[2] == 'gaussian':
                self._zfil = pygaussian.gaussianstuff.gaussian( self._nz, self._periodic[2] )
        
        # Initialize the data reshaper.
        self._data_reshaper = data_reshaper.DataReshaper(self._dim, data_order='F')


    def filter_x(self, data, data_filtered=None, component_idx=None, bc=(0,0)):
        """
        Method to filter data in the first direction.

        data : input 3D numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        data_filtered : optional output 3D numpy array in Fortran contiguous layout. This array must be consistent with the 3D
                        decomposition and the problem dimension. This method will return data_filtered if data_filtered is None
        component_idx : index of component in data for filtering. None if there is only one component in the data
        bc : integer iterable of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        data_shape = data.shape

        if component_idx is not None:
            data_shape = data_shape[0:-1]

        if self._dim == 1:
            if len(data_shape) != 1:
                raise RuntimeError("Make sure data is 1D!")
            if data_shape[0] != self._chunk_3d_size[0]:
                raise RuntimeError("Make sure data is of the same size as in grid_partition!")

        elif self._dim == 2:
            if len(data_shape) != 2:
                raise RuntimeError("Make sure data is 2D!")
            if data_shape[0] != self._chunk_3d_size[0] or \
               data_shape[1] != self._chunk_3d_size[1]:
                raise RuntimeError("Make sure data is of the same size as in grid_partition!")

        else:
            if len(data_shape) != 3:
                raise RuntimeError("Make sure data is 3D!")
            if data_shape[0] != self._chunk_3d_size[0] or \
               data_shape[1] != self._chunk_3d_size[1] or \
               data_shape[2] != self._chunk_3d_size[2]:
                raise RuntimeError("Make sure data is of the same size as in grid_partition!")

        return_data_filtered = True
        if data_filtered is None:
            data_filtered = numpy.empty(self._chunk_3d_size[0:self._dim], dtype=numpy.float64, order='F')
        else:
            return_data_filtered = False

        data_x          = numpy.empty( self._chunk_x_size, dtype=numpy.float64, order='F' )
        data_filtered_x = numpy.empty( self._chunk_x_size, dtype=numpy.float64, order='F' )

        data_3d = []
        if component_idx is None:
            data_3d = self._data_reshaper.reshapeTo3d(data)
        else:
            data_3d = self._data_reshaper.reshapeTo3d(data, component_idx)

        self._grid_partition.transpose_3d_to_x( data_3d, data_x )

        self._xfil.filter1(data_x, data_filtered_x, self._chunk_x_size[1], self._chunk_x_size[2], bc1_=bc[0], bcn_=bc[1])

        data_filtered_3d = self._data_reshaper.reshapeTo3d(data_filtered)
        self._grid_partition.transpose_x_to_3d( data_filtered_x, data_filtered_3d )
        data_filtered = self._data_reshaper.reshapeFrom3d(data_filtered_3d)

        if return_data_filtered:
            return data_filtered


    def filter_y(self, data, data_filtered=None, component_idx=None, bc=(0,0)):
        """
        Function to filter data in the second direction.

        data : input 3D numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        data_filtered : optional output 3D numpy array in Fortran contiguous layout. This array must be consistent with the 3D
                        decomposition and the problem dimension. This method will return data_filtered if data_filtered is None
        component_idx : index of component in data for filtering. None if there is only one component in the data
        bc : integer tuple of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        data_shape = data.shape

        if component_idx is not None:
            data_shape = data_shape[0:-1]

        if self._dim == 1:
            raise RuntimeError("There is no ddy for 1D problem!")

        elif self._dim == 2:
            if len(data_shape) != 2:
                raise RuntimeError("Make sure data is 2D!")
            if data_shape[0] != self._chunk_3d_size[0] or \
               data_shape[1] != self._chunk_3d_size[1]:
                raise RuntimeError("Make sure data is of the same size as in grid_partition!")

        else:
            if len(data_shape) != 3:
                raise RuntimeError("Make sure data is 3D!")
            if data_shape[0] != self._chunk_3d_size[0] or \
               data_shape[1] != self._chunk_3d_size[1] or \
               data_shape[2] != self._chunk_3d_size[2]:
                raise RuntimeError("Make sure data is of the same size as in grid_partition!")

        return_data_filtered = True
        if data_filtered is None:
            data_filtered = numpy.empty(self._chunk_3d_size[0:self._dim], dtype=numpy.float64, order='F')
        else:
            return_data_filtered = False

        data_y          = numpy.empty( self._chunk_y_size, dtype=numpy.float64, order='F' )
        data_filtered_y = numpy.empty( self._chunk_y_size, dtype=numpy.float64, order='F' )

        data_3d = []
        if component_idx is None:
            data_3d = self._data_reshaper.reshapeTo3d(data)
        else:
            data_3d = self._data_reshaper.reshapeTo3d(data, component_idx)

        self._grid_partition.transpose_3d_to_y( data_3d, data_y )
        
        self._yfil.filter2(data_y, data_filtered_y, self._chunk_y_size[0], self._chunk_y_size[2], bc1_=bc[0], bcn_=bc[1])
        
        data_filtered_3d = self._data_reshaper.reshapeTo3d(data_filtered)
        self._grid_partition.transpose_y_to_3d( data_filtered_y, data_filtered_3d )
        data_filtered = self._data_reshaper.reshapeFrom3d(data_filtered_3d)

        if return_data_filtered:
            return data_filtered


    def filter_z(self, data, data_filtered=None, component_idx=None, bc=(0,0)):
        """
        Function to filter data in the third direction.

        data : input 3D numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        data_filtered : optional output 3D numpy array in Fortran contiguous layout. This array must be consistent with the 3D
                        decomposition and the problem dimension. This method will return data_filtered if data_filtered is None
        component_idx : index of component in data for filtering. None if there is only one component in the data
        bc : integer tuple of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        data_shape = data.shape

        if component_idx is not None:
            data_shape = data_shape[0:-1]

        if self._dim == 1:
            raise RuntimeError("There is no ddz for 1D problem!")

        elif self._dim == 2:
            raise RuntimeError("There is no ddz for 2D problem!")

        else:
            if len(data_shape) != 3:
                raise RuntimeError("Make sure data is 3D!")
            if data_shape[0] != self._chunk_3d_size[0] or \
               data_shape[1] != self._chunk_3d_size[1] or \
               data_shape[2] != self._chunk_3d_size[2]:
                raise RuntimeError("Make sure data is of the same size as in grid_partition!")

        return_data_filtered = True
        if data_filtered is None:
            data_filtered = numpy.empty(self._chunk_3d_size[0:self._dim], dtype=numpy.float64, order='F')
        else:
            return_data_filtered = False

        data_z          = numpy.empty( self._chunk_z_size, dtype=numpy.float64, order='F' )
        data_filtered_z = numpy.empty( self._chunk_z_size, dtype=numpy.float64, order='F' )

        data_3d = []
        if component_idx is None:
            data_3d = self._data_reshaper.reshapeTo3d(data)
        else:
            data_3d = self._data_reshaper.reshapeTo3d(data, component_idx)

        self._grid_partition.transpose_3d_to_z( data_3d, data_z )
        
        self._zfil.filter3(data_z, data_filtered_z, self._chunk_z_size[0], self._chunk_z_size[1], bc1_=bc[0], bcn_=bc[1])
        
        self._grid_partition.transpose_z_to_3d( data_filtered_z, data_filtered )

        if return_data_filtered:
            return data_filtered


    def filter_all(self, data, data_filtered=None, component_idx=None, x_bc=(0,0), y_bc=(0,0), z_bc=(0,0), ntimes=1):
        """
        Function to filter data in all directions.

        data : input 3D numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        data_filtered : optional output 3D numpy array in Fortran contiguous layout. This array must be consistent with the 3D
                        decomposition and the problem dimension. This method will return data_filtered if data_filtered is None
        component_idx : index of component in data for filtering. None if there is only one component in the data
        *_bc : integer tuple of size 2 with the boundary condition at the left and right.
               0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        data_shape = data.shape

        if component_idx is not None:
            data_shape = data_shape[0:-1]

        if self._dim == 1:
            if len(data_shape) != 1:
                raise RuntimeError("Make sure data is 1D!")
            if data_shape[0] != self._chunk_3d_size[0]:
                raise RuntimeError("Make sure data is of the same size as in grid_partition!")

        elif self._dim == 2:
            if len(data_shape) != 2:
                raise RuntimeError("Make sure data is 2D!")
            if data_shape[0] != self._chunk_3d_size[0] or \
               data_shape[1] != self._chunk_3d_size[1]:
                raise RuntimeError("Make sure data is of the same size as in grid_partition!")

        else:
            if len(data_shape) != 3:
                raise RuntimeError("Make sure data is 3D!")
            if data_shape[0] != self._chunk_3d_size[0] or \
               data_shape[1] != self._chunk_3d_size[1] or \
               data_shape[2] != self._chunk_3d_size[2]:
                raise RuntimeError("Make sure data is of the same size as in grid_partition!")

        return_data_filtered = True
        if data_filtered is None:
            data_filtered = numpy.empty(self._chunk_3d_size[0:self._dim], dtype=numpy.float64, order='F')
        else:
            return_data_filtered = False

        tmp           = numpy.empty( data_shape, dtype=numpy.float64, order='F' )
        data_filtered = numpy.empty( data_shape, dtype=numpy.float64, order='F' )

        self.filter_x(data, data_filtered, bc=x_bc)
        for i in range(ntimes-1):
            tmp = numpy.copy(data_filtered, order='F')
            self.filter_x(tmp, data_filtered, bc=x_bc)

        if self._dim > 1:
            self.filter_y(data_filtered, tmp, bc=y_bc)
            for i in range(ntimes-1):
                data_filtered = numpy.copy(tmp, order='F')
                self.filter_y(data_filtered, tmp, bc=y_bc)
                data_filtered = numpy.copy(tmp, order='F')
            data_filtered = numpy.copy(tmp, order='F')

        if self._dim > 2:
            self.filter_z(tmp, data_filtered, bc=z_bc)
            for i in range(ntimes-1):
                tmp = numpy.copy(data_filtered, order='F')
                self.filter_z(tmp, data_filtered, bc=z_bc)

        if return_data_filtered:
            return data_filtered

