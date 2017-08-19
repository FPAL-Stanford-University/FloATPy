import numpy
import compact.pycd06 as pycd06
import compact.pycd10 as pycd10

from floatpy.parallel import _t3dmod
from floatpy.utilities import data_reshaper

class CompactDifferentiator(object):
    """
    Class to perform parallel compact finite difference routines.
    """
    
    def __init__(self, grid_partition, grid_spacing, order, dimension=3, periodic_dimensions=(False, False, False)):
        """
        Constructor of the class.
        
        grid_partition : the grid_partition property (t3d object) of the parallel data reader class
        grid_spacing : iterable of floats for the grid spacing in each direction
        order : integer iterable with each value in {6, 10} representing the order of accuracy of the derivatives
                in each direction
        dimension : dimension of problem
        periodic_dimensions : iterable of boolean descibing whether the periodicity in each direction
        """
        
        if not isinstance(grid_partition, _t3dmod.t3d):
            raise RuntimeError("The given grid partition object is not an instance of the t3d class!")
        
        if dimension < 1 or dimension > 3:
            raise RuntimeError("Class only works with data with number of dimensions between 1 and 3!")
        
        self._dim = dimension
        
        if len(grid_spacing) < self._dim:
            raise RuntimeError("Size of 'grid_spacing' is smaller than problem dimension!")
        
        if len(order) < self._dim:
            raise RuntimeError("Size of 'order' is smaller than problem dimension!")
        
        if len(periodic_dimensions) < self._dim:
            raise RuntimeError("Size of 'periodic_dimensions' is smaller than problem dimension!")
        
        for i in range(self._dim):
            if order[i] not in [6,10]:
                raise RuntimeError("order[%d] has to be one of {6,10}" %i)
        
        self._order = tuple(order)
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
        
        self.der_x = None
        self.der_y = None
        self.der_z = None
        
        if self._order[0] == 6:
            self._der_x = pycd06.cd06stuff.cd06( self._nx, grid_spacing[0], self._periodic[0], 0, 0 )
        elif self._order[0] == 10:
            self._der_x = pycd10.cd10stuff.cd10( self._nx, grid_spacing[0], self._periodic[0], 0, 0 )
        
        if self._dim > 1:
            if self._order[1] == 6:
                self._der_y = pycd06.cd06stuff.cd06( self._ny, grid_spacing[1], self._periodic[1], 0, 0 )
            elif self._order[1] == 10:
                self._der_y = pycd10.cd10stuff.cd10( self._ny, grid_spacing[1], self._periodic[1], 0, 0 )
        
        if self._dim > 2:
            if self._order[2] == 6:
                self._der_z = pycd06.cd06stuff.cd06( self._nz, grid_spacing[2], self._periodic[2], 0, 0 )
            elif self._order[2] == 10:
                self._der_z = pycd10.cd10stuff.cd10( self._nz, grid_spacing[2], self._periodic[2], 0, 0 )
        
        # Initialize the data reshaper.
        self._data_reshaper = data_reshaper.DataReshaper(self._dim, data_order='F')
    
    
    def ddx(self, data, der=None, component_idx=None, bc=(0,0)):
        """
        Function to compute the first order derivative of data in first direction.
        
        data : input numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        der : optional output numpy array in Fortran contiguous layout. This array must be consistent with the 3D
              decomposition and the problem dimension. This method will return der if der is None
        component_idx : index of component in data to for taking derivative. None if there is only one component in the
                        data
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
        
        return_der = True
        if der is None:
            der = numpy.empty(self._chunk_3d_size[0:self._dim], dtype=numpy.float64, order='F')
        else:
            return_der = False
        
        data_x = numpy.empty(self._chunk_x_size, dtype=numpy.float64, order='F')
        der_x  = numpy.empty(self._chunk_x_size, dtype=numpy.float64, order='F')
        
        data_3d = []
        if component_idx is None:
            data_3d = self._data_reshaper.reshapeTo3d(data)
        else:
            data_3d = self._data_reshaper.reshapeTo3d(data, component_idx)
        
        self._grid_partition.transpose_3d_to_x(data_3d, data_x)
        
        der_3d = self._data_reshaper.reshapeTo3d(der)
        
        if self._order[0] == 6:
            # symmetry BC only supported in 10th order for now
            self._der_x.dd1(data_x, der_x, self._chunk_x_size[1], self._chunk_x_size[2])
        elif self._order[0] == 10:
            self._der_x.dd1(data_x, der_x, self._chunk_x_size[1], self._chunk_x_size[2], bc1_=bc[0], bcn_=bc[1])
        
        self._grid_partition.transpose_x_to_3d(der_x, der_3d)
        der = self._data_reshaper.reshapeFrom3d(der_3d)
        
        if return_der:
            return der
    
    
    def ddy(self, data, der=None, component_idx=None, bc=(0,0)):
        """
        Function to compute the first order derivative of data in second direction.
        
        data : input numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        der : optional output numpy array in Fortran contiguous layout. This array must be consistent with the 3D
              decomposition and the problem dimension. This method will return der if der is None
        component_idx : index of component in data to for taking derivative. None if there is only one component in the
                        data
        bc : integer iterable of size 2 with the boundary condition at the left and right.
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
        
        return_der = True
        if der is None:
            der = numpy.empty(self._chunk_3d_size[0:self._dim], dtype=numpy.float64, order='F')
        else:
            return_der = False
        
        data_y = numpy.empty(self._chunk_y_size, dtype=numpy.float64, order='F')
        der_y  = numpy.empty(self._chunk_y_size, dtype=numpy.float64, order='F')
        
        data_3d = []
        if component_idx is None:
            data_3d = self._data_reshaper.reshapeTo3d(data)
        else:
            data_3d = self._data_reshaper.reshapeTo3d(data, component_idx)
        
        self._grid_partition.transpose_3d_to_y(data_3d, data_y)
        
        der_3d = self._data_reshaper.reshapeTo3d(der)
        
        if self._order[0] == 6:
            # symmetry BC only supported in 10th order for now
            self._der_y.dd2(data_y, der_y, self._chunk_y_size[0], self._chunk_y_size[2])
        elif self._order[0] == 10:
            self._der_y.dd2(data_y, der_y, self._chunk_y_size[0], self._chunk_y_size[2], bc1_=bc[0], bcn_=bc[1])
        
        self._grid_partition.transpose_y_to_3d(der_y, der_3d)
        der = self._data_reshaper.reshapeFrom3d(der_3d)
        
        if return_der:
            return der
    
    
    def ddz(self, data, der=None, component_idx=None, bc=(0,0)):
        """
        Function to compute the first order derivative of data in third direction.
        
        data : input numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        der : optional output numpy array in Fortran contiguous layout. This array must be consistent with the 3D
              decomposition and the problem dimension. This method will return der if der is None
        component_idx : index of component in data to for taking derivative. None if there is only one component in the
                        data
        bc : integer iterable of size 2 with the boundary condition at the left and right.
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
        
        return_der = True
        if der is None:
            der = numpy.empty(self._chunk_3d_size, dtype=numpy.float64, order='F')
        else:
            return_der = False
        
        data_z = numpy.empty(self._chunk_z_size, dtype=numpy.float64, order='F')
        der_z  = numpy.empty(self._chunk_z_size, dtype=numpy.float64, order='F')
        
        data_3d = []
        if component_idx is None:
            data_3d = data
        else:
            data_3d = data[:, :, :, component_idx]
        
        self._grid_partition.transpose_3d_to_z(data_3d, data_z)
        
        if self._order[0] == 6:
            # symmetry BC only supported in 10th order for now
            self._der_z.dd3(data_z, der_z, self._chunk_z_size[0], self._chunk_z_size[1])
        elif self._order[0] == 10:
            self._der_z.dd3(data_z, der_z, self._chunk_z_size[0], self._chunk_z_size[1], bc1_=bc[0], bcn_=bc[1])
        
        self._grid_partition.transpose_z_to_3d(der_z, der)
        
        if return_der:
            return der
    
    
    def d2dx2(self, data, der=None, component_idx=None, bc=(0,0)):
        """
        Function to compute the second order derivative of data in first direction.
        
        data : input numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        der : optional output numpy array in Fortran contiguous layout. This array must be consistent with the 3D
              decomposition and the problem dimension. This method will return der if der is None
        component_idx : index of component in data to for taking derivative. None if there is only one component in the
                        data
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
        
        return_der = True
        if der is None:
            der = numpy.empty(self._chunk_3d_size[0:self._dim], dtype=numpy.float64, order='F')
        else:
            return_der = False
        
        data_x = numpy.empty(self._chunk_x_size, dtype=numpy.float64, order='F')
        der_x  = numpy.empty(self._chunk_x_size, dtype=numpy.float64, order='F')
        
        data_3d = []
        if component_idx is None:
            data_3d = self._data_reshaper.reshapeTo3d(data)
        else:
            data_3d = self._data_reshaper.reshapeTo3d(data, component_idx)
        
        self._grid_partition.transpose_3d_to_x(data_3d, data_x)
        
        der_3d = self._data_reshaper.reshapeTo3d(der)
        
        if self._order[0] == 6:
            raise NotImplementedError("6th order 2nd derivatives are not implemented yet. Sorry!")
        elif self._order[0] == 10:
            self._der_x.d2d1(data_x, der_x, self._chunk_x_size[1], self._chunk_x_size[2], bc1_=bc[0], bcn_=bc[1])
        
        self._grid_partition.transpose_x_to_3d(der_x, der_3d)
        der = self._data_reshaper.reshapeFrom3d(der_3d)
        
        if return_der:
            return der
    
    
    def d2dy2(self, data, der=None, component_idx=None, bc=(0,0)):
        """
        Function to compute the second order derivative of data in second direction.
        
        data : input numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        der : optional output numpy array in Fortran contiguous layout. This array must be consistent with the 3D
              decomposition and the problem dimension. This method will return der if der is None
        component_idx : index of component in data to for taking derivative. None if there is only one component in the
                        data
        bc : integer iterable of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """
        
        data_shape = data.shape
        
        if component_idx is not None:
            data_shape = data_shape[0:-1]
        
        if self._dim == 1:
            raise RuntimeError("There is no d2dy2 for 1D problem!")
        
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
        
        return_der = True
        if der is None:
            der = numpy.empty(self._chunk_3d_size[0:self._dim], dtype=numpy.float64, order='F')
        else:
            return_der = False
        
        data_y = numpy.empty(self._chunk_y_size, dtype=numpy.float64, order='F')
        der_y  = numpy.empty(self._chunk_y_size, dtype=numpy.float64, order='F')
        
        data_3d = []
        if component_idx is None:
            data_3d = self._data_reshaper.reshapeTo3d(data)
        else:
            data_3d = self._data_reshaper.reshapeTo3d(data, component_idx)
        
        self._grid_partition.transpose_3d_to_y(data_3d, data_y)
        
        der_3d = self._data_reshaper.reshapeTo3d(der)
        
        if self._order[0] == 6:
            raise NotImplementedError("6th order 2nd derivatives are not implemented yet. Sorry!")
        elif self._order[0] == 10:
            self._der_y.d2d2(data_y, der_y, self._chunk_y_size[0], self._chunk_y_size[2], bc1_=bc[0], bcn_=bc[1])
        
        self._grid_partition.transpose_y_to_3d(der_y, der_3d)
        der = self._data_reshaper.reshapeFrom3d(der_3d)
        
        if return_der:
            return der
    
    
    def d2dz2(self, data, der=None, component_idx=None, bc=(0,0)):
        """
        Function to compute the second order derivative of data in third direction.
        
        data : input numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        der : optional output numpy array in Fortran contiguous layout. This array must be consistent with the 3D
              decomposition and the problem dimension. This method will return der if der is None
        component_idx : index of component in data to for taking derivative. None if there is only one component in the
                        data
        bc : integer iterable of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """
        
        data_shape = data.shape
        
        if component_idx is not None:
            data_shape = data_shape[0:-1]
        
        if self._dim == 1:
            raise RuntimeError("There is no d2dz2 for 1D problem!")
        
        elif self._dim == 2:
            raise RuntimeError("There is no d2dz2 for 1D problem!")
        
        else:
            if len(data_shape) != 3:
                raise RuntimeError("Make sure data is 3D!")
            if data_shape[0] != self._chunk_3d_size[0] or \
               data_shape[1] != self._chunk_3d_size[1] or \
               data_shape[2] != self._chunk_3d_size[2]:
                raise RuntimeError("Make sure data is of the same size as in grid_partition!")
        
        return_der = True
        if der is None:
            der = numpy.empty(self._chunk_3d_size, dtype=numpy.float64, order='F')
        else:
            return_der = False
        
        data_z = numpy.empty(self._chunk_z_size, dtype=numpy.float64, order='F')
        der_z  = numpy.empty(self._chunk_z_size, dtype=numpy.float64, order='F')
        
        data_3d = []
        if component_idx is None:
            data_3d = data
        else:
            data_3d = data[:, :, :, component_idx]
        
        self._grid_partition.transpose_3d_to_z(data_3d, data_z)
        
        if self._order[0] == 6:
            raise NotImplementedError("6th order 2nd derivatives are not implemented yet. Sorry!")
        elif self._order[0] == 10:
            self._der_z.d2d3(data_z, der_z, self._chunk_z_size[0], self._chunk_z_size[1], bc1_=bc[0], bcn_=bc[1])
        
        self._grid_partition.transpose_z_to_3d(der_z, der)
        
        if return_der:
            return der
    
    
    def gradient(self, data, component_idx=None, x_bc=(0,0), y_bc=(0,0), z_bc=(0,0)):
        """
        Function to compute the gradient of data.

        data : input numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        component_idx : index of component in data to for taking derivative. None if there is only one component in the
                        data
        *_bc : integer tuple of size 2 with the boundary condition at the left and right.
               0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        gradient_* : retruned output numpy array in Fortran contiguous layout. This array must be consistent with the
                     3D decomposition and the problem dimension. This method will return der if der is None
        """
        
        if self._dim == 1:
            return self.ddx(data, component_idx=component_idx, bc=x_bc)
        
        elif self._dim == 2:
            return self.ddx(data, component_idx=component_idx, bc=x_bc), \
                   self.ddy(data, component_idx=component_idx, bc=y_bc)
        
        else:
            return self.ddx(data, component_idx=component_idx, bc=x_bc), \
                   self.ddy(data, component_idx=component_idx, bc=y_bc), \
                   self.ddz(data, component_idx=component_idx, bc=z_bc)
    
    
    def divergence(self, data, x_bc=(0,0), y_bc=(0,0), z_bc=(0,0)):
        """
        Function to compute the gradient of a vector (u,v,w)

        data : input numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension. The number of components should be as same as the number of dimensions
        divergence : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
                     This array is in the 3D decomposition
        *_bc : integer tuple of size 2 with the boundary condition at the left and right.
               0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """
        
        data_shape = data.shape
        
        if self._dim == 1 and len(data_shape) != 2:
            raise RuntimeError("Make sure data is 1D and has enough number of components!")
        
        if self._dim == 2 and len(data_shape) != 3:
            raise RuntimeError("Make sure data is 2D and has enough number of components!")
        
        if self._dim == 3 and len(data_shape) != 4:
            raise RuntimeError("Make sure data is 3D and has enough number of components!")
        
        divergence = None
        
        if self._dim >= 1:
            divergence = self.ddx(data, component_idx=0, bc=x_bc)
        if self._dim >= 2:
            divergence = divergence + self.ddy(data, component_idx=1, bc=y_bc)
        if self._dim == 3:
            divergence = divergence + self.ddz(data, component_idx=2, bc=z_bc)
        
        return divergence
    
    
    def laplacian(self, data, component_idx=None, x_bc=(0,0), y_bc=(0,0), z_bc=(0,0)):
        """
        Function to compute the laplacian of data.

        data : input numpy array in Fortran contiguous layout. This array must be consistent with the 3D decomposition
               and the problem dimension
        component_idx : index of component in data to for taking derivative. None if there is only one component in the
                        data
        *_bc : integer tuple of size 2 with the boundary condition at the left and right.
               0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        laplacian_* : retruned output numpy array in Fortran contiguous layout. This array must be consistent with the
                      3D decomposition and the problem dimension. This method will return der if der is None
        """
        
        laplacian = None
        
        if self._dim >= 1:
            laplacian = self.d2dx2(data, component_idx=component_idx, bc=x_bc)
        if self._dim >= 2:
            laplacian = laplacian + self.d2dy2(data, component_idx=component_idx, bc=y_bc)
        if self._dim == 3:
            laplacian = laplacian + self.d2dz2(data, component_idx=component_idx, bc=z_bc)
        
        return laplacian
