import numpy

from floatpy.parallel import _t3dmod

class TransposeWrapper(object):
    """
    Class to transpose data to/from pencil with parallel communication. Only data in Fortran order can be used.
    """
    
    def __init__(self, grid_partition, direction, dimension=3):
        """
        Constructor of the class.
        
        grid_partition : t3d object or the grid_partition property of the parallel data reader class
        direction : direction of pencil
        """
        
        if not isinstance(grid_partition, _t3dmod.t3d):
            raise RuntimeError("The given grid partition object is not an instance of the t3d class!")
        
        self._grid_partition = grid_partition
        
        if direction < 0 or direction > 2:
            raise RuntimeError('Direction < 0 or > 2 is invalid!')
        
        if dimension < 2 or dimension > 3:
            raise RuntimeError('Only data with dimension of 2 or 3 can be transposed!')
        
        # Get the dimension of data.
        self._dim = dimension
        
        if direction >= self._dim:
            raise RuntimeError('Direction to transpose is not allowed with the dimensinality of data!')
        
        self._direction = direction
        
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
        
        # Get size of pencil of this processor and lo and hi of the pencil.
        self._pencil_size = numpy.empty(3, dtype=numpy.int32)
        self._pencil_lo   = numpy.empty(3, dtype=numpy.int32)
        self._pencil_hi   = numpy.empty(3, dtype=numpy.int32)
        
        if direction == 0:
            self._grid_partition.get_szx(self._pencil_size)
            self._grid_partition.get_stx(self._pencil_lo)
            self._grid_partition.get_enx(self._pencil_hi)
        
        elif direction == 1:
            self._grid_partition.get_szy(self._pencil_size)
            self._grid_partition.get_sty(self._pencil_lo)
            self._grid_partition.get_eny(self._pencil_hi)
        
        else:
            self._grid_partition.get_szz(self._pencil_size)
            self._grid_partition.get_stz(self._pencil_lo)
            self._grid_partition.get_enz(self._pencil_hi)
        
        # Convert to 0 based indexing.
        self._pencil_lo = self._pencil_lo - 1
        self._pencil_hi = self._pencil_hi - 1
    
    
    @property
    def full_pencil(self):
        """
        Return two tuples containing the full chunk of pencil as a lower bound (lo) and an upper bound (hi).
        """
        
        return tuple(self._pencil_lo[0:self._dim]), tuple(self._pencil_hi[0:self._dim])
    
    
    @property
    def full_pencil_size(self):
        """
        Return a tuple containing the size of the full chunk of pencil.
        """
        
        return tuple(self._pencil_size[0:self._dim])
    
    
    def transposeToPencil(self, data):
        """
        Transpose data to pencil.
        
        data : data to transpose
        """
        
        if not numpy.all(numpy.isreal(data)):
            raise ValueError("The given data is complex! Only real data can be transposed.")
        
        num_components = 1
        if data.ndim == self._dim + 1:
            num_components = data.shape[self._dim]
        
        shape_3d = data.shape[0:self._dim]
            
        if self._dim == 2:
            shape_3d = numpy.append(shape_3d, 1)
        
        data_out = []
        
        if num_components == 1:
            data_to_transpose = numpy.empty(self._pencil_size, dtype=data.dtype, order='F')
            
            data_3d = numpy.reshape(data, shape_3d, order='F')
            
            if self._direction == 0:
                self._grid_partition.transpose_3d_to_x(data_3d, data_to_transpose)
            elif self._direction == 1:
                self._grid_partition.transpose_3d_to_y(data_3d, data_to_transpose)
            else:
                self._grid_partition.transpose_3d_to_z(data_3d, data_to_transpose)
            
            if self._dim == 2:
                data_out = numpy.reshape(data_to_transpose, (self._pencil_size[0], self._pencil_size[1]), order='F')
            else:
                data_out = data_to_transpose
        
        
        else:
            data_to_transpose = numpy.empty(numpy.append(shape_3d, num_components), dtype=data.dtype, order='F')
            print data_to_transpose.shape
            
            for ic in range(num_components):
                data_3d = []
                if self._dim == 2:
                    data_3d = numpy.reshape(data[:, :, ic], shape_3d, order='F')
                else:
                    data_3d = data[:, :, :, ic]
                
                if self._direction == 0:
                    self._grid_partition.transpose_3d_to_x(data_3d, data_to_transpose[:, :, :, ic])
                elif self._direction == 1:
                    self._grid_partition.transpose_3d_to_y(data_3d, data_to_transpose[:, :, :, ic])
                else:
                    self._grid_partition.transpose_3d_to_z(data_3d, data_to_transpose[:, :, :, ic])
            
            if self._dim == 2:
                data_out = numpy.reshape(data_to_transpose, (self._pencil_size[0], self._pencil_size[1], num_components), \
                                         order='F')
            else:
                data_out = data_to_transpose
        
        return data_out
    
    
    def transposeFromPencil(self, data):
        """
        Transpose data from pencil.
        
        data : data to transpose
        """
        
        if not numpy.all(numpy.isreal(data)):
            raise ValueError("The given data is complex! Only real data can be transposed.")
        
        num_components = 1
        if data.ndim == self._dim + 1:
            num_components = data.shape[self._dim]
        
        shape_3d = data.shape[0:self._dim]
            
        if self._dim == 2:
            shape_3d = numpy.append(shape_3d, 1)
        
        data_out = []
        
        if num_components == 1:
            if self._dim == 2:
                data_out = numpy.empty((self._3d_size[0], self._3d_size[1]), dtype=data.dtype, order='F')
            else:
                data_out = numpy.empty((self._3d_size[0], self._3d_size[1], self._3d_size[2]), \
                                       dtype=data.dtype, order='F')
            
            data_3d = numpy.reshape(data, shape_3d, order='F')
            data_to_transpose = numpy.reshape(data_out, self._3d_size, order='F')
            
            if self._direction == 0:
                self._grid_partition.transpose_x_to_3d(data_3d, data_to_transpose)
            elif self._direction == 1:
                self._grid_partition.transpose_y_to_3d(data_3d, data_to_transpose)
            else:
                self._grid_partition.transpose_z_to_3d(data_3d, data_to_transpose)
        
        else:
            if self._dim == 2:
                data_out = numpy.empty((self._3d_size[0], self._3d_size[1], num_components), dtype=data.dtype, \
                                       order='F')
            else:
                data_out = numpy.empty((self._3d_size[0], self._3d_size[1], self._3d_size[2], num_components), \
                                       dtype=data.dtype, order='F')
            
            for ic in range(num_components):
                data_3d = []
                data_to_transpose = []
                
                if self._dim == 2:
                    data_3d = numpy.reshape(data[:, :, ic], shape_3d, order='F')
                    data_to_transpose = numpy.reshape(data_out[:, :, ic], self._3d_size, order='F')
                else:
                    data_3d = data[:, :, :, ic]
                    data_to_transpose = data_out[:, :, :, ic]
                
                if self._direction == 0:
                    self._grid_partition.transpose_x_to_3d(data_3d, data_to_transpose)
                elif self._direction == 1:
                    self._grid_partition.transpose_y_to_3d(data_3d, data_to_transpose)
                else:
                    self._grid_partition.transpose_z_to_3d(data_3d, data_to_transpose)
        
        return data_out
