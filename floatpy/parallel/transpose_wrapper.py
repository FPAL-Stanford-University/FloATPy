import numpy

from floatpy.readers.parallel_reader import ParallelDataReader

class TransposeWrapper(object):
    """
    Class to transpose data with parallel communication. Only data in Fortran order can be used.
    """
    
    def __init__(self, parallel_reader, direction):
        """
        Constructor of the class.
        
        parallel_reader : a concrete object of ParallelDataReader
        direction : direction to take transpose
        """
        
        if not isinstance(parallel_reader, ParallelDataReader):
            raise RuntimeError("The given data reader is not instance of the parallel data reader!")
        
        self._parallel_reader = parallel_reader
        
        if self._parallel_reader.dimension < 2:
            raise RuntimeError('Only data with dimension greater than 1 can be transposed!.')
        
        if direction < 0 or direction > 2:
            raise RuntimeError('Direction < 0 or > 2 is invalid!')
        
        # Get the dimension of data.
        self._dim = self._parallel_reader.dimension
        
        if direction >= self._dim:
            raise RuntimeError('Direction to transpose is not allowed with the dimensinality of data!')
        
        self._direction = direction
        
        # Get size of sub-domain of this processor and lo and hi of the sub-domain.
        self._size = numpy.empty(3, dtype=numpy.int32)
        self._lo   = numpy.empty(3, dtype=numpy.int32)
        self._hi   = numpy.empty(3, dtype=numpy.int32)
        
        self._grid_partition = self._parallel_reader.grid_partition
        
        if direction == 0:
            self._grid_partition.get_szx(self._size)
            self._grid_partition.get_stx(self._lo)
            self._grid_partition.get_enx(self._hi)
        
        elif direction == 1:
            self._grid_partition.get_szy(self._size)
            self._grid_partition.get_sty(self._lo)
            self._grid_partition.get_eny(self._hi)
        
        else:
            self._grid_partition.get_szz(self._size)
            self._grid_partition.get_stz(self._lo)
            self._grid_partition.get_enz(self._hi)
        
        # Convert to 0 based indexing.
        self._lo = self._lo - 1
        self._hi = self._hi - 1
    
    
    @property
    def full_chunk(self):
        """
        Return two tuples containing the full chunk of sub-domain after transpose used in the parallel reader as
        a lower bound (lo) and an upper bound (hi).
        """
        
        return tuple(self._lo[0:self._dim]), tuple(self._hi[0:self._dim])
    
    
    @property
    def full_chunk_size(self):
        """
        Return a tuple containing the size of the full chunk of sub-domain after transpose used in the parallel reader.
        """
        
        return tuple(self._size[0:self._dim])
    
    
    def transpose(self, data):
        """
        Transpose data.
        
        data : data to transpose
        """
        
        if not numpy.all(numpy.isreal(data)):
            raise ValueError("The given data is complex! Only real data can be transposed.")
        
        num_components = 1
        if data.ndim == self._dim + 1:
            num_components = data.shape[self._dim]
        
        shape_3D = data.shape[0:self._dim]
            
        if self._dim == 2:
            shape_3D = numpy.append(shape_3D, 1)
        
        data_out = []
        
        if num_components == 1:
            if self._dim == 2:
                data_out = numpy.empty((self._size[0], self._size[1]), dtype=data.dtype, order='F')
            else:
                data_out = numpy.empty((self._size[0], self._size[1], self._size[2]), dtype=data.dtype, order='F')
            
            data_3D = numpy.reshape(data, shape_3D, order='F')
            data_transposed = numpy.reshape(data_out, self._size, order='F')
            
            if self._direction == 0:
                self._grid_partition.transpose_3d_to_x(data_3D, data_transposed)
            elif self._direction == 1:
                self._grid_partition.transpose_3d_to_y(data_3D, data_transposed)
            else:
                self._grid_partition.transpose_3d_to_z(data_3D, data_transposed)
        
        else:
            if self._dim == 2:
                data_out = numpy.empty((self._size[0], self._size[1], num_components), dtype=data.dtype, order='F')
            else:
                data_out = numpy.empty((self._size[0], self._size[1], self._size[2], num_components), dtype=data.dtype, order='F')
            
            for ic in range(num_components):
                data_3D = []
                data_transposed = []
                
                if self._dim == 2:
                    data_3D = numpy.reshape(data[:, :, ic], shape_3D, order='F')
                    data_transposed = numpy.reshape(data_out[:, :, ic], self._size, order='F')
                else:
                    data_3D = data[:, :, :, ic]
                    data_transposed = data_out[:, :, :, ic]
                
                if self._direction == 0:
                    self._grid_partition.transpose_3d_to_x(data_3D, data_transposed)
                elif self._direction == 1:
                    self._grid_partition.transpose_3d_to_y(data_3D, data_transposed)
                else:
                    self._grid_partition.transpose_3d_to_z(data_3D, data_transposed)
        
        return data_out
