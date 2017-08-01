from floatpy.readers.parallel_reader import ParallelDataReader

class Transpose():
    """
    Class to transpose data with parallel communication.
    """
    
    def __init__(self, parallel_reader):
        """
        Constructor of the class.
        
        parallel_reader: a concrete object of ParallelDataReader
        """
        
        if not isinstance(parallel_reader, ParallelDataReader):
            raise RuntimeError("The given data reader is not instance of the parallel data reader!")
        
        self._parallel_reader = parallel_reader
        
        if self._parallel_reader.dimension < 2:
            raise RuntimeError('Only data with dimension greater than 1 can be transposed!.')
    
    
    def transpose(self, data, direction=0):
        """
        Transpose data.
        
        data:      data to transpose
        direction: direction to take transpose
        """
        
        if not numpy.all(numpy.isreal(data)):
            raise ValueError("The given data is complex! Only real data can be transposed.")
        
        if not numpy.isfortran(data):
            raise RuntimeError("The given data is not in Fortran order! Only data in Fortran order can be transposed.")
        
        if direction < 0 or direction > 2:
            raise RuntimeError('Direction < 0 or > 2 is invalid!')
        
        # Get the dimension of data.
        dim = self._parallel_reader.dimension
        
        if direction >= dim:
            raise RuntimeError('Direction to transpose is not allowed with the dimensinality of data!')
        
        # Get size of sub-domain of this processor and lo and up of the sub-domain.
        size = numpy.empty(3, dtype=numpy.int32)
        lo   = numpy.empty(3, dtype=numpy.int32)
        up   = numpy.empty(3, dtype=numpy.int32)
        
        gp = self._parallel_reader._grid_partition
        
        if direction == 0:
            gp.get_szx(size)
            gp.get_stx(lo)
            gp.get_enx(up)
        
        elif direction == 1:
            gp.get_szy(size)
            gp.get_sty(lo)
            gp.get_eny(up)
        
        else:
            gp.get_szz(size)
            gp.get_stz(lo)
            gp.get_enz(up)
        
        # Convert to 0 based indexing.
        lo = lo - 1
        up = up - 1
        
        num_components = 1
        if data.ndim == dim + 1:
            num_components = data.shape[dim]
        
        shape_3D = data.shape[0:dim]
            
        if dim == 2:
            shape_3D = numpy.append(shape_3D, 1)
        
        data_out = []
        
        if num_components == 1:
            if dim == 2:
                data_out = numpy.empty((size[0], size[1]), dtype=data.dtype, order='F')
            else:
                data_out = numpy.empty((size[0], size[1], size[2]), dtype=data.dtype, order='F')
            
            data_3D = numpy.reshape(data, shape_3D)
            data_transposed = numpy.reshape(data_out, size)
            
            if direction == 0:
                gp.transpose_3d_to_x(data_3D, data_transposed)
            elif direction == 1:
                gp.transpose_3d_to_y(data_3D, data_transposed)
            else:
                gp.transpose_3d_to_z(data_3D, data_transposed)
        
        else:
            if dim == 2:
                data_out = numpy.empty((size[0], size[1], num_components), dtype=data.dtype, order='F')
            else:
                data_out = numpy.empty((size[0], size[1], size[2], num_components), dtype=data.dtype, order='F')
            
            for ic in range(num_components):
                data_3D = []
                data_transposed = []
                
                if dim == 2:
                    data_3D = numpy.reshape(data[:, :, ic], shape_3D)
                    data_transposed = numpy.reshape(data_out[:, :, ic], size)
                else:
                    data_3D = data[:, :, :, ic]
                    data_transposed = data_out[:, :, :, ic]
                
                if direction == 0:
                    gp.transpose_3d_to_x(data_3D, data_transposed)
                elif direction == 1:
                    gp.transpose_3d_to_y(data_3D, data_transposed)
                else:
                    gp.transpose_3d_to_z(data_3D, data_transposed)
        
        return data_out, lo, up
