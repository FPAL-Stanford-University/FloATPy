import numpy

class DataReshaper(object):
    
    def __init__(self, dimension, data_order='F'):
        
        if dimension < 1 or dimension > 3:
            raise ValueError("Number of dimensions should be between 1 and 3!")
        
        self._dim = dimension
        
        if data_order != 'C' and \
           data_order != 'F':
            raise RuntimeError("Invalid data order! Data order can only be 'C' or 'F'.")
        
        self._data_order = data_order
    
        
    def reshapeTo3d(self, data, component_idx=None, data_output=None):
        """
        Reshape data to 3D.
        
        data : numpy data to reshape
        component_idx : integer representing component index of data to reshape. None if there is only one component
        data_output : optional output numpy array. This method will return data_output if it is not given
        """
        
        # Get the shape of data.
        
        data_shape = numpy.array(data.shape)
        
        # Check whether the shape of data is valid.
        
        if component_idx is None:
            if (self._dim == 1 and data.ndim != 1) or \
               (self._dim == 2 and data.ndim != 2) or \
               (self._dim == 3 and data.ndim != 3):
                raise RuntimeError('Dimension of data is invalid!')
        else:
            if (self._dim == 1 and data.ndim != 2) or \
               (self._dim == 2 and data.ndim != 3) or \
               (self._dim == 3 and data.ndim != 4):
                raise RuntimeError('Dimension of data is invalid!')
            
            # Check whether the component_idx is valid and get the shape of the component's data.
            
            if self._data_order == 'C':
                if component_idx >= data.shape[0] or component_idx < 0:
                    raise RuntimeError('Component index is invalid!')
                
                data_shape = numpy.array(data_shape[1:])
            
            else:
                if component_idx >= data.shape[-1] or component_idx < 0:
                    raise RuntimeError('Component index is invalid!')
                
                data_shape = numpy.array(data_shape[:-1])
        
        shape_3d = numpy.ones(3)
        shape_3d[:self._dim] = data_shape
        
        # Get the component's data.
        
        data_component = None
       
        if component_idx is None:
            data_component = data
        else:
            if self._data_order == 'C':
                if self._dim == 1:
                    data_component = data[component_idx, :]
                elif self._dim == 2:
                    data_component = data[component_idx, :, :]
                elif self._dim == 3:
                    data_component = data[component_idx, :, :, :]
            else:
                if self._dim == 1:
                    data_component = data[:, component_idx]
                elif self._dim == 2:
                    data_component = data[:, :, component_idx]
                elif self._dim == 3:
                    data_component = data[:, :, :, component_idx]
        
        if data_output is None:
            return numpy.reshape(data_component, shape_3d, order=self._data_order)
        else:
            data_output = numpy.reshape(data_component, shape_3d, order=self._data_order)
    
    
    def reshapeFrom3d(self, data, data_output=None):
        """
        Reshape data into low dimension from 3D.
        
        data : numpy data to reshape
        data_output : optional output numpy array. This method will return data_output if it is not given
        """
        
        # Get the shape of data.
        
        shape_low_dim = numpy.array(data.shape)
        
        # Check whether the shape of data is valid.
        
        if component_idx is None:
            if data.ndim != 3:
                raise RuntimeError('Dimension of data is invalid!')
            
            if not (numpy.all(data.shape[self._dim:] == 1)):
                raise RuntimeError('Shape of data is invalid!')
            
            shape_low_dim = numpy.array(shape_low_dim[:self._dim])
        
        else:
            if data.ndim != 4:
                raise RuntimeError('Dimension of data is invalid!')
            
            # Check whether the component_idx is valid and get the shape of the component's data.
            
            if self._data_order == 'C':
                if component_idx >= data.shape[0] or component_idx < 0:
                    raise RuntimeError('Component index is invalid!')
                
                if not (numpy.all(data.shape[self._dim+1:] == 1)):
                    raise RuntimeError('Shape of data is invalid!')
                
                shape_low_dim = numpy.array(shape_low_dim[:1+self._dim])
            
            else:
                if component_idx >= data.shape[-1] or component_idx < 0:
                    raise RuntimeError('Component index is invalid!')
                
                if not (numpy.all(data.shape[self._dim:-1] == 1)):
                    raise RuntimeError('Shape of data is invalid!')
                
                shape_low_dim = numpy.append(shape_low_dim[:self._dim], shape_low_dim[-1])
        
        if data_output is None:
            return numpy.reshape(data, shape_low_dim, order=self._data_order)
        else:
            data_output = numpy.reshape(data, shape_low_dim, order=self._data_order)
