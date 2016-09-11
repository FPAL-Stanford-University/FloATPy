"""
Module for computing first deriatives with finite differencing.
"""

import numpy

def computeSecondOrderFirstDerivative(data, dx, direction = 0, component_idx = 0, uses_one_sided = True):
    """
    Computing first derivative using explicit second order finite differencing.
    """
    
    data_shape = data.shape
    
    # Check whether the direction is valid.
    
    if direction < 0 or direction > 2:
        raise RuntimeError('Direction < 0 or > 2 is invalid!')
    
    # Check whether the dimension of data is valid.
    
    if data.ndim < 2:
        raise RuntimeError('Shape of data is invalid!')
    
    # Check whether the component_idx is valid.
    
    if component_idx >= data_shape[0] or component_idx < 0:
        raise RuntimeError('Component index is invalid!')
    
    # Check whether data size is large enough for second order first derivative.
    
    if direction == 0:
        if data_shape[1] < 3:
            raise RuntimeError('First dimension of data is not large enough!')
    
    elif direction == 1:
        if data_shape[2] < 3:
            raise RuntimeError('Second dimension of data is not large enough!')
        
    elif direction == 2:
        if data_shape[3] < 3:
            raise RuntimeError('Third dimension of data is not large enough!')
    
    # Initialize container to store the derivatives. The elements in the container
    # are initialized as NAN values.
    data_shape = numpy.delete(data_shape, [0])
    diff_data = numpy.empty(data_shape)
    diff_data[:] = numpy.NAN
    
    # Compute the derivatives in the interior of the domain.
    
    if direction == 0:
        if diff_data.ndim == 1:
            diff_data[1:-1] = (-1.0/2.0*data[component_idx, 0:-2] + 1.0/2.0*data[component_idx, 2:])/dx
        
        elif diff_data.ndim == 2:
            diff_data[1:-1, :] = (-1.0/2.0*data[component_idx, 0:-2, :] + 1.0/2.0*data[component_idx, 2:, :])/dx
        
        elif diff_data.ndim == 3:
            diff_data[1:-1, :, :] = (-1.0/2.0*data[component_idx, 0:-2, :, :] + 1.0/2.0*data[component_idx, 2:, :, :])/dx
        
        else:
            raise RuntimeError('Data dimension > 3 not supported!')
    
    elif direction == 1:
        if diff_data.ndim < 2:
            raise IOError('There is no second direction in data with less than two dimensions!')
        
        elif diff_data.ndim == 2:
            diff_data[:, 1:-1] = (-1.0/2.0*data[component_idx, :, 0:-2] + 1.0/2.0*data[component_idx, :, 2:])/dx
        
        elif diff_data.ndim == 3:
            diff_data[:, 1:-1, :] = (-1.0/2.0*data[component_idx, :, 0:-2, :] + 1.0/2.0*data[component_idx, :, 2:, :])/dx
        
        else:
            raise RuntimeError('Data dimension > 3 not supported!')
    
    elif direction == 2:
        if diff_data.ndim < 3:
            raise IOError('There is no third direction in data with less than three dimensions!')
        
        elif diff_data.ndim == 3:
            diff_data[:, :, 1:-1] = (-1.0/2.0*data[component_idx, :, :, 0:-2] + 1.0/2.0*data[component_idx, :, :, 2:])/dx
        
        else:
            raise RuntimeError('Data dimension > 3 not supported!')
    
    # Compute the derivatives at the boundaries.
    
    if uses_one_sided == True:
        if direction == 0:
            if diff_data.ndim == 1:
                diff_data[0] = (-3.0/2.0*data[component_idx, 0] + 2.0*data[component_idx, 1] \
                                - 1.0/2.0*data[component_idx, 2])/dx
                
                diff_data[-1] = (1.0/2.0*data[component_idx, -3] - 2.0*data[component_idx, -2] \
                                 + 3.0/2.0*data[component_idx, -1])/dx
            
            elif diff_data.ndim == 2:
                diff_data[0, :] = (-3.0/2.0*data[component_idx, 0, :] + 2.0*data[component_idx, 1, :] \
                                   - 1.0/2.0*data[component_idx, 2, :])/dx
                
                diff_data[-1, :] = (1.0/2.0*data[component_idx, -3, :] - 2.0*data[component_idx, -2, :] \
                                    + 3.0/2.0*data[component_idx, -1, :])/dx
            
            elif diff_data.ndim == 3:
                diff_data[0, :, :] = (-3.0/2.0*data[component_idx, 0, :, :] + 2.0*data[component_idx, 1, :, :] \
                                      - 1.0/2.0*data[component_idx, 2, :, :])/dx
                
                diff_data[-1:, :, :] = (1.0/2.0*data[component_idx, -3, :, :] - 2.0*data[component_idx, -2, :, :] \
                                        + 3.0/2.0*data[component_idx, -1, :, :])/dx
            
            else:
                raise RuntimeError('Data dimension > 3 not supported!')
        
        elif direction == 1:
            if diff_data.ndim < 2:
                raise RuntimeError('There is no second direction in data with less than two dimensions!')
            
            elif diff_data.ndim == 2:
                diff_data[:, 0] = (-3.0/2.0*data[component_idx, :, 0] + 2.0*data[component_idx, :, 1] \
                                   - 1.0/2.0*data[component_idx, :, 2])/dx
                
                diff_data[:, -1] = (1.0/2.0*data[component_idx, :, -3] - 2.0*data[component_idx, :, -2] \
                                    + 3.0/2.0*data[component_idx, :, -1])/dx
            
            elif diff_data.ndim == 3:
                diff_data[:, 0, :] = (-3.0/2.0*data[component_idx, :, 0, :] + 2.0*data[component_idx, :, 1, :] \
                                      - 1.0/2.0*data[component_idx, :, 2, :])/dx
                
                diff_data[:, -1, :] = (1.0/2.0*data[component_idx, :, -3, :] - 2.0*data[component_idx, :, -2, :] \
                                       + 3.0/2.0*data[component_idx, :, -1, :])/dx
            
            else:
                raise RuntimeError('Data dimension > 3 not supported!')
        
        elif direction == 2:
            if diff_data.ndim < 3:
                raise IOError('There is no third direction in data with less than three dimensions!')
            
            elif diff_data.ndim == 3:
                diff_data[:, :, 0] = (-3.0/2.0*data[component_idx, :, :, 0] + 2.0*data[component_idx, :, :, 1] \
                                      - 1.0/2.0*data[component_idx, :, :, 2])/dx
                
                diff_data[:, :, -1] = (1.0/2.0*data[component_idx, :, :, -3] - 2.0*data[component_idx, :, :, -2] \
                                       + 3.0/2.0*data[component_idx, :, :, -1])/dx
            
            else:
                raise RuntimeError('Data dimension > 3 not supported!')
    
    return diff_data


def computeFourthOrderFirstDerivative(data, dx, direction = 0, component_idx = 0, uses_one_sided = True):
    """
    Computing first derivative using explicit fourth order finite differencing.
    """
    
    data_shape = data.shape
    
    # Check whether the direction is valid.
    
    if direction < 0 or direction > 2:
        raise RuntimeError('Direction < 0 or > 2 is invalid!')
    
    # Check whether the dimension of data is valid.
    
    if data.ndim < 2:
        raise RuntimeError('Shape of data is invalid!')
    
    # Check whether the component_idx is valid.
    
    if component_idx >= data_shape[0] or component_idx < 0:
        raise RuntimeError('Component index is invalid!')
    
    # Check whether data size is large enough for fourth order first derivative.
    
    if direction == 0:
        if data_shape[1] < 5:
            raise RuntimeError('First dimension of data is not large enough!')
    
    elif direction == 1:
        if data_shape[2] < 5:
            raise RuntimeError('Second dimension of data is not large enough!')
        
    elif direction == 2:
        if data_shape[3] < 5:
            raise RuntimeError('Third dimension of data is not large enough!')
    
    # Initialize container to store the derivatives. The elements in the container
    # are initialized as NAN values.
    data_shape = numpy.delete(data_shape, [0])
    diff_data = numpy.empty(data_shape)
    diff_data[:] = numpy.NAN
    
    # Compute the derivatives in the interior of the domain.
    
    if direction == 0:
        if diff_data.ndim == 1:
            diff_data[2:-2] = (1.0/12.0*data[component_idx, 0:-4] - 2.0/3.0*data[component_idx, 1:-3] \
                               + 2.0/3.0*data[component_idx, 3:-1] - 1.0/12.0*data[component_idx, 4:])/dx
        
        elif diff_data.ndim == 2:
            diff_data[2:-2, :] = (1.0/12.0*data[component_idx, 0:-4, :] - 2.0/3.0*data[component_idx, 1:-3, :] \
                                  + 2.0/3.0*data[component_idx, 3:-1, :] - 1.0/12.0*data[component_idx, 4:, :])/dx
        
        elif diff_data.ndim == 3:
            diff_data[2:-2, :, :] = (1.0/12.0*data[component_idx, 0:-4, :, :] - 2.0/3.0*data[component_idx, 1:-3, :, :] \
                                     + 2.0/3.0*data[component_idx, 3:-1, :, :] - 1.0/12.0*data[component_idx, 4:, :, :])/dx
        else:
            raise RuntimeError('Data dimension > 3 not supported!')
    
    elif direction == 1:
        if diff_data.ndim < 2:
            raise IOError('There is no second direction in data with less than two dimensions!')
        
        elif diff_data.ndim == 2:
            diff_data[:, 2:-2] = (1.0/12.0*data[component_idx, :, 0:-4] - 2.0/3.0*data[component_idx, :, 1:-3] \
                                  + 2.0/3.0*data[component_idx, :, 3:-1] - 1.0/12.0*data[component_idx, :, 4:])/dx
        
        elif diff_data.ndim == 3:
            diff_data[:, 2:-2, :] = (1.0/12.0*data[component_idx, :, 0:-4, :] - 2.0/3.0*data[component_idx, :, 1:-3, :] \
                                     + 2.0/3.0*data[component_idx, :, 3:-1, :] - 1.0/12.0*data[component_idx, :, 4:, :])/dx
        
        else:
            raise RuntimeError('Data dimension > 3 not supported!')
    
    elif direction == 2:
        if diff_data.ndim < 3:
            raise IOError('There is no third direction in data with less than three dimensions!')
        
        elif diff_data.ndim == 3:
            diff_data[:, :, 2:-2] = (1.0/12.0*data[component_idx, :, :, 0:-4] - 2.0/3.0*data[component_idx, :, :, 1:-3] \
                                     + 2.0/3.0*data[component_idx, :, :, 3:-1] - 1.0/12.0*data[component_idx, :, :, 4:])/dx
        
        else:
            raise RuntimeError('Data dimension > 3 not supported!')
    
    # Compute the derivatives at the boundaries.
    
    if uses_one_sided == True:
        if direction == 0:
            if diff_data.ndim == 1:
                diff_data[0] = (-25.0/12.0*data[component_idx, 0] + 4.0*data[component_idx, 1] \
                                - 3.0*data[component_idx, 2] + 4.0/3.0*data[component_idx, 3] \
                                - 1.0/4.0*data[component_idx, 4])/dx
                
                diff_data[1] = (-1.0/4.0*data[component_idx, 0] - 5.0/6.0*data[component_idx, 1] \
                                + 3.0/2.0*data[component_idx, 2] - 1.0/2.0*data[component_idx, 3] \
                                + 1.0/12.0*data[component_idx, 4])/dx
                
                diff_data[-2] = (-1.0/12.0*data[component_idx, -5] + 1.0/2.0*data[component_idx, -4] \
                                 - 3.0/2.0*data[component_idx, -3] + 5.0/6.0*data[component_idx, -2] \
                                 + 1.0/4.0*data[component_idx, -1])/dx
                
                diff_data[-1] = (1.0/4.0*data[component_idx, -5] - 4.0/3.0*data[component_idx, -4] \
                                 + 3.0*data[component_idx, -3] - 4.0*data[component_idx, -2] \
                                 + 25.0/12.0*data[component_idx, -1])/dx
            
            elif diff_data.ndim == 2:
                diff_data[0, :] = (-25.0/12.0*data[component_idx, 0, :] + 4.0*data[component_idx, 1, :] \
                                   - 3.0*data[component_idx, 2, :] + 4.0/3.0*data[component_idx, 3, :] \
                                   - 1.0/4.0*data[component_idx, 4, :])/dx
                
                diff_data[1, :] = (-1.0/4.0*data[component_idx, 0, :] - 5.0/6.0*data[component_idx, 1, :] \
                                   + 3.0/2.0*data[component_idx, 2, :] - 1.0/2.0*data[component_idx, 3, :] \
                                   + 1.0/12.0*data[component_idx, 4, :])/dx
                
                diff_data[-2, :] = (-1.0/12.0*data[component_idx, -5, :] + 1.0/2.0*data[component_idx, -4, :] \
                                    - 3.0/2.0*data[component_idx, -3, :] + 5.0/6.0*data[component_idx, -2, :] \
                                    + 1.0/4.0*data[component_idx, -1, :])/dx
                
                diff_data[-1, :] = (1.0/4.0*data[component_idx, -5, :] - 4.0/3.0*data[component_idx, -4, :] \
                                    + 3.0*data[component_idx, -3, :] - 4.0*data[component_idx, -2, :] \
                                    + 25.0/12.0*data[component_idx, -1, :])/dx
            
            elif diff_data.ndim == 3:
                diff_data[0, :, :] = (-25.0/12.0*data[component_idx, 0, :, :] + 4.0*data[component_idx, 1, :, :] \
                                      - 3.0*data[component_idx, 2, :, :] + 4.0/3.0*data[component_idx, 3, :, :] \
                                      - 1.0/4.0*data[component_idx, 4, :, :])/dx
                
                diff_data[1, :, :] = (-1.0/4.0*data[component_idx, 0, :, :] - 5.0/6.0*data[component_idx, 1, :, :] \
                                      + 3.0/2.0*data[component_idx, 2, :, :] - 1.0/2.0*data[component_idx, 3, :, :] \
                                      + 1.0/12.0*data[component_idx, 4, :, :])/dx
                
                diff_data[-2, :, :] = (-1.0/12.0*data[component_idx, -5, :, :] + 1.0/2.0*data[component_idx, -4, :, :] \
                                       - 3.0/2.0*data[component_idx, -3, :, :] + 5.0/6.0*data[component_idx, -2, :, :] \
                                       + 1.0/4.0*data[component_idx, -1, :, :])/dx
                
                diff_data[-1, :, :] = (1.0/4.0*data[component_idx, -5, :, :] - 4.0/3.0*data[component_idx, -4, :, :] \
                                       + 3.0*data[component_idx, -3, :, :] - 4.0*data[component_idx, -2, :, :] \
                                       + 25.0/12.0*data[component_idx, -1, :, :])/dx
            
            else:
                raise RuntimeError('Data dimension > 3 not supported!')
        
        elif direction == 1:
            if diff_data.ndim < 2:
                raise RuntimeError('There is no second direction in data with less than two dimensions!')
            
            elif diff_data.ndim == 2:
                diff_data[:, 0] = (-25.0/12.0*data[component_idx, :, 0] + 4.0*data[component_idx, :, 1] \
                                   - 3.0*data[component_idx, :, 2] + 4.0/3.0*data[component_idx, :, 3] \
                                   - 1.0/4.0*data[component_idx, :, 4])/dx
                
                diff_data[:, 1] = (-1.0/4.0*data[component_idx, :, 0] - 5.0/6.0*data[component_idx, :, 1] \
                                   + 3.0/2.0*data[component_idx, :, 2] - 1.0/2.0*data[component_idx, :, 3] \
                                   + 1.0/12.0*data[component_idx, :, 4])/dx
                
                diff_data[:, -2] = (-1.0/12.0*data[component_idx, :, -5] + 1.0/2.0*data[component_idx, :, -4] \
                                    - 3.0/2.0*data[component_idx, :, -3] + 5.0/6.0*data[component_idx, :, -2] \
                                    + 1.0/4.0*data[component_idx, :, -1])/dx
                
                diff_data[:, -1] = (1.0/4.0*data[component_idx, :, -5] - 4.0/3.0*data[component_idx, :, -4] \
                                    + 3.0*data[component_idx, :, -3] - 4.0*data[component_idx, :, -2] \
                                    + 25.0/12.0*data[component_idx, :, -1])/dx
            
            elif diff_data.ndim == 3:
                diff_data[:, 0, :] = (-25.0/12.0*data[component_idx, :, 0, :] + 4.0*data[component_idx, :, 1, :] \
                                      - 3.0*data[component_idx, :, 2, :] + 4.0/3.0*data[component_idx, :, 3, :] \
                                      - 1.0/4.0*data[component_idx, :, 4, :])/dx
                
                diff_data[:, 1, :] = (-1.0/4.0*data[component_idx, :, 0, :] - 5.0/6.0*data[component_idx, :, 1, :] \
                                      + 3.0/2.0*data[component_idx, :, 2, :] - 1.0/2.0*data[component_idx, :, 3, :] \
                                      + 1.0/12.0*data[component_idx, :, 4, :])/dx
                
                diff_data[:, -2, :] = (-1.0/12.0*data[component_idx, :, -5, :] + 1.0/2.0*data[component_idx, :, -4, :] \
                                       - 3.0/2.0*data[component_idx, :, -3, :] + 5.0/6.0*data[component_idx, :, -2, :] \
                                       + 1.0/4.0*data[component_idx, :, -1, :])/dx
                
                diff_data[:, -1, :] = (1.0/4.0*data[component_idx, :, -5, :] - 4.0/3.0*data[component_idx, :, -4, :] \
                                       + 3.0*data[component_idx, :, -3, :] - 4.0*data[component_idx, :, -2, :] \
                                       + 25.0/12.0*data[component_idx, :, -1, :])/dx
            
            else:
                raise RuntimeError('Data dimension > 3 not supported!')
        
        elif direction == 2:
            if diff_data.ndim < 3:
                raise IOError('There is no third direction in data with less than three dimensions!')
            
            elif diff_data.ndim == 3:
                diff_data[:, :, 0] = (-25.0/12.0*data[component_idx, :, :, 0] + 4.0*data[component_idx, :, :, 1] \
                                      - 3.0*data[component_idx, :, :, 2] + 4.0/3.0*data[component_idx, :, :, 3] \
                                      - 1.0/4.0*data[component_idx, :, :, 4])/dx
                
                diff_data[:, :, 1] = (-1.0/4.0*data[component_idx, :, :, 0] - 5.0/6.0*data[component_idx, :, :, 1] \
                                      + 3.0/2.0*data[component_idx, :, :, 2] - 1.0/2.0*data[component_idx, :, :, 3] \
                                      + 1.0/12.0*data[component_idx, :, :, 4])/dx
                
                diff_data[:, :, -2] = (-1.0/12.0*data[component_idx, :, :, -5] + 1.0/2.0*data[component_idx, :, :, -4] \
                                       - 3.0/2.0*data[component_idx, :, :, -3] + 5.0/6.0*data[component_idx, :, :, -2] \
                                       + 1.0/4.0*data[component_idx, :, :, -1])/dx
                
                diff_data[:, :, -1] = (1.0/4.0*data[component_idx, :, :, -5] - 4.0/3.0*data[component_idx, :, :, -4] \
                                       + 3.0*data[component_idx, :, :, -3] - 4.0*data[component_idx, :, :, -2] \
                                       + 25.0/12.0*data[component_idx, :, :, -1])/dx
            
            else:
                raise RuntimeError('Data dimension > 3 not supported!')
    
    return diff_data


def computeSixthOrderFirstDerivative(data, dx, direction = 0, component_idx = 0, uses_one_sided = True):
    """
    Computing first derivative using explicit sixth order finite differencing.
    """
    
    data_shape = data.shape
    
    # Check whether the direction is valid.
    
    if direction < 0 or direction > 2:
        raise RuntimeError('Direction < 0 or > 2 is invalid!')
    
    # Check whether the dimension of data is valid.
    
    if data.ndim < 2:
        raise RuntimeError('Shape of data is invalid!')
    
    # Check whether the component_idx is valid.
    
    if component_idx >= data_shape[0] or component_idx < 0:
        raise RuntimeError('Component index is invalid!')
    
    # Check whether data size is large enough for sixth order first derivative.
    
    if direction == 0:
        if data_shape[1] < 7:
            raise RuntimeError('First dimension of data is not large enough!')
    
    elif direction == 1:
        if data_shape[2] < 7:
            raise RuntimeError('Second dimension of data is not large enough!')
        
    elif direction == 2:
        if data_shape[3] < 7:
            raise RuntimeError('Third dimension of data is not large enough!')
    
    # Initialize container to store the derivatives. The elements in the container
    # are initialized as NAN values.
    data_shape = numpy.delete(data_shape, [0])
    diff_data = numpy.empty(data_shape)
    diff_data[:] = numpy.NAN
    
    # Compute the derivatives in the interior of the domain.
    
    if direction == 0:
        if diff_data.ndim == 1:
            diff_data[3:-3] = (-1.0/60.0*data[component_idx, 0:-6] + 3.0/20.0*data[component_idx, 1:-5] \
                               - 3.0/4.0*data[component_idx, 2:-4] + 3.0/4.0*data[component_idx, 4:-2] \
                               - 3.0/20.0*data[component_idx, 5:-1] + 1.0/60.0*data[component_idx, 6:])/dx
        
        elif diff_data.ndim == 2:
            diff_data[3:-3, :] = (-1.0/60.0*data[component_idx, 0:-6, :] + 3.0/20.0*data[component_idx, 1:-5, :] \
                                  - 3.0/4.0*data[component_idx, 2:-4, :] + 3.0/4.0*data[component_idx, 4:-2, :] \
                                  - 3.0/20.0*data[component_idx, 5:-1, :] + 1.0/60.0*data[component_idx, 6:, :])/dx
        
        elif diff_data.ndim == 3:
            diff_data[3:-3, :, :] = (-1.0/60.0*data[component_idx, 0:-6, :, :] + 3.0/20.0*data[component_idx, 1:-5, :, :] \
                                     - 3.0/4.0*data[component_idx, 2:-4, :, :] + 3.0/4.0*data[component_idx, 4:-2, :, :] \
                                     - 3.0/20.0*data[component_idx, 5:-1, :, :] + 1.0/60.0*data[component_idx, 6:, :, :])/dx
        
        else:
            raise RuntimeError('Data dimension > 3 not supported!')
    
    elif direction == 1:
        if diff_data.ndim < 2:
            raise IOError('There is no second direction in data with less than two dimensions!')
        
        elif diff_data.ndim == 2:
            diff_data[:, 3:-3] = (-1.0/60.0*data[component_idx, :, 0:-6] + 3.0/20.0*data[component_idx, :, 1:-5] \
                                  - 3.0/4.0*data[component_idx, :, 2:-4] + 3.0/4.0*data[component_idx, :, 4:-2] \
                                  - 3.0/20.0*data[component_idx, :, 5:-1] + 1.0/60.0*data[component_idx, :, 6:])/dx
        
        elif diff_data.ndim == 3:
            diff_data[:, 3:-3, :] = (-1.0/60.0*data[component_idx, :, 0:-6, :] + 3.0/20.0*data[component_idx, :, 1:-5, :] \
                                     - 3.0/4.0*data[component_idx, :, 2:-4, :] + 3.0/4.0*data[component_idx, :, 4:-2, :] \
                                     - 3.0/20.0*data[component_idx, :, 5:-1, :] + 1.0/60.0*data[component_idx, :, 6:, :])/dx
        
        else:
            raise RuntimeError('Data dimension > 3 not supported!')
    
    elif direction == 2:
        if diff_data.ndim < 3:
            raise IOError('There is no third direction in data with less than three dimensions!')
        
        elif diff_data.ndim == 3:
            diff_data[:, :, 3:-3] = (-1.0/60.0*data[component_idx, :, :, 0:-6] + 3.0/20.0*data[component_idx, :, :, 1:-5] \
                                     - 3.0/4.0*data[component_idx, :, :, 2:-4] + 3.0/4.0*data[component_idx, :, :, 4:-2] \
                                     - 3.0/20.0*data[component_idx, :, :, 5:-1] + 1.0/60.0*data[component_idx, :, :, 6:])/dx
        
        else:
            raise RuntimeError('Data dimension > 3 not supported!')
    
    # Compute the derivatives at the boundaries.
    
    if uses_one_sided == True:
        if direction == 0:
            if diff_data.ndim == 1:
                diff_data[0] = (-49.0/20.0*data[component_idx, 0] + 6.0*data[component_idx, 1] \
                                - 15.0/2.0*data[component_idx, 2] + 20.0/3.0*data[component_idx, 3] \
                                - 15.0/4.0*data[component_idx, 4] + 6.0/5.0*data[component_idx, 5] \
                                - 1.0/6.0*data[component_idx, 6])/dx
                
                diff_data[1] = (-1.0/6.0*data[component_idx, 0] - 77.0/60.0*data[component_idx, 1] \
                                + 5.0/2.0*data[component_idx, 2] - 5.0/3.0*data[component_idx, 3] \
                                + 5.0/6.0*data[component_idx, 4] - 1.0/4.0*data[component_idx, 5] \
                                + 1.0/30.0*data[component_idx, 6])/dx
                
                diff_data[2] = (1.0/30.0*data[component_idx, 0] - 2.0/5.0*data[component_idx, 1] \
                                - 7.0/12.0*data[component_idx, 2] + 4.0/3.0*data[component_idx, 3] \
                                - 1.0/2.0*data[component_idx, 4] + 2.0/15.0*data[component_idx, 5] \
                                - 1.0/60.0*data[component_idx, 6])/dx
                
                diff_data[-3] = (1.0/60.0*data[component_idx, -7] - 2.0/15.0*data[component_idx, -6] \
                                 + 1.0/2.0*data[component_idx, -5] - 4.0/3.0*data[component_idx, -4] \
                                 + 7.0/12.0*data[component_idx, -3] + 2.0/5.0*data[component_idx, -2] \
                                 - 1.0/30.0*data[component_idx, -1])/dx
                
                diff_data[-2] = (-1.0/30.0*data[component_idx, -7] + 1.0/4.0*data[component_idx, -6] \
                                 - 5.0/6.0*data[component_idx, -5] + 5.0/3.0*data[component_idx, -4] \
                                 - 5.0/2.0*data[component_idx, -3] + 77.0/60.0*data[component_idx, -2] \
                                 + 1.0/6.0*data[component_idx, -1])/dx
                
                diff_data[-1] = (1.0/6.0*data[component_idx, -7] - 6.0/5.0*data[component_idx, -6] \
                                 + 15.0/4.0*data[component_idx, -5] - 20.0/3.0*data[component_idx, -4] \
                                 + 15.0/2.0*data[component_idx, -3] - 6.0*data[component_idx, -2] \
                                 + 49.0/20.0*data[component_idx, -1])/dx
            
            elif diff_data.ndim == 2:
                diff_data[0, :] = (-49.0/20.0*data[component_idx, 0, :] + 6.0*data[component_idx, 1, :] \
                                   - 15.0/2.0*data[component_idx, 2, :] + 20.0/3.0*data[component_idx, 3, :] \
                                   - 15.0/4.0*data[component_idx, 4, :] + 6.0/5.0*data[component_idx, 5, :] \
                                   - 1.0/6.0*data[component_idx, 6, :])/dx
                
                diff_data[1, :] = (-1.0/6.0*data[component_idx, 0, :] - 77.0/60.0*data[component_idx, 1, :] \
                                   + 5.0/2.0*data[component_idx, 2, :] - 5.0/3.0*data[component_idx, 3, :] \
                                   + 5.0/6.0*data[component_idx, 4, :] - 1.0/4.0*data[component_idx, 5, :] \
                                   + 1.0/30.0*data[component_idx, 6, :])/dx
                
                diff_data[2, :] = (1.0/30.0*data[component_idx, 0, :] - 2.0/5.0*data[component_idx, 1, :] \
                                   - 7.0/12.0*data[component_idx, 2, :] + 4.0/3.0*data[component_idx, 3, :] \
                                   - 1.0/2.0*data[component_idx, 4, :] + 2.0/15.0*data[component_idx, 5, :] \
                                   - 1.0/60.0*data[component_idx, 6, :])/dx
                
                diff_data[-3, :] = (1.0/60.0*data[component_idx, -7, :] - 2.0/15.0*data[component_idx, -6, :] \
                                    + 1.0/2.0*data[component_idx, -5, :] - 4.0/3.0*data[component_idx, -4, :] \
                                    + 7.0/12.0*data[component_idx, -3, :] + 2.0/5.0*data[component_idx, -2, :] \
                                    - 1.0/30.0*data[component_idx, -1, :])/dx
                
                diff_data[-2, :] = (-1.0/30.0*data[component_idx, -7, :] + 1.0/4.0*data[component_idx, -6, :] \
                                    - 5.0/6.0*data[component_idx, -5, :] + 5.0/3.0*data[component_idx, -4, :] \
                                    - 5.0/2.0*data[component_idx, -3, :] + 77.0/60.0*data[component_idx, -2, :] \
                                    + 1.0/6.0*data[component_idx, -1, :])/dx
                
                diff_data[-1, :] = (1.0/6.0*data[component_idx, -7, :] - 6.0/5.0*data[component_idx, -6, :] \
                                    + 15.0/4.0*data[component_idx, -5, :] - 20.0/3.0*data[component_idx, -4, :] \
                                    + 15.0/2.0*data[component_idx, -3, :] - 6.0*data[component_idx, -2, :] \
                                    + 49.0/20.0*data[component_idx, -1, :])/dx
            
            elif diff_data.ndim == 3:
                diff_data[0, :, :] = (-49.0/20.0*data[component_idx, 0, :, :] + 6.0*data[component_idx, 1, :, :] \
                                      - 15.0/2.0*data[component_idx, 2, :, :] + 20.0/3.0*data[component_idx, 3, :, :] \
                                      - 15.0/4.0*data[component_idx, 4, :, :] + 6.0/5.0*data[component_idx, 5, :, :] \
                                      - 1.0/6.0*data[component_idx, 6, :, :])/dx
                
                diff_data[1, :, :] = (-1.0/6.0*data[component_idx, 0, :, :] - 77.0/60.0*data[component_idx, 1, :, :] \
                                      + 5.0/2.0*data[component_idx, 2, :, :] - 5.0/3.0*data[component_idx, 3, :, :] \
                                      + 5.0/6.0*data[component_idx, 4, :, :] - 1.0/4.0*data[component_idx, 5, :, :] \
                                      + 1.0/30.0*data[component_idx, 6, :, :])/dx
                
                diff_data[2, :, :] = (1.0/30.0*data[component_idx, 0, :, :] - 2.0/5.0*data[component_idx, 1, :, :] \
                                      - 7.0/12.0*data[component_idx, 2, :, :] + 4.0/3.0*data[component_idx, 3, :, :] \
                                      - 1.0/2.0*data[component_idx, 4, :, :] + 2.0/15.0*data[component_idx, 5, :, :] \
                                      - 1.0/60.0*data[component_idx, 6, :, :])/dx
                
                diff_data[-3, :, :] = (1.0/60.0*data[component_idx, -7, :, :] - 2.0/15.0*data[component_idx, -6, :, :] \
                                       + 1.0/2.0*data[component_idx, -5, :, :] - 4.0/3.0*data[component_idx, -4, :, :] \
                                       + 7.0/12.0*data[component_idx, -3, :, :] + 2.0/5.0*data[component_idx, -2, :, :] \
                                       - 1.0/30.0*data[component_idx, -1, :, :])/dx
                
                diff_data[-2, :, :] = (-1.0/30.0*data[component_idx, -7, :, :] + 1.0/4.0*data[component_idx, -6, :, :] \
                                       - 5.0/6.0*data[component_idx, -5, :, :] + 5.0/3.0*data[component_idx, -4, :, :] \
                                       - 5.0/2.0*data[component_idx, -3, :, :] + 77.0/60.0*data[component_idx, -2, :, :] \
                                       + 1.0/6.0*data[component_idx, -1, :, :])/dx
                
                diff_data[-1, :, :] = (1.0/6.0*data[component_idx, -7, :, :] - 6.0/5.0*data[component_idx, -6, :, :] \
                                       + 15.0/4.0*data[component_idx, -5, :, :] - 20.0/3.0*data[component_idx, -4, :, :] \
                                       + 15.0/2.0*data[component_idx, -3, :, :] - 6.0*data[component_idx, -2, :, :] \
                                       + 49.0/20.0*data[component_idx, -1, :, :])/dx
            
            else:
                raise RuntimeError('Data dimension > 3 not supported!')
        
        elif direction == 1:
            if diff_data.ndim < 2:
                raise RuntimeError('There is no second direction in data with less than two dimensions!')
            
            elif diff_data.ndim == 2:
                diff_data[:, 0] = (-49.0/20.0*data[component_idx, :, 0] + 6.0*data[component_idx, :, 1] \
                                   - 15.0/2.0*data[component_idx, :, 2] + 20.0/3.0*data[component_idx, :, 3] \
                                   - 15.0/4.0*data[component_idx, :, 4] + 6.0/5.0*data[component_idx, :, 5] \
                                   - 1.0/6.0*data[component_idx, :, 6])/dx
                
                diff_data[:, 1] = (-1.0/6.0*data[component_idx, :, 0] - 77.0/60.0*data[component_idx, :, 1] \
                                   + 5.0/2.0*data[component_idx, :, 2] - 5.0/3.0*data[component_idx, :, 3] \
                                   + 5.0/6.0*data[component_idx, :, 4] - 1.0/4.0*data[component_idx, :, 5] \
                                   + 1.0/30.0*data[component_idx, :, 6])/dx
                
                diff_data[:, 2] = (1.0/30.0*data[component_idx, :, 0] - 2.0/5.0*data[component_idx, :, 1] \
                                   - 7.0/12.0*data[component_idx, :, 2] + 4.0/3.0*data[component_idx, :, 3] \
                                   - 1.0/2.0*data[component_idx, :, 4] + 2.0/15.0*data[component_idx, :, 5] \
                                   - 1.0/60.0*data[component_idx, :, 6])/dx
                
                diff_data[:, -3] = (1.0/60.0*data[component_idx, :, -7] - 2.0/15.0*data[component_idx, :, -6] \
                                    + 1.0/2.0*data[component_idx, :, -5] - 4.0/3.0*data[component_idx, :, -4] \
                                    + 7.0/12.0*data[component_idx, :, -3] + 2.0/5.0*data[component_idx, :, -2] \
                                    - 1.0/30.0*data[component_idx, :, -1])/dx
                
                diff_data[:, -2] = (-1.0/30.0*data[component_idx, :, -7] + 1.0/4.0*data[component_idx, :, -6] \
                                    - 5.0/6.0*data[component_idx, :, -5] + 5.0/3.0*data[component_idx, :, -4] \
                                    - 5.0/2.0*data[component_idx, :, -3] + 77.0/60.0*data[component_idx, :, -2] \
                                    + 1.0/6.0*data[component_idx, :, -1])/dx
                
                diff_data[:, -1] = (1.0/6.0*data[component_idx, :, -7] - 6.0/5.0*data[component_idx, :, -6] \
                                    + 15.0/4.0*data[component_idx, :, -5] - 20.0/3.0*data[component_idx, :, -4] \
                                    + 15.0/2.0*data[component_idx, :, -3] - 6.0*data[component_idx, :, -2] \
                                    + 49.0/20.0*data[component_idx, :, -1])/dx
            
            elif diff_data.ndim == 3:
                diff_data[:, 0, :] = (-49.0/20.0*data[component_idx, :, 0, :] + 6.0*data[component_idx, :, 1, :] \
                                      - 15.0/2.0*data[component_idx, :, 2, :] + 20.0/3.0*data[component_idx, :, 3, :] \
                                      - 15.0/4.0*data[component_idx, :, 4, :] + 6.0/5.0*data[component_idx, :, 5, :] \
                                      - 1.0/6.0*data[component_idx, :, 6, :])/dx
                
                diff_data[:, 1, :] = (-1.0/6.0*data[component_idx, :, 0, :] - 77.0/60.0*data[component_idx, :, 1, :] \
                                      + 5.0/2.0*data[component_idx, :, 2, :] - 5.0/3.0*data[component_idx, :, 3, :] \
                                      + 5.0/6.0*data[component_idx, :, 4, :] - 1.0/4.0*data[component_idx, :, 5, :] \
                                      + 1.0/30.0*data[component_idx, :, 6, :])/dx
                
                diff_data[:, 2, :] = (1.0/30.0*data[component_idx, :, 0, :] - 2.0/5.0*data[component_idx, :, 1, :] \
                                      - 7.0/12.0*data[component_idx, :, 2, :] + 4.0/3.0*data[component_idx, :, 3, :] \
                                      - 1.0/2.0*data[component_idx, :, 4, :] + 2.0/15.0*data[component_idx, :, 5, :] \
                                      - 1.0/60.0*data[component_idx, :, 6, :])/dx
                
                diff_data[:, -3, :] = (1.0/60.0*data[component_idx, :, -7, :] - 2.0/15.0*data[component_idx, :, -6, :] \
                                       + 1.0/2.0*data[component_idx, :, -5, :] - 4.0/3.0*data[component_idx, :, -4, :] \
                                       + 7.0/12.0*data[component_idx, :, -3, :] + 2.0/5.0*data[component_idx, :, -2, :] \
                                       - 1.0/30.0*data[component_idx, :, -1, :])/dx
                
                diff_data[:, -2, :] = (-1.0/30.0*data[component_idx, :, -7, :] + 1.0/4.0*data[component_idx, :, -6, :] \
                                       - 5.0/6.0*data[component_idx, :, -5, :] + 5.0/3.0*data[component_idx, :, -4, :] \
                                       - 5.0/2.0*data[component_idx, :, -3, :] + 77.0/60.0*data[component_idx, :, -2, :] \
                                       + 1.0/6.0*data[component_idx, :, -1, :])/dx
                
                diff_data[:, -1, :] = (1.0/6.0*data[component_idx, :, -7, :] - 6.0/5.0*data[component_idx, :, -6, :] \
                                       + 15.0/4.0*data[component_idx, :, -5, :] - 20.0/3.0*data[component_idx, :, -4, :] \
                                       + 15.0/2.0*data[component_idx, :, -3, :] - 6.0*data[component_idx, :, -2, :] \
                                       + 49.0/20.0*data[component_idx, :, -1, :])/dx
            
            else:
                raise RuntimeError('Data dimension > 3 not supported!')
        
        elif direction == 2:
            if diff_data.ndim < 3:
                raise IOError('There is no third direction in data with less than three dimensions!')
            
            elif diff_data.ndim == 3:
                diff_data[:, :, 0] = (-49.0/20.0*data[component_idx, :, :, 0] + 6.0*data[component_idx, :, :, 1] \
                                      - 15.0/2.0*data[component_idx, :, :, 2] + 20.0/3.0*data[component_idx, :, :, 3] \
                                      - 15.0/4.0*data[component_idx, :, :, 4] + 6.0/5.0*data[component_idx, :, :, 5] \
                                      - 1.0/6.0*data[component_idx, :, :, 6])/dx
                
                diff_data[:, :, 1] = (-1.0/6.0*data[component_idx, :, :, 0] - 77.0/60.0*data[component_idx, :, :, 1] \
                                      + 5.0/2.0*data[component_idx, :, :, 2] - 5.0/3.0*data[component_idx, :, :, 3] \
                                      + 5.0/6.0*data[component_idx, :, :, 4] - 1.0/4.0*data[component_idx, :, :, 5] \
                                      + 1.0/30.0*data[component_idx, :, :, 6])/dx
                
                diff_data[:, :, 2] = (1.0/30.0*data[component_idx, :, :, 0] - 2.0/5.0*data[component_idx, :, :, 1] \
                                      - 7.0/12.0*data[component_idx, :, :, 2] + 4.0/3.0*data[component_idx, :, :, 3] \
                                      - 1.0/2.0*data[component_idx, :, :, 4] + 2.0/15.0*data[component_idx, :, :, 5] \
                                      - 1.0/60.0*data[component_idx, :, :, 6])/dx
                
                diff_data[:, :, -3] = (1.0/60.0*data[component_idx, :, :, -7] - 2.0/15.0*data[component_idx, :, :, -6] \
                                       + 1.0/2.0*data[component_idx, :, :, -5] - 4.0/3.0*data[component_idx, :, :, -4] \
                                       + 7.0/12.0*data[component_idx, :, :, -3] + 2.0/5.0*data[component_idx, :, :, -2] \
                                       - 1.0/30.0*data[component_idx, :, :, -1])/dx
                
                diff_data[:, :, -2] = (-1.0/30.0*data[component_idx, :, :, -7] + 1.0/4.0*data[component_idx, :, :, -6] \
                                       - 5.0/6.0*data[component_idx, :, :, -5] + 5.0/3.0*data[component_idx, :, :, -4] \
                                       - 5.0/2.0*data[component_idx, :, :, -3] + 77.0/60.0*data[component_idx, :, :, -2] \
                                       + 1.0/6.0*data[component_idx, :, :, -1])/dx
                
                diff_data[:, :, -1] = (1.0/6.0*data[component_idx, :, :, -7] - 6.0/5.0*data[component_idx, :, :, -6] \
                                       + 15.0/4.0*data[component_idx, :, :, -5] - 20.0/3.0*data[component_idx, :, :, -4] \
                                       + 15.0/2.0*data[component_idx, :, :, -3] - 6.0*data[component_idx, :, :, -2] \
                                       + 49.0/20.0*data[component_idx, :, :, -1])/dx
            
            else:
                raise RuntimeError('Data dimension > 3 not supported!')   
    
    return diff_data

