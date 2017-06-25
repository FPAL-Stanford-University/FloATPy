"""
Module for upsampling data.
"""

import numpy

def getNumberOfGhostCellsUpsampling(method = 'constant'):
    """
    Determine the number of ghost cells needed for upsampling in the interior of domain.
    """
    
    if method == 'constant':
        return 0
    elif method == 'second_order_Lagrange':
        return 1
    elif method == 'fourth_order_Lagrange':
        return 2
    elif method == 'sixth_order_Lagrange':
        return 3
    else:
        raise RuntimeError("Unknown method '" + method + "' for number of ghost cells!")


def upsample(data, refine_ratio, component_idx = 0, method = 'constant'):
    """
    Upsampling the data.
    """
    
    if method == 'constant':
        return upsampleConstant(data, refine_ratio, component_idx)
    elif method == 'second_order_Lagrange':
        return upsampleSecondOrderLagrange(data, refine_ratio, component_idx)
    elif method == 'fourth_order_Lagrange':
        return upsampleFourthOrderLagrange(data, refine_ratio, component_idx)
    elif method == 'sixth_order_Lagrange':
        return upsampleSixthOrderLagrange(data, refine_ratio, component_idx)
    else:
        raise RuntimeError("Unknown method '" + method + "' for upsampling!")


def upsampleConstant(data, refine_ratio, component_idx = 0):
    """
    Upsampling the data using constant interpolation.
    """
    
    r = refine_ratio
    
    # Get the order and shape of data.
    
    data_order = 'C'
    if numpy.isfortran(data):
        data_order = 'F'
    
    # Check whether the dimension of data is valid.
    
    if data.ndim < 2 or data.ndim > 4:
        raise RuntimeError('Shape of data is invalid!')
    
    # Check whether the component_idx is valid and get the shape of the component's data.
    
    if data_order == 'C':
        if component_idx >= data.shape[0] or component_idx < 0:
            raise RuntimeError('Component index is invalid!')
        
        data_shape = numpy.array(data.shape[1:])
    
    else:
        if component_idx >= data.shape[-1] or component_idx < 0:
            raise RuntimeError('Component index is invalid!')
        
        data_shape = numpy.array(data.shape[:-1])
    
    # Get the dimension of data.
    
    dim = data_shape.shape[0]
    
    # Initialize container to store the upsampled data.
    
    upsampled_data = None
    
    # Get the component's data.
    
    data_component = None
    
    if data_order == 'C':
        if dim == 1:
            data_component = data[component_idx, :]
        elif dim == 2:
            data_component = data[component_idx, :, :]
        elif dim == 3:
            data_component = data[component_idx, :, :, :]
    else:
        if dim == 1:
            data_component = data[:, component_idx]
        elif dim == 2:
            data_component = data[:, :, component_idx]
        elif dim == 3:
            data_component = data[:, :, :, component_idx]
    
    # Upsample the data with constant interpolation.
    
    if dim == 1:
        upsampled_data = numpy.repeat(data_component, r[0], axis = 0)
    
    elif dim == 2:
        if data_order == 'C':
            upsampled_data = numpy.repeat(data_component, r[0], axis = 0)
            upsampled_data = numpy.repeat(upsampled_data, r[1], axis = 1)
        
        else:
            upsampled_data = numpy.ravel(data_component, order = 'F')
            upsampled_data = upsampled_data.reshape(
                (data_shape[1], data_shape[0]), order = 'C')
            
            upsampled_data = numpy.repeat(upsampled_data, r[1], axis = 0)
            upsampled_data = numpy.repeat(upsampled_data, r[0], axis = 1)
            
            upsampled_data = numpy.ravel(upsampled_data, order = 'C')
            upsampled_data = upsampled_data.reshape(
                data_shape*r[0:dim], order = 'F')
    
    elif dim == 3:
        if data_order == 'C':
            upsampled_data = numpy.repeat(data_component, r[0], axis = 0)
            upsampled_data = numpy.repeat(upsampled_data, r[1], axis = 1)
            upsampled_data = numpy.repeat(upsampled_data, r[2], axis = 2)
        
        else:
            upsampled_data = numpy.ravel(data_component, order = 'F')
            upsampled_data = upsampled_data.reshape(
                (data_shape[2], data_shape[1], data_shape[0]), order = 'C')
            
            upsampled_data = numpy.repeat(upsampled_data, r[2], axis = 0)
            upsampled_data = numpy.repeat(upsampled_data, r[1], axis = 1)
            upsampled_data = numpy.repeat(upsampled_data, r[0], axis = 2)
            
            upsampled_data = numpy.ravel(upsampled_data, order = 'C')
            upsampled_data = upsampled_data.reshape(
                data_shape*r[0:dim], order = 'F')
    
    return upsampled_data


def upsampleSecondOrderLagrange(data, refine_ratio, component_idx = 0):
    """
    Upsampling the data using second order Lagrange interpolation.
    """
    
    r = refine_ratio
    
    stencil_size = 2
    half_stencil_size = stencil_size/2
    
    # Get the order and shape of data.
    
    data_order = 'C'
    if numpy.isfortran(data):
        data_order = 'F'
    
    # Check whether the dimension of data is valid.
    
    if data.ndim < 2 or data.ndim > 4:
        raise RuntimeError('Shape of data is invalid!')
    
    # Check whether the component_idx is valid and get the shape of the component's data.
    
    if data_order == 'C':
        if component_idx >= data.shape[0] or component_idx < 0:
            raise RuntimeError('Component index is invalid!')
        
        data_shape = numpy.array(data.shape[1:])
    
    else:
        if component_idx >= data.shape[-1] or component_idx < 0:
            raise RuntimeError('Component index is invalid!')
        
        data_shape = numpy.array(data.shape[:-1])
    
    # Get the dimension of data.
    
    dim = data_shape.shape[0]
    
    # Check whether data size is large enough for second order Lagrange interpolation.
    
    if data_shape[0] < 3:
        raise RuntimeError('First dimension of data is not large enough!')
    
    if dim > 1:
        if data_shape[1] < 3:
            raise RuntimeError('Second dimension of data is not large enough!')
        
    elif dim > 2:
        if data_shape[2] < 3:
            raise RuntimeError('Third dimension of data is not large enough!')
    
    # Compute the coefficients for second order Lagrange interpolation.
    
    c_0 = None
    c_1 = None
    c_2 = None
    
    c_0 = numpy.ones([r[0], stencil_size], dtype=data.dtype)
    if dim == 2:
        c_1 = numpy.ones([r[1], stencil_size], dtype=data.dtype)
    if dim == 3:
        c_2 = numpy.ones([r[2], stencil_size], dtype=data.dtype)
    
    if r[0] % 2 == 0: # If the upsampled nodes don't overlap with the original nodes.
        
        delta = 1.0/r[0]
        for s in range(r[0]):
            idx_r = (half_stencil_size - 1) + (s + 0.5)*delta
            for i in range(stencil_size):
                for j in range(stencil_size):
                    if (i != j):
                        c_0[s, i] = c_0[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 2:
            delta = 1.0/r[1]
            for s in range(r[1]):
                idx_r = (half_stencil_size - 1) + (s + 0.5)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_1[s, i] = c_1[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 3:
            delta = 1.0/r[2]
            for s in range(r[2]):
                idx_r = (half_stencil_size - 1) + (s + 0.5)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_2[s, i] = c_2[s, i] * (idx_r - j)/(i - j)
    
    else: # If some of the upsampled nodes overlap with the original nodes.
        
        delta = 1.0/r[0]
        for s in range(r[0]):
            idx_r = (half_stencil_size - 1) + (s + 1.0)*delta
            for i in range(stencil_size):
                for j in range(stencil_size):
                    if (i != j):
                        c_0[s, i] = c_0[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 2:
            delta = 1.0/r[1]
            for s in range(r[1]):
                idx_r = (half_stencil_size - 1) + (s + 1.0)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_1[s, i] = c_1[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 3:
            delta = 1.0/r[2]
            for s in range(r[2]):
                idx_r = (half_stencil_size - 1) + (s + 1.0)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_2[s, i] = c_2[s, i] * (idx_r - j)/(i - j)
    
    # Initialize containers to store the upsampled data. The elements in the containers
    # are initialized as NAN values.
    
    upsampled_data_shape_0 = None
    upsampled_data_shape_1 = None
    upsampled_data_shape = None
    
    upsampled_data_0 = None
    upsampled_data_1 = None
    upsampled_data = None
    
    if dim == 1:
        upsampled_data_shape = numpy.multiply(data_shape, r[0])
        upsampled_data = numpy.empty(upsampled_data_shape, dtype = data.dtype, order = data_order)
        upsampled_data[:] = numpy.NAN
    
    elif dim == 2:
        upsampled_data_shape_0 = numpy.copy(data_shape)
        upsampled_data_shape_0[0] = upsampled_data_shape_0[0]*r[0]
        
        upsampled_data_shape = numpy.copy(upsampled_data_shape_0)
        upsampled_data_shape[1] = upsampled_data_shape[1]*r[1]
        
        upsampled_data_0 = numpy.empty(upsampled_data_shape_0, dtype = data.dtype, order = data_order)
        upsampled_data = numpy.empty(upsampled_data_shape, dtype = data.dtype, order = data_order)
        
        upsampled_data_0[:] = numpy.NAN
        upsampled_data[:] = numpy.NAN
        
    elif dim == 3:
        upsampled_data_shape_0 = numpy.copy(data_shape)
        upsampled_data_shape_0[0] = upsampled_data_shape_0[0]*r[0]
        
        upsampled_data_shape_1 = numpy.copy(upsampled_data_shape_0)
        upsampled_data_shape_1[1] = upsampled_data_shape_1[1]*r[1]
        
        upsampled_data_shape = numpy.copy(upsampled_data_shape_1)
        upsampled_data_shape[2] = upsampled_data_shape[2]*r[2]
        
        upsampled_data_0 = numpy.empty(upsampled_data_shape_0, dtype = data.dtype, order = data_order)
        upsampled_data_1 = numpy.empty(upsampled_data_shape_1, dtype = data.dtype, order = data_order)
        upsampled_data = numpy.empty(upsampled_data_shape, dtype = data.dtype, order = data_order)
        
        upsampled_data_0[:] = numpy.NAN
        upsampled_data_1[:] = numpy.NAN
        upsampled_data[:] = numpy.NAN
    
    # Get the component's data.
    
    data_component = None
    
    if data_order == 'C':
        if dim == 1:
            data_component = data[component_idx, :]
        elif dim == 2:
            data_component = data[component_idx, :, :]
        elif dim == 3:
            data_component = data[component_idx, :, :, :]
    else:
        if dim == 1:
            data_component = data[:, component_idx]
        elif dim == 2:
            data_component = data[:, :, component_idx]
        elif dim == 3:
            data_component = data[:, :, :, component_idx]
    
    # Upsample the data with second order Lagrange interpolation.
    
    if dim == 1:
        start_idx_fine = (half_stencil_size - 1)*r[0] + (r[0] + 1)/2
        end_idx_fine = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
        
        for s in range(r[0]):
            upsampled_data[(start_idx_fine + s):end_idx_fine:r[0]] = c_0[s, 0]*data_component[0:(-stencil_size + 1)] \
                + c_0[s, 1]*data_component[1:]
    
    elif dim == 2:
        start_idx_fine_0 = (half_stencil_size - 1)*r[0] + (r[0] + 1)/2
        end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
        
        start_idx_fine_1 = (half_stencil_size - 1)*r[1] + (r[1] + 1)/2
        end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
        
        for s in range(r[0]):
            upsampled_data_0[(start_idx_fine_0 + s):end_idx_fine_0:r[0], :] = \
                c_0[s, 0]*data_component[0:(-stencil_size + 1), :] + c_0[s, 1]*data_component[1:, :]
        
        for s in range(r[1]):
            upsampled_data[start_idx_fine_0:end_idx_fine_0, (start_idx_fine_1 + s):end_idx_fine_1:r[1]] = \
                c_1[s, 0]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 0:(-stencil_size + 1)] \
                + c_1[s, 1]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 1:]
    
    elif dim == 3:
        start_idx_fine_0 = (half_stencil_size - 1)*r[0] + (r[0] + 1)/2
        end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
        
        start_idx_fine_1 = (half_stencil_size - 1)*r[1] + (r[1] + 1)/2
        end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
        
        start_idx_fine_2 = (half_stencil_size - 1)*r[2] + (r[2] + 1)/2
        end_idx_fine_2 = upsampled_data_shape[2] - (half_stencil_size - 1)*r[2] - r[2]/2
        
        for s in range(r[0]):
            upsampled_data_0[(start_idx_fine_0 + s):end_idx_fine_0:r[0], : :] = \
                c_0[s, 0]*data_component[0:(-stencil_size + 1), :, :] \
                + c_0[s, 1]*data_component[1:, :, :]
        
        for s in range(r[1]):
            upsampled_data_1[start_idx_fine_0:end_idx_fine_0, (start_idx_fine_1 + s):end_idx_fine_1:r[1], :] = \
                c_1[s, 0]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 0:(-stencil_size + 1), :] \
                + c_1[s, 1]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 1:, :]
        
        for s in range(r[2]):
            upsampled_data[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, (start_idx_fine_2 + s):end_idx_fine_2:r[2]] = \
                c_2[s, 0]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 0:(-stencil_size + 1)] \
                + c_2[s, 1]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 1:]
    
    # If some of the upsampled nodes overlap with the original nodes, copy the data from the original nodes to the first
    # sets of overlapped nodes since they are not upsampled yet.
    
    if r[0] % 2 == 1: # If some of the upsampled nodes overlap with the original nodes.
        overlap_idx_coarse = half_stencil_size - 1
        overlap_idx_fine = (half_stencil_size - 1)*r[0] + r[0]/2
        
        if dim == 1:
            upsampled_data[overlap_idx_fine] = data_component[overlap_idx_coarse]
        
        elif dim == 2:
            start_idx_fine_0 = (half_stencil_size - 1)*r[0] + r[0]/2
            end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
            
            start_idx_fine_1 = (half_stencil_size - 1)*r[1] + r[1]/2
            end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
            
            upsampled_line_data_0 = numpy.repeat(data_component[overlap_idx_coarse, :], \
                r[1])
            
            upsampled_line_data_1 = numpy.repeat(data_component[:, overlap_idx_coarse], \
                r[0])
            
            upsampled_data[overlap_idx_fine, start_idx_fine_1:end_idx_fine_1] = \
                upsampled_line_data_0[start_idx_fine_1:end_idx_fine_1]
            
            upsampled_data[start_idx_fine_0:end_idx_fine_0, overlap_idx_fine] = \
                upsampled_line_data_1[start_idx_fine_0:end_idx_fine_0]
        
        elif dim == 3:
            start_idx_fine_0 = (half_stencil_size - 1)*r[0] + r[0]/2
            end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
            
            start_idx_fine_1 = (half_stencil_size - 1)*r[1] + r[1]/2
            end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
            
            start_idx_fine_2 = (half_stencil_size - 1)*r[2] + r[2]/2
            end_idx_fine_2 = upsampled_data_shape[2] - (half_stencil_size - 1)*r[2] - r[2]/2
            
            upsampled_face_data_0 = numpy.repeat(data_component[overlap_idx_coarse, :, :], \
                r[1], axis = 0)
            upsampled_face_data_0 = numpy.repeat(upsampled_face_data_0, r[2], axis = 1)
            
            upsampled_face_data_1 = numpy.repeat(data_component[:, overlap_idx_coarse, :], \
                r[0], axis = 0)
            upsampled_face_data_1 = numpy.repeat(upsampled_face_data_1, r[2], axis = 1)
            
            upsampled_face_data_2 = numpy.repeat(data_component[:, :, overlap_idx_coarse], \
                r[0], axis = 0)
            upsampled_face_data_2 = numpy.repeat(upsampled_face_data_2, r[1], axis = 1)
            
            upsampled_data[overlap_idx_fine, start_idx_fine_1:end_idx_fine_1, start_idx_fine_2:end_idx_fine_2] = \
                upsampled_face_data_0[start_idx_fine_1:end_idx_fine_1, start_idx_fine_2:end_idx_fine_2]
            
            upsampled_data[start_idx_fine_0:end_idx_fine_0, overlap_idx_fine, start_idx_fine_2:end_idx_fine_2] = \
                upsampled_face_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_2:end_idx_fine_2]
            
            upsampled_data[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, overlap_idx_fine] = \
                upsampled_face_data_2[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1]
    
    return upsampled_data


def upsampleFourthOrderLagrange(data, refine_ratio, component_idx = 0):
    """
    Upsampling the data using fourth order Lagrange interpolation.
    """
    
    r = refine_ratio
    
    stencil_size = 4
    half_stencil_size = stencil_size/2
    
    # Get the order and shape of data.
    
    data_order = 'C'
    if numpy.isfortran(data):
        data_order = 'F'
    
    # Check whether the dimension of data is valid.
    
    if data.ndim < 2 or data.ndim > 4:
        raise RuntimeError('Shape of data is invalid!')
    
    # Check whether the component_idx is valid and get the shape of the component's data.
    
    if data_order == 'C':
        if component_idx >= data.shape[0] or component_idx < 0:
            raise RuntimeError('Component index is invalid!')
        
        data_shape = numpy.array(data.shape[1:])
    
    else:
        if component_idx >= data.shape[-1] or component_idx < 0:
            raise RuntimeError('Component index is invalid!')
        
        data_shape = numpy.array(data.shape[:-1])
    
    # Get the dimension of data.
    
    dim = data_shape.shape[0]
    
    # Check whether data size is large enough for fourth order Lagrange interpolation.
    
    if data_shape[0] < 5:
        raise RuntimeError('First dimension of data is not large enough!')
    
    if dim > 1:
        if data_shape[1] < 5:
            raise RuntimeError('Second dimension of data is not large enough!')
        
    elif dim > 2:
        if data_shape[2] < 5:
            raise RuntimeError('Third dimension of data is not large enough!')
    
    # Compute the coefficients for fourth order Lagrange interpolation.
    
    c_0 = None
    c_1 = None
    c_2 = None
    
    c_0 = numpy.ones([r[0], stencil_size], dtype=data.dtype)
    if dim == 2:
        c_1 = numpy.ones([r[1], stencil_size], dtype=data.dtype)
    if dim == 3:
        c_2 = numpy.ones([r[2], stencil_size], dtype=data.dtype)
    
    if r[0] % 2 == 0: # If the upsampled nodes don't overlap with the original nodes.
        
        delta = 1.0/r[0]
        for s in range(r[0]):
            idx_r = (half_stencil_size - 1) + (s + 0.5)*delta
            for i in range(stencil_size):
                for j in range(stencil_size):
                    if (i != j):
                        c_0[s, i] = c_0[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 2:
            delta = 1.0/r[1]
            for s in range(r[1]):
                idx_r = (half_stencil_size - 1) + (s + 0.5)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_1[s, i] = c_1[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 3:
            delta = 1.0/r[2]
            for s in range(r[2]):
                idx_r = (half_stencil_size - 1) + (s + 0.5)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_2[s, i] = c_2[s, i] * (idx_r - j)/(i - j)
    
    else: # If some of the upsampled nodes overlap with the original nodes.
        
        delta = 1.0/r[0]
        for s in range(r[0]):
            idx_r = (half_stencil_size - 1) + (s + 1.0)*delta
            for i in range(stencil_size):
                for j in range(stencil_size):
                    if (i != j):
                        c_0[s, i] = c_0[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 2:
            delta = 1.0/r[1]
            for s in range(r[1]):
                idx_r = (half_stencil_size - 1) + (s + 1.0)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_1[s, i] = c_1[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 3:
            delta = 1.0/r[2]
            for s in range(r[2]):
                idx_r = (half_stencil_size - 1) + (s + 1.0)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_2[s, i] = c_2[s, i] * (idx_r - j)/(i - j)
        
    # Initialize containers to store the upsampled data. The elements in the containers
    # are initialized as NAN values.
    
    upsampled_data_shape_0 = None
    upsampled_data_shape_1 = None
    upsampled_data_shape = None
    
    upsampled_data_0 = None
    upsampled_data_1 = None
    upsampled_data = None
    
    if dim == 1:
        upsampled_data_shape = numpy.multiply(data_shape, r[0])
        upsampled_data = numpy.empty(upsampled_data_shape, dtype = data.dtype, order = data_order)
        upsampled_data[:] = numpy.NAN
    
    elif dim == 2:
        upsampled_data_shape_0 = numpy.copy(data_shape)
        upsampled_data_shape_0[0] = upsampled_data_shape_0[0]*r[0]
        
        upsampled_data_shape = numpy.copy(upsampled_data_shape_0)
        upsampled_data_shape[1] = upsampled_data_shape[1]*r[1]
        
        upsampled_data_0 = numpy.empty(upsampled_data_shape_0, dtype = data.dtype, order = data_order)
        upsampled_data = numpy.empty(upsampled_data_shape, dtype = data.dtype, order = data_order)
        
        upsampled_data_0[:] = numpy.NAN
        upsampled_data[:] = numpy.NAN
        
    elif dim == 3:
        upsampled_data_shape_0 = numpy.copy(data_shape)
        upsampled_data_shape_0[0] = upsampled_data_shape_0[0]*r[0]
        
        upsampled_data_shape_1 = numpy.copy(upsampled_data_shape_0)
        upsampled_data_shape_1[1] = upsampled_data_shape_1[1]*r[1]
        
        upsampled_data_shape = numpy.copy(upsampled_data_shape_1)
        upsampled_data_shape[2] = upsampled_data_shape[2]*r[2]
        
        upsampled_data_0 = numpy.empty(upsampled_data_shape_0, dtype = data.dtype, order = data_order)
        upsampled_data_1 = numpy.empty(upsampled_data_shape_1, dtype = data.dtype, order = data_order)
        upsampled_data = numpy.empty(upsampled_data_shape, dtype = data.dtype, order = data_order)
        
        upsampled_data_0[:] = numpy.NAN
        upsampled_data_1[:] = numpy.NAN
        upsampled_data[:] = numpy.NAN
    
    # Get the component's data.
    
    data_component = None
    
    if data_order == 'C':
        if dim == 1:
            data_component = data[component_idx, :]
        elif dim == 2:
            data_component = data[component_idx, :, :]
        elif dim == 3:
            data_component = data[component_idx, :, :, :]
    else:
        if dim == 1:
            data_component = data[:, component_idx]
        elif dim == 2:
            data_component = data[:, :, component_idx]
        elif dim == 3:
            data_component = data[:, :, :, component_idx]
    
    # Upsample the data with fourth order Lagrange interpolation.
    
    if dim == 1:
        start_idx_fine = (half_stencil_size - 1)*r[0] + (r[0] + 1)/2
        end_idx_fine = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
        
        for s in range(r[0]):
            upsampled_data[(start_idx_fine + s):end_idx_fine:r[0]] = c_0[s, 0]*data_component[0:(-stencil_size + 1)] \
                + c_0[s, 1]*data_component[1:(-stencil_size + 2)] + c_0[s, 2]*data_component[2:(-stencil_size + 3)] \
                + c_0[s, 3]*data_component[3:]
    
    elif dim == 2:
        start_idx_fine_0 = (half_stencil_size - 1)*r[0] + (r[0] + 1)/2
        end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
        
        start_idx_fine_1 = (half_stencil_size - 1)*r[1] + (r[1] + 1)/2
        end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
        
        for s in range(r[0]):
            upsampled_data_0[(start_idx_fine_0 + s):end_idx_fine_0:r[0], :] = \
                c_0[s, 0]*data_component[0:(-stencil_size + 1), :] + c_0[s, 1]*data_component[1:(-stencil_size + 2), :] \
                + c_0[s, 2]*data_component[2:(-stencil_size + 3), :] + c_0[s, 3]*data_component[3:, :]
        
        for s in range(r[1]):
            upsampled_data[start_idx_fine_0:end_idx_fine_0, (start_idx_fine_1 + s):end_idx_fine_1:r[1]] = \
                c_1[s, 0]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 0:(-stencil_size + 1)] \
                + c_1[s, 1]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 1:(-stencil_size + 2)] \
                + c_1[s, 2]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 2:(-stencil_size + 3)] \
                + c_1[s, 3]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 3:]
    
    elif dim == 3:
        start_idx_fine_0 = (half_stencil_size - 1)*r[0] + (r[0] + 1)/2
        end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
        
        start_idx_fine_1 = (half_stencil_size - 1)*r[1] + (r[1] + 1)/2
        end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
        
        start_idx_fine_2 = (half_stencil_size - 1)*r[2] + (r[2] + 1)/2
        end_idx_fine_2 = upsampled_data_shape[2] - (half_stencil_size - 1)*r[2] - r[2]/2
        
        for s in range(r[0]):
            upsampled_data_0[(start_idx_fine_0 + s):end_idx_fine_0:r[0], : :] = \
                c_0[s, 0]*data_component[0:(-stencil_size + 1), :, :] \
                + c_0[s, 1]*data_component[1:(-stencil_size + 2), :, :] \
                + c_0[s, 2]*data_component[2:(-stencil_size + 3), :, :] \
                + c_0[s, 3]*data_component[3:, :, :]
        
        for s in range(r[1]):
            upsampled_data_1[start_idx_fine_0:end_idx_fine_0, (start_idx_fine_1 + s):end_idx_fine_1:r[1], :] = \
                c_1[s, 0]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 0:(-stencil_size + 1), :] \
                + c_1[s, 1]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 1:(-stencil_size + 2), :] \
                + c_1[s, 2]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 2:(-stencil_size + 3), :] \
                + c_1[s, 3]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 3:, :]
        
        for s in range(r[2]):
            upsampled_data[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, (start_idx_fine_2 + s):end_idx_fine_2:r[2]] = \
                c_2[s, 0]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 0:(-stencil_size + 1)] \
                + c_2[s, 1]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 1:(-stencil_size + 2)] \
                + c_2[s, 2]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 2:(-stencil_size + 3)] \
                + c_2[s, 3]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 3:]
    
    # If some of the upsampled nodes overlap with the original nodes, copy the data from the original nodes to the first
    # sets of overlapped nodes since they are not upsampled yet.
    
    if r[0] % 2 == 1: # If some of the upsampled nodes overlap with the original nodes.
        overlap_idx_coarse = half_stencil_size - 1
        overlap_idx_fine = (half_stencil_size - 1)*r[0] + r[0]/2
        
        if dim == 1:
            upsampled_data[overlap_idx_fine] = data_component[overlap_idx_coarse]
        
        elif dim == 2:
            start_idx_fine_0 = (half_stencil_size - 1)*r[0] + r[0]/2
            end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
            
            start_idx_fine_1 = (half_stencil_size - 1)*r[1] + r[1]/2
            end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
            
            upsampled_line_data_0 = numpy.repeat(data_component[overlap_idx_coarse, :], \
                r[1])
            
            upsampled_line_data_1 = numpy.repeat(data_component[:, overlap_idx_coarse], \
                r[0])
            
            upsampled_data[overlap_idx_fine, start_idx_fine_1:end_idx_fine_1] = \
                upsampled_line_data_0[start_idx_fine_1:end_idx_fine_1]
            
            upsampled_data[start_idx_fine_0:end_idx_fine_0, overlap_idx_fine] = \
                upsampled_line_data_1[start_idx_fine_0:end_idx_fine_0]
        
        elif dim == 3:
            start_idx_fine_0 = (half_stencil_size - 1)*r[0] + r[0]/2
            end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
            
            start_idx_fine_1 = (half_stencil_size - 1)*r[1] + r[1]/2
            end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
            
            start_idx_fine_2 = (half_stencil_size - 1)*r[2] + r[2]/2
            end_idx_fine_2 = upsampled_data_shape[2] - (half_stencil_size - 1)*r[2] - r[2]/2
            
            upsampled_face_data_0 = numpy.repeat(data_component[overlap_idx_coarse, :, :], \
                r[1], axis = 0)
            upsampled_face_data_0 = numpy.repeat(upsampled_face_data_0, r[2], axis = 1)
            
            upsampled_face_data_1 = numpy.repeat(data_component[:, overlap_idx_coarse, :], \
                r[0], axis = 0)
            upsampled_face_data_1 = numpy.repeat(upsampled_face_data_1, r[2], axis = 1)
            
            upsampled_face_data_2 = numpy.repeat(data_component[:, :, overlap_idx_coarse], \
                r[0], axis = 0)
            upsampled_face_data_2 = numpy.repeat(upsampled_face_data_2, r[1], axis = 1)
            
            upsampled_data[overlap_idx_fine, start_idx_fine_1:end_idx_fine_1, start_idx_fine_2:end_idx_fine_2] = \
                upsampled_face_data_0[start_idx_fine_1:end_idx_fine_1, start_idx_fine_2:end_idx_fine_2]
            
            upsampled_data[start_idx_fine_0:end_idx_fine_0, overlap_idx_fine, start_idx_fine_2:end_idx_fine_2] = \
                upsampled_face_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_2:end_idx_fine_2]
            
            upsampled_data[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, overlap_idx_fine] = \
                upsampled_face_data_2[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1]
    
    return upsampled_data


def upsampleSixthOrderLagrange(data, refine_ratio, component_idx = 0):
    """
    Upsampling the data using sixth order Lagrange interpolation.
    """
    
    r = refine_ratio
    
    stencil_size = 6
    half_stencil_size = stencil_size/2
    
    # Get the order and shape of data.
    
    data_order = 'C'
    if numpy.isfortran(data):
        data_order = 'F'
    
    # Check whether the dimension of data is valid.
    
    if data.ndim < 2 or data.ndim > 4:
        raise RuntimeError('Shape of data is invalid!')
    
    # Check whether the component_idx is valid and get the shape of the component's data.
    
    if data_order == 'C':
        if component_idx >= data.shape[0] or component_idx < 0:
            raise RuntimeError('Component index is invalid!')
        
        data_shape = numpy.array(data.shape[1:])
    
    else:
        if component_idx >= data.shape[-1] or component_idx < 0:
            raise RuntimeError('Component index is invalid!')
        
        data_shape = numpy.array(data.shape[:-1])
    
    # Get the dimension of data.
    
    dim = data_shape.shape[0]
    
    # Check whether data size is large enough for sixth order Lagrange interpolation.
    
    if data_shape[0] < 7:
        raise RuntimeError('First dimension of data is not large enough!')
    
    if dim > 1:
        if data_shape[1] < 7:
            raise RuntimeError('Second dimension of data is not large enough!')
        
    elif dim > 2:
        if data_shape[2] < 7:
            raise RuntimeError('Third dimension of data is not large enough!')
    
    # Compute the coefficients for sixth order Lagrange interpolation.
    
    c_0 = None
    c_1 = None
    c_2 = None
    
    c_0 = numpy.ones([r[0], stencil_size], dtype=data.dtype)
    if dim == 2:
        c_1 = numpy.ones([r[1], stencil_size], dtype=data.dtype)
    if dim == 3:
        c_2 = numpy.ones([r[2], stencil_size], dtype=data.dtype)
    
    if r[0] % 2 == 0: # If the upsampled nodes don't overlap with the original nodes.
        
        delta = 1.0/r[0]
        for s in range(r[0]):
            idx_r = (half_stencil_size - 1) + (s + 0.5)*delta
            for i in range(stencil_size):
                for j in range(stencil_size):
                    if (i != j):
                        c_0[s, i] = c_0[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 2:
            delta = 1.0/r[1]
            for s in range(r[1]):
                idx_r = (half_stencil_size - 1) + (s + 0.5)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_1[s, i] = c_1[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 3:
            delta = 1.0/r[2]
            for s in range(r[2]):
                idx_r = (half_stencil_size - 1) + (s + 0.5)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_2[s, i] = c_2[s, i] * (idx_r - j)/(i - j)
    
    else: # If some of the upsampled nodes overlap with the original nodes.
        
        delta = 1.0/r[0]
        for s in range(r[0]):
            idx_r = (half_stencil_size - 1) + (s + 1.0)*delta
            for i in range(stencil_size):
                for j in range(stencil_size):
                    if (i != j):
                        c_0[s, i] = c_0[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 2:
            delta = 1.0/r[1]
            for s in range(r[1]):
                idx_r = (half_stencil_size - 1) + (s + 1.0)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_1[s, i] = c_1[s, i] * (idx_r - j)/(i - j)
        
        if dim >= 3:
            delta = 1.0/r[2]
            for s in range(r[2]):
                idx_r = (half_stencil_size - 1) + (s + 1.0)*delta
                for i in range(stencil_size):
                    for j in range(stencil_size):
                        if (i != j):
                            c_2[s, i] = c_2[s, i] * (idx_r - j)/(i - j)
    
    # Initialize containers to store the upsampled data. The elements in the containers
    # are initialized as NAN values.
    
    upsampled_data_shape_0 = None
    upsampled_data_shape_1 = None
    upsampled_data_shape = None
    
    upsampled_data_0 = None
    upsampled_data_1 = None
    upsampled_data = None
    
    if dim == 1:
        upsampled_data_shape = numpy.multiply(data_shape, r[0])
        upsampled_data = numpy.empty(upsampled_data_shape, dtype = data.dtype, order = data_order)
        upsampled_data[:] = numpy.NAN
    
    elif dim == 2:
        upsampled_data_shape_0 = numpy.copy(data_shape)
        upsampled_data_shape_0[0] = upsampled_data_shape_0[0]*r[0]
        
        upsampled_data_shape = numpy.copy(upsampled_data_shape_0)
        upsampled_data_shape[1] = upsampled_data_shape[1]*r[1]
        
        upsampled_data_0 = numpy.empty(upsampled_data_shape_0, dtype = data.dtype, order = data_order)
        upsampled_data = numpy.empty(upsampled_data_shape, dtype = data.dtype, order = data_order)
        
        upsampled_data_0[:] = numpy.NAN
        upsampled_data[:] = numpy.NAN
        
    elif dim == 3:
        upsampled_data_shape_0 = numpy.copy(data_shape)
        upsampled_data_shape_0[0] = upsampled_data_shape_0[0]*r[0]
        
        upsampled_data_shape_1 = numpy.copy(upsampled_data_shape_0)
        upsampled_data_shape_1[1] = upsampled_data_shape_1[1]*r[1]
        
        upsampled_data_shape = numpy.copy(upsampled_data_shape_1)
        upsampled_data_shape[2] = upsampled_data_shape[2]*r[2]
        
        upsampled_data_0 = numpy.empty(upsampled_data_shape_0, dtype = data.dtype, order = data_order)
        upsampled_data_1 = numpy.empty(upsampled_data_shape_1, dtype = data.dtype, order = data_order)
        upsampled_data = numpy.empty(upsampled_data_shape, dtype = data.dtype, order = data_order)
        
        upsampled_data_0[:] = numpy.NAN
        upsampled_data_1[:] = numpy.NAN
        upsampled_data[:] = numpy.NAN
    
    # Get the component's data.
    
    data_component = None
    
    if data_order == 'C':
        if dim == 1:
            data_component = data[component_idx, :]
        elif dim == 2:
            data_component = data[component_idx, :, :]
        elif dim == 3:
            data_component = data[component_idx, :, :, :]
    else:
        if dim == 1:
            data_component = data[:, component_idx]
        elif dim == 2:
            data_component = data[:, :, component_idx]
        elif dim == 3:
            data_component = data[:, :, :, component_idx]
    
    # Upsample the data with sixth order Lagrange interpolation.
    
    if dim == 1:
        start_idx_fine = (half_stencil_size - 1)*r[0] + (r[0] + 1)/2
        end_idx_fine = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
        
        for s in range(r[0]):
            upsampled_data[(start_idx_fine + s):end_idx_fine:r[0]] = c_0[s, 0]*data_component[0:(-stencil_size + 1)] \
                + c_0[s, 1]*data_component[1:(-stencil_size + 2)] + c_0[s, 2]*data_component[2:(-stencil_size + 3)] \
                + c_0[s, 3]*data_component[3:(-stencil_size + 4)] + c_0[s, 4]*data_component[4:(-stencil_size + 5)] \
                + c_0[s, 5]*data_component[5:]
    
    elif dim == 2:
        start_idx_fine_0 = (half_stencil_size - 1)*r[0] + (r[0] + 1)/2
        end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
        
        start_idx_fine_1 = (half_stencil_size - 1)*r[1] + (r[1] + 1)/2
        end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
        
        for s in range(r[0]):
            upsampled_data_0[(start_idx_fine_0 + s):end_idx_fine_0:r[0], :] = \
                c_0[s, 0]*data_component[0:(-stencil_size + 1), :] + c_0[s, 1]*data_component[1:(-stencil_size + 2), :] \
                + c_0[s, 2]*data_component[2:(-stencil_size + 3), :] + c_0[s, 3]*data_component[3:(-stencil_size + 4), :] \
                + c_0[s, 4]*data_component[4:(-stencil_size + 5), :] + c_0[s, 5]*data_component[5:, :]
        
        for s in range(r[1]):
            upsampled_data[start_idx_fine_0:end_idx_fine_0, (start_idx_fine_1 + s):end_idx_fine_1:r[1]] = \
                c_1[s, 0]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 0:(-stencil_size + 1)] \
                + c_1[s, 1]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 1:(-stencil_size + 2)] \
                + c_1[s, 2]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 2:(-stencil_size + 3)] \
                + c_1[s, 3]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 3:(-stencil_size + 4)] \
                + c_1[s, 4]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 4:(-stencil_size + 5)] \
                + c_1[s, 5]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 5:]
    
    elif dim == 3:
        start_idx_fine_0 = (half_stencil_size - 1)*r[0] + (r[0] + 1)/2
        end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
        
        start_idx_fine_1 = (half_stencil_size - 1)*r[1] + (r[1] + 1)/2
        end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
        
        start_idx_fine_2 = (half_stencil_size - 1)*r[2] + (r[2] + 1)/2
        end_idx_fine_2 = upsampled_data_shape[2] - (half_stencil_size - 1)*r[2] - r[2]/2
        
        for s in range(r[0]):
            upsampled_data_0[(start_idx_fine_0 + s):end_idx_fine_0:r[0], : :] = \
                c_0[s, 0]*data_component[0:(-stencil_size + 1), :, :] \
                + c_0[s, 1]*data_component[1:(-stencil_size + 2), :, :] \
                + c_0[s, 2]*data_component[2:(-stencil_size + 3), :, :] \
                + c_0[s, 3]*data_component[3:(-stencil_size + 4), :, :] \
                + c_0[s, 4]*data_component[4:(-stencil_size + 5), :, :] \
                + c_0[s, 5]*data_component[5:, :, :]
        
        for s in range(r[1]):
            upsampled_data_1[start_idx_fine_0:end_idx_fine_0, (start_idx_fine_1 + s):end_idx_fine_1:r[1], :] = \
                c_1[s, 0]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 0:(-stencil_size + 1), :] \
                + c_1[s, 1]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 1:(-stencil_size + 2), :] \
                + c_1[s, 2]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 2:(-stencil_size + 3), :] \
                + c_1[s, 3]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 3:(-stencil_size + 4), :] \
                + c_1[s, 4]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 4:(-stencil_size + 5), :] \
                + c_1[s, 5]*upsampled_data_0[start_idx_fine_0:end_idx_fine_0, 5:, :]
        
        for s in range(r[2]):
            upsampled_data[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, (start_idx_fine_2 + s):end_idx_fine_2:r[2]] = \
                c_2[s, 0]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 0:(-stencil_size + 1)] \
                + c_2[s, 1]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 1:(-stencil_size + 2)] \
                + c_2[s, 2]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 2:(-stencil_size + 3)] \
                + c_2[s, 3]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 3:(-stencil_size + 4)] \
                + c_2[s, 4]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 4:(-stencil_size + 5)] \
                + c_2[s, 5]*upsampled_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, 5:]
    
    # If some of the upsampled nodes overlap with the original nodes, copy the data from the original nodes to the first
    # sets of overlapped nodes since they are not upsampled yet.
    
    if r[0] % 2 == 1: # If some of the upsampled nodes overlap with the original nodes.
        overlap_idx_coarse = half_stencil_size - 1
        overlap_idx_fine = (half_stencil_size - 1)*r[0] + r[0]/2
        
        if dim == 1:
            upsampled_data[overlap_idx_fine] = data_component[overlap_idx_coarse]
        
        elif dim == 2:
            start_idx_fine_0 = (half_stencil_size - 1)*r[0] + r[0]/2
            end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
            
            start_idx_fine_1 = (half_stencil_size - 1)*r[1] + r[1]/2
            end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
            
            upsampled_line_data_0 = numpy.repeat(data_component[overlap_idx_coarse, :], \
                r[1])
            
            upsampled_line_data_1 = numpy.repeat(data_component[:, overlap_idx_coarse], \
                r[0])
            
            upsampled_data[overlap_idx_fine, start_idx_fine_1:end_idx_fine_1] = \
                upsampled_line_data_0[start_idx_fine_1:end_idx_fine_1]
            
            upsampled_data[start_idx_fine_0:end_idx_fine_0, overlap_idx_fine] = \
                upsampled_line_data_1[start_idx_fine_0:end_idx_fine_0]
        
        elif dim == 3:
            start_idx_fine_0 = (half_stencil_size - 1)*r[0] + r[0]/2
            end_idx_fine_0 = upsampled_data_shape[0] - (half_stencil_size - 1)*r[0] - r[0]/2
            
            start_idx_fine_1 = (half_stencil_size - 1)*r[1] + r[1]/2
            end_idx_fine_1 = upsampled_data_shape[1] - (half_stencil_size - 1)*r[1] - r[1]/2
            
            start_idx_fine_2 = (half_stencil_size - 1)*r[2] + r[2]/2
            end_idx_fine_2 = upsampled_data_shape[2] - (half_stencil_size - 1)*r[2] - r[2]/2
            
            upsampled_face_data_0 = numpy.repeat(data_component[overlap_idx_coarse, :, :], \
                r[1], axis = 0)
            upsampled_face_data_0 = numpy.repeat(upsampled_face_data_0, r[2], axis = 1)
            
            upsampled_face_data_1 = numpy.repeat(data_component[:, overlap_idx_coarse, :], \
                r[0], axis = 0)
            upsampled_face_data_1 = numpy.repeat(upsampled_face_data_1, r[2], axis = 1)
            
            upsampled_face_data_2 = numpy.repeat(data_component[:, :, overlap_idx_coarse], \
                r[0], axis = 0)
            upsampled_face_data_2 = numpy.repeat(upsampled_face_data_2, r[1], axis = 1)
            
            upsampled_data[overlap_idx_fine, start_idx_fine_1:end_idx_fine_1, start_idx_fine_2:end_idx_fine_2] = \
                upsampled_face_data_0[start_idx_fine_1:end_idx_fine_1, start_idx_fine_2:end_idx_fine_2]
            
            upsampled_data[start_idx_fine_0:end_idx_fine_0, overlap_idx_fine, start_idx_fine_2:end_idx_fine_2] = \
                upsampled_face_data_1[start_idx_fine_0:end_idx_fine_0, start_idx_fine_2:end_idx_fine_2]
            
            upsampled_data[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1, overlap_idx_fine] = \
                upsampled_face_data_2[start_idx_fine_0:end_idx_fine_0, start_idx_fine_1:end_idx_fine_1]
    
    return upsampled_data

