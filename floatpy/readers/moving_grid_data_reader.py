import h5py
import math
import numpy

from base_reader import BaseReader

class MovingGridDataReader(BaseReader):
    """
    Class to read data.
    """
    
    def __init__(self, data_directory_path, data_order = 'C'):
        """
        Constructor of the class.
        """
        
        self._data_directory_path = data_directory_path
        
        # Get the summary file.
        
        full_summary_path = self._data_directory_path + '/summary.h5'
        summary_f = h5py.File(full_summary_path, 'r')
        
        # Get the periodicity.
        
        self._periodic_dimensions = summary_f['periodic_dimensions'].value
        
        # Get the dimension.
        
        self._dimension = self._periodic_dimensions.shape[0]
        
        # Get the domain size.
        
        self._domain_size = summary_f['domain_size'].value
        
        # Get the variable names.
        
        self._var_names = summary_f['var_names'].value
        
        # Get the number of components in each variable.
        
        self._var_number_components = summary_f['var_number_components'].value
        
        # Get all the steps.
        
        self._steps = summary_f['steps'].value
        
        # Set the current time step to be first available step.
        
        self._step = self._steps[0]
        
        # Get the grid partition.
        
        self._grid_partition = summary_f['grid_partition'].value
        
        # Set the data order.
        
        if data_order == 'C':
            self._data_order = 'C'
        elif data_order == 'F':
            self._data_order = 'F'
        else:
            raise RuntimeError("Unknown data order '" + data_order + "'!")
        
        summary_f.close()
        
        full_step_summary_path = self._data_directory_path + '/n_' + str(self._step).zfill(5) + '/step_summary.h5'
        
        step_summary_f = h5py.File(full_step_summary_path, 'r')
        
        # Get the time.
        
        self._time = step_summary_f['time'].value[0]
        
        # Get the physical domain.
        
        self._x_lo = step_summary_f['x_lo'].value
        self._x_hi = step_summary_f['x_hi'].value
        
        step_summary_f.close()
        
        # Get the sub-domain.
        
        self._subdomain_lo = numpy.zeros(self._dimension, dtype=int)
        self._subdomain_hi = self._domain_size - 1
    
    
    def setStep(self, step):
        """
        Set the step.
        """
        
        assert (step in self._steps), "Step to read in is not available in the dataset."
        self._step = step
        
        full_step_summary_path = self._data_directory_path + '/n_' + str(self._step).zfill(5) + '/step_summary.h5'
        
        step_summary_f = h5py.File(full_step_summary_path, 'r')
        
        # Get the time.
        
        self._time = step_summary_f['time'].value[0]
        
        # Get the physical domain.
        
        self._x_lo = step_summary_f['x_lo'].value
        self._x_hi = step_summary_f['x_hi'].value
        
        step_summary_f.close()
    
    
    def getStep(self):
        """
        Get the step.
        """
        
        return self._step
    
    
    step = property(getStep, setStep)
    
    
    def setSubDomain(self, lo_and_hi):
        """
        Set the sub-domain for reading coordinates and data.
        """
        
        try:
            lo, hi = lo_and_hi
        except ValueError:
            raise ValueError("Pass an iterable with two items!")
        
        assert (len(lo) == self._dimension), "Length of lo is not equal to the dimension."
        assert (len(hi) == self._dimension), "Length of hi is not equal to the dimension."
        
        for i in range(self._dimension):
            if lo[i] < 0 or lo[i] >= self._domain_size[i]:
                raise ValueError('Invalid indices in chunk. Cannot be < 0 or >= domain size!')
            if hi[i] < 0 or hi[i] >= self._domain_size[i]:
                raise ValueError('Invalid indices in chunk. Cannot be < 0 or >= domain size!')
            if hi[i] < lo[i]:
                raise ValueError('Invalid indices in chunk. Upper bound cannot be smaller than lower bound!')
        
        for i in range(self._dimension):
            self._subdomain_lo[i] = lo[i]
            self._subdomain_hi[i] = hi[i]
    
    
    def getSubDomain(self):
        """
        Return two tuples containing the sub-domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """
        
        lo = tuple(self._subdomain_lo)
        hi = tuple(self._subdomain_hi)
        
        return lo, hi
    
    
    sub_domain = property(getSubDomain, setSubDomain)
    
    
    @property
    def physical_domain(self):
        """
        Return two tuples containing the physical domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """
        
        lo = tuple(self._x_lo)
        hi = tuple(self._x_hi)
        
        return lo, hi
    
    
    @property
    def domain_size(self):
        """
        Return a tuple containing the full domain size of this dataset.
        """
        
        return tuple(self._domain_size)
    
    
    @property
    def dimension(self):
        """
        Return the dimension of the domain.
        """
        
        return self._dimension
    
    
    @property
    def periodic_dimensions(self):
        """
        Return a tuple indicating if data is periodic in each dimension.
        """
        
        return tuple(self._periodic_dimensions)
    
    
    @property
    def time(self):
        """
        Return the simulation time at current time step.
        """
        
        return self._time
    
    
    @property
    def data_order(self):
        """
        Return the data order.
        """
        
        return self._data_order
    
    
    @property
    def steps(self):
        """
        Return all of the steps.
        """
        
        return list(self._steps)
    
    
    def readCoordinates(self):
        """
        Method to return the X, Y and Z coordinates.
        """
        
        subdomain_size = self._subdomain_hi - self._subdomain_lo + 1
        
        if self._dimension == 1:
            dx = float(self._x_hi[0] - self._x_lo[0])/float(self._domain_size[0])
            
            x_c = numpy.linspace(self._subdomain_lo[0]*dx + self._x_lo[0] + 0.5*dx, \
                                 self._subdomain_hi[0]*dx + self._x_lo[0] + 0.5*dx, \
                                 num=subdomain_size[0])
            
            return x_c
        
        elif self._dimension == 2:
            dx = float(self._x_hi[0] - self._x_lo[0])/float(self._domain_size[0])
            dy = float(self._x_hi[1] - self._x_lo[1])/float(self._domain_size[1])
            
            x_c = numpy.linspace(self._subdomain_lo[0]*dx + self._x_lo[0] + 0.5*dx, \
                                 self._subdomain_hi[0]*dx + self._x_lo[0] + 0.5*dx, \
                                 num=subdomain_size[0])
            
            y_c = numpy.linspace(self._subdomain_lo[1]*dy + self._x_lo[1] + 0.5*dy, \
                                 self._subdomain_hi[1]*dy + self._x_lo[1] + 0.5*dy, \
                                 num=subdomain_size[1])
            
            X_c, Y_c = numpy.meshgrid(x_c, y_c, sparse=False, indexing='ij')
            
            return X_c. Y_c
        
        elif self._dimension == 3:
            dx = float(self._x_hi[0] - self._x_lo[0])/float(self._domain_size[0])
            dy = float(self._x_hi[1] - self._x_lo[1])/float(self._domain_size[1])
            dz = float(self._x_hi[2] - self._x_lo[2])/float(self._domain_size[2])
            
            x_c = numpy.linspace(self._subdomain_lo[0]*dx + self._x_lo[0] + 0.5*dx, \
                                 self._subdomain_hi[0]*dx + self._x_lo[0] + 0.5*dx, \
                                 num=subdomain_size[0])
            
            y_c = numpy.linspace(self._subdomain_lo[1]*dy + self._x_lo[1] + 0.5*dy, \
                                 self._subdomain_hi[1]*dy + self._x_lo[1] + 0.5*dy, \
                                 num=subdomain_size[1])
            
            z_c = numpy.linspace(self._subdomain_lo[2]*dz + self._x_lo[2] + 0.5*dz, \
                                 self._subdomain_hi[2]*dz + self._x_lo[2] + 0.5*dz, \
                                 num=subdomain_size[2])
            
            X_c, Y_c, Z_c = numpy.meshgrid(x_c, y_c, z_c, sparse=False, indexing='ij')
            
            return X_c, Y_c, Z_c
        ppend(subdomain_size, var_num_components)
    
    def readData(self, var_names, data=None):
        """
        Method to read a subdomain of data for variables at current step.
        """
        
        # If a simple string is passed in, convert to a tuple.
        if isinstance(var_names, basestring):
            var_names = (var_names,)
        
        subdomain_size = self._subdomain_hi - self._subdomain_lo + 1
        
        _data = []
        
        for var_i in range(len(var_names)):
            var_name = var_names[var_i]
            assert (var_name in self._var_names), "Variable to read is not available in the dataset."
            
            var_idx = numpy.where(self._var_names == var_name)[0][0]
            var_num_components = self._var_number_components[var_idx]
            
            if var_num_components > 1:
                var_size = numpy.insert(subdomain_size, 0, var_num_components)
            else:
                var_size = subdomain_size
            
            data_var = numpy.empty(var_size, order = 'C')
            
            if self._dimension == 1:
                chunk_size_x = int(math.ceil(float(self._domain_size[0])/float(self._grid_partition[0])))
                
                for gp_j in range(self._grid_partition[1]):
                    for gp_i in range(self._grid_partition[0]):
                        gp_linear_idx = gp_i
                        chunk_lo = numpy.array([chunk_size_x*gp_i])
                        chunk_hi = numpy.array([chunk_size_x*(gp_i+1)-1])
                        
                        chunk_hi = numpy.minimum(chunk_hi, self._domain_size-1)
                        
                        if self._subdomain_lo[0] > chunk_hi[0] or self._subdomain_hi[0] < chunk_lo[0]:
                            continue
                        
                        full_proc_path = self._data_directory_path + '/n_' + str(self._step).zfill(5) \
                                         + '/proc_' + str(gp_linear_idx).zfill(5) + '/data.h5'
                        
                        proc_f = h5py.File(full_proc_path, 'r')
                        
                        data_f = proc_f['data']
                        
                        if var_num_components > 1:
                            data_var[ :,
                                      max(0, chunk_lo[0]-self._subdomain_lo[0]):
                                      min(self._subdomain_hi[0]-self._subdomain_lo[0], chunk_hi[0]-self._subdomain_lo[0])+1 ] = \
                            data_f[var_name][ :,
                                              max(0, self._subdomain_lo[0]-chunk_lo[0]):
                                              min(chunk_hi[0]-chunk_lo[0], self._subdomain_hi[0]-chunk_lo[0])+1 ]
                        
                        else:
                            data_var[ max(0, chunk_lo[0]-self._subdomain_lo[0]):
                                      min(self._subdomain_hi[0]-self._subdomain_lo[0], chunk_hi[0]-self._subdomain_lo[0])+1 ] = \
                            data_f[var_name][ max(0, self._subdomain_lo[0]-chunk_lo[0]):
                                              min(chunk_hi[0]-chunk_lo[0], self._subdomain_hi[0]-chunk_lo[0])+1 ]
                        
                        proc_f.close()
                
                if self._data_order == 'F':
                    if data == None:
                        if var_num_components > 1:
                            data_var_F = numpy.empty(numpy.append(subdomain_size, var_num_components), order='F')
                            for component_i in range(var_num_components):
                                data_var_F[:, component_i] = numpy.asfortranarray(data_var[component_i, :])
                            _data.append(data_var_F)
                        else:
                            _data.append(numpy.asfortranarray(data_var))
                    else:
                        if var_num_components > 1:
                            for component_i in range(var_num_components):
                                data[var_i][:, component_i] = numpy.asfortranarray(data_var[component_i, :])
                        else:
                            data[var_i] = numpy.asfortranarray(data_var)
                else:
                    if data == None:
                        _data.append(data_var)
                    else:
                        data[var_i] = data_var
            
            elif self._dimension == 2:
                chunk_size_x = int(math.ceil(float(self._domain_size[0])/float(self._grid_partition[0])))
                chunk_size_y = int(math.ceil(float(self._domain_size[1])/float(self._grid_partition[1])))
                
                for gp_j in range(self._grid_partition[1]):
                    for gp_i in range(self._grid_partition[0]):
                        gp_linear_idx = gp_i \
                                        + gp_j*self._grid_partition[0]
                        chunk_lo = numpy.array([chunk_size_x*gp_i, chunk_size_y*gp_j])
                        chunk_hi = numpy.array([chunk_size_x*(gp_i+1)-1, chunk_size_y*(gp_j+1)-1])
                        
                        chunk_hi = numpy.minimum(chunk_hi, self._domain_size-1)
                        
                        if self._subdomain_lo[0] > chunk_hi[0] or self._subdomain_hi[0] < chunk_lo[0]:
                            continue
                        if self._subdomain_lo[1] > chunk_hi[1] or self._subdomain_hi[1] < chunk_lo[1]:
                            continue
                        
                        full_proc_path = self._data_directory_path + '/n_' + str(self._step).zfill(5) \
                                         + '/proc_' + str(gp_linear_idx).zfill(5) + '/data.h5'
                        
                        proc_f = h5py.File(full_proc_path, 'r')
                        
                        data_f = proc_f['data']
                        
                        if var_num_components > 1:
                            data_var[ :,
                                      max(0, chunk_lo[0]-self._subdomain_lo[0]):
                                      min(self._subdomain_hi[0]-self._subdomain_lo[0], chunk_hi[0]-self._subdomain_lo[0])+1,
                                      max(0, chunk_lo[1]-self._subdomain_lo[1]):
                                      min(self._subdomain_hi[1]-self._subdomain_lo[1], chunk_hi[1]-self._subdomain_lo[1])+1 ] = \
                            data_f[var_name][ :,
                                              max(0, self._subdomain_lo[0]-chunk_lo[0]):
                                              min(chunk_hi[0]-chunk_lo[0], self._subdomain_hi[0]-chunk_lo[0])+1,
                                              max(0, self._subdomain_lo[1]-chunk_lo[1]):
                                              min(chunk_hi[1]-chunk_lo[1], self._subdomain_hi[1]-chunk_lo[1])+1 ]
                        
                        else:
                            data_var[ max(0, chunk_lo[0]-self._subdomain_lo[0]):
                                      min(self._subdomain_hi[0]-self._subdomain_lo[0], chunk_hi[0]-self._subdomain_lo[0])+1,
                                      max(0, chunk_lo[1]-self._subdomain_lo[1]):
                                      min(self._subdomain_hi[1]-self._subdomain_lo[1], chunk_hi[1]-self._subdomain_lo[1])+1 ] = \
                            data_f[var_name][ max(0, self._subdomain_lo[0]-chunk_lo[0]):
                                              min(chunk_hi[0]-chunk_lo[0], self._subdomain_hi[0]-chunk_lo[0])+1,
                                              max(0, self._subdomain_lo[1]-chunk_lo[1]):
                                              min(chunk_hi[1]-chunk_lo[1], self._subdomain_hi[1]-chunk_lo[1])+1 ]
                        
                        proc_f.close()
                
                if self._data_order == 'F':
                    if data == None:
                        if var_num_components > 1:
                            data_var_F = numpy.empty(numpy.append(subdomain_size, var_num_components), order='F')
                            for component_i in range(var_num_components):
                                data_var_F[:, :, component_i] = numpy.asfortranarray(data_var[component_i, :, :])
                            _data.append(data_var_F)
                        else:
                            _data.append(numpy.asfortranarray(data_var))
                    else:
                        if var_num_components > 1:
                            for component_i in range(var_num_components):
                                data[var_i][:, :, component_i] = numpy.asfortranarray(data_var[component_i, :, :])
                        else:
                            data[var_i] = numpy.asfortranarray(data_var)
                else:
                    if data == None:
                        _data.append(data_var)
                    else:
                        data[var_i] = data_var
            
            elif self._dimension == 3:
                chunk_size_x = int(math.ceil(float(self._domain_size[0])/float(self._grid_partition[0])))
                chunk_size_y = int(math.ceil(float(self._domain_size[1])/float(self._grid_partition[1])))
                chunk_size_z = int(math.ceil(float(self._domain_size[2])/float(self._grid_partition[2])))
                
                for gp_k in range(self._grid_partition[2]):
                    for gp_j in range(self._grid_partition[1]):
                        for gp_i in range(self._grid_partition[0]):
                            gp_linear_idx = gp_i \
                                            + gp_j*self._grid_partition[0] \
                                            + gp_k*self._grid_partition[0]*self._grid_partition[1]
                            chunk_lo = numpy.array([chunk_size_x*gp_i, chunk_size_y*gp_j, chunk_size_z*gp_k])
                            chunk_hi = numpy.array([chunk_size_x*(gp_i+1)-1, chunk_size_y*(gp_j+1)-1, chunk_size_z*(gp_k+1)-1])
                            
                            chunk_hi = numpy.minimum(chunk_hi, self._domain_size-1)
                            
                            if self._subdomain_lo[0] > chunk_hi[0] or self._subdomain_hi[0] < chunk_lo[0]:
                                continue
                            if self._subdomain_lo[1] > chunk_hi[1] or self._subdomain_hi[1] < chunk_lo[1]:
                                continue
                            if self._subdomain_lo[2] > chunk_hi[2] or self._subdomain_hi[2] < chunk_lo[2]:
                                continue
                            
                            full_proc_path = self._data_directory_path + '/n_' + str(self._step).zfill(5) \
                                             + '/proc_' + str(gp_linear_idx).zfill(5) + '/data.h5'
                            
                            proc_f = h5py.File(full_proc_path, 'r')
                            
                            data_f = proc_f['data']
                            
                            if var_num_components > 1:
                                data_var[ :,
                                          max(0, chunk_lo[0]-self._subdomain_lo[0]):
                                          min(self._subdomain_hi[0]-self._subdomain_lo[0], chunk_hi[0]-self._subdomain_lo[0])+1,
                                          max(0, chunk_lo[1]-self._subdomain_lo[1]):
                                          min(self._subdomain_hi[1]-self._subdomain_lo[1], chunk_hi[1]-self._subdomain_lo[1])+1,
                                          max(0, chunk_lo[2]-self._subdomain_lo[2]):
                                          min(self._subdomain_hi[2]-self._subdomain_lo[2], chunk_hi[2]-self._subdomain_lo[2])+1 ] = \
                                data_f[var_name][ :,
                                                  max(0, self._subdomain_lo[0]-chunk_lo[0]):
                                                  min(chunk_hi[0]-chunk_lo[0], self._subdomain_hi[0]-chunk_lo[0])+1,
                                                  max(0, self._subdomain_lo[1]-chunk_lo[1]):
                                                  min(chunk_hi[1]-chunk_lo[1], self._subdomain_hi[1]-chunk_lo[1])+1,
                                                  max(0, self._subdomain_lo[2]-chunk_lo[2]):
                                                  min(chunk_hi[2]-chunk_lo[2], self._subdomain_hi[2]-chunk_lo[2])+1 ]
                            
                            else:
                                data_var[ max(0, chunk_lo[0]-self._subdomain_lo[0]):
                                          min(self._subdomain_hi[0]-self._subdomain_lo[0], chunk_hi[0]-self._subdomain_lo[0])+1,
                                          max(0, chunk_lo[1]-self._subdomain_lo[1]):
                                          min(self._subdomain_hi[1]-self._subdomain_lo[1], chunk_hi[1]-self._subdomain_lo[1])+1,
                                          max(0, chunk_lo[2]-self._subdomain_lo[2]):
                                          min(self._subdomain_hi[2]-self._subdomain_lo[2], chunk_hi[2]-self._subdomain_lo[2])+1 ] = \
                                data_f[var_name][ max(0, self._subdomain_lo[0]-chunk_lo[0]):
                                                  min(chunk_hi[0]-chunk_lo[0], self._subdomain_hi[0]-chunk_lo[0])+1,
                                                  max(0, self._subdomain_lo[1]-chunk_lo[1]):
                                                  min(chunk_hi[1]-chunk_lo[1], self._subdomain_hi[1]-chunk_lo[1])+1,
                                                  max(0, self._subdomain_lo[2]-chunk_lo[2]):
                                                  min(chunk_hi[2]-chunk_lo[2], self._subdomain_hi[2]-chunk_lo[2])+1 ]
                            
                            proc_f.close()
                
                if self._data_order == 'F':
                    if data == None:
                        if var_num_components > 1:
                            data_var_F = numpy.empty(numpy.append(subdomain_size, var_num_components), order='F')
                            for component_i in range(var_num_components):
                                data_var_F[:, :, :, component_i] = numpy.asfortranarray(data_var[component_i, :, :, :])
                            _data.append(data_var_F)
                        else:
                            _data.append(numpy.asfortranarray(data_var))
                    else:
                        if var_num_components > 1:
                            for component_i in range(var_num_components):
                                data[var_i][:, :, :, component_i] = numpy.asfortranarray(data_var[component_i, :, :, :])
                        else:
                            data[var_i] = numpy.asfortranarray(data_var)
                else:
                    if data == None:
                        _data.append(data_var)
                    else:
                        data[var_i] = data_var
        
        if data == None:
            return _data


BaseReader.register(MovingGridDataReader)

