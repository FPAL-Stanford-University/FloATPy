from base_reader import BaseReader
from mpi4py import MPI
import numpy
import os
import sys

cwd = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(cwd, '../parallel/pyt3d')))
import pyt3d

class ParallelDataReader():
    """
    Class to read data and exchange data across nodes with MPI.
    """
    
    def __init__(self, comm, serial_reader, sub_domain=None, num_ghosts=None):
        """
        Constructor of the class.
        
        comm : mpi4py communicator object
        serial_reader : a concrete object that extends BaseReader (Do not use this outside of this class)
        sub_domain : Iterable of size 2 with the first entry being lo and second entry being hi
        num_ghosts : numpy integer array of size 3 with the no. of ghost values in the x, y and z directions respectively
        """
        
        if not isinstance(serial_reader, BaseReader):
            raise RuntimeError("The given serial data reader is not instance of the base reader!")
       
        if serial_reader.data_order != 'F':
            raise RuntimeError("The data order should be 'F'!")
        
        # Set the communicator and it's Fortran value.
        self._comm  = comm
        self._fcomm = comm.py2f()
        
        # Set the serial data reader to use.
        self._serial_reader = serial_reader
        
        # Dimensionality of the data set (1D, 2D or 3D)
        self._dim = serial_reader.dimension
        
        if self._dim < 0 or self._dim > 3:
            raise RuntimeError('Dimension of data should be between 1 and 3!')
        
        self._periodic_dimensions = None
        self._domain_size = None
        
        if num_ghosts is None:
            self._num_ghosts = numpy.array([0, 0, 0], dtype=numpy.int32)
        else:
            if self._dim == 1:
                if len(num_ghosts) < 1 or len(num_ghosts) > 3:
                    raise RuntimeError('Dimension of num_ghosts should be between 1 and 3!')
                self._num_ghosts = numpy.array([num_ghosts[0], 0, 0], dtype=numpy.int32)
                
            elif self._dim == 2:
                if len(num_ghosts) < 2 or len(num_ghosts) > 3:
                    raise RuntimeError('Dimension of num_ghosts should be between 2 and 3!')
                self._num_ghosts = numpy.array([num_ghosts[0], num_ghosts[1], 0], dtype=numpy.int32)
                
            else:
                if len(num_ghosts) != 3:
                    raise RuntimeError('Dimension of num_ghosts should be 3!')
                self._num_ghosts = numpy.asarray(num_ghosts, dtype=numpy.int32)
        
        if self._dim == 1:
            self._periodic_dimensions = numpy.array([serial_reader.dimension[0], True, True])
            self._domain_size = numpy.array([serial_reader.domain_size[0], 1, 1])
            self._num_ghosts[1:] = 0
        
        elif self._dim == 2:
            self._periodic_dimensions = numpy.array([serial_reader.dimension[0], serial_reader.dimension[1] , True])
            self._domain_size = numpy.array([serial_reader.domain_size[0], serial_reader.domain_size[1], 1])
            self._num_ghosts[2] = 0
        
        else:
            self._periodic_dimensions = numpy.asarray(serial_reader.periodic_dimensions)
            self._domain_size = numpy.asarray(serial_reader.domain_size)
       
        if sub_domain == None:
            self._subdomain_lo = numpy.array([0, 0, 0], dtype=self._domain_size.dtype)
            self._subdomain_hi = self._domain_size - 1
            self._subdomain_size = self._domain_size
        else:
            # Need to change periodic_dimensions if only reading in a sub-domain
            raise NotImplementedError('Reading in only a partial sub domain is not yet implemented. Sorry!')
            
            try:
                lo, hi = sub_domain
            except ValueError:
                raise ValueError("Pass an iterable of sub_domain with two items!")
            
            for i in range(self._dim):
                if lo[i] < 0 or lo[i] > self._domain_size[i]:
                    raise ValueError('Invalid indices in sub-domain. Cannot be < 0 or > domain size!')
                if hi[i] < 0 or hi[i] > self._domain_size[i]:
                    raise ValueError('Invalid indices in sub-domain. Cannot be < 0 or > domain size!')
                if hi[i] < lo[i]:
                    raise ValueError('Invalid indices in sub-domain. Upper bound cannot be smaller than lower bound!')
            
            self._subdomain_lo = numpy.asarray(lo)
            self._subdomain_hi = numpy.asarray(hi)
            self._subdomain_size = self._subdomain_hi - self._subdomain_lo + 1
        
        # Create the parallel grid partition object that handles all the communication stuff.
        self._grid_partition = pyt3d.t3dmod.t3d(self._fcomm, \
                                               self._subdomain_size[0], self._subdomain_size[1], self._subdomain_size[2], \
                                               self._periodic_dimensions, nghosts=self._num_ghosts )
        
        # Size of the interior chunk of this process.
        self._interior_chunk_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        # Indices of the start and end of the interior chunk of this process.
        self._interior_chunk_lo = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._interior_chunk_hi = numpy.zeros(3, dtype=numpy.int32, order='F')
        
        self._grid_partition.get_sz3d(self._interior_chunk_size)
        self._grid_partition.get_st3d(self._interior_chunk_lo)
        self._grid_partition.get_en3d(self._interior_chunk_hi)
        self._interior_chunk_lo = self._interior_chunk_lo - 1 # Convert to 0 based indexing
        self._interior_chunk_hi = self._interior_chunk_hi - 1 # Convert to 0 based indexing
        
        # Size of the full chunk of this process.
        self._full_chunk_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        # Indices of the start and end of the full chunk of this process.
        self._full_chunk_lo = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._full_chunk_hi = numpy.zeros(3, dtype=numpy.int32, order='F')
        
        self._grid_partition.get_sz3dg(self._full_chunk_size)
        self._grid_partition.get_st3dg(self._full_chunk_lo)
        self._grid_partition.get_en3dg(self._full_chunk_hi)
        self._full_chunk_lo = self._full_chunk_lo - 1 # Convert to 0 based indexing
        self._full_chunk_hi = self._full_chunk_hi - 1 # Convert to 0 based indexing
        
        # Set the sub domain to read in using the serial data reader.
        self._serial_reader.sub_domain = ( tuple(self._interior_chunk_lo), tuple(self._interior_chunk_hi) )
        
        self._interior = ( slice(self._num_ghosts[0],self._full_chunk_size[0] - self._num_ghosts[0]),
                           slice(self._num_ghosts[1],self._full_chunk_size[1] - self._num_ghosts[1]),
                           slice(self._num_ghosts[2],self._full_chunk_size[2] - self._num_ghosts[2]) )
    
    
    @property
    def serial_reader(self):
        """
        Return the serial data reader.
        """
        
        return self._serial_reader
    
    
    @property
    def domain_size(self):
        """
        Return a tuple containing the full domain size of this dataset.
        """
        
        return tuple(self._domain_size[0:self._dim])
    
    
    @property
    def periodic_dimensions(self):
        """
        Return a tuple indicating if data is periodic in each dimension.
        """
        
        return tuple(self._periodic_dimensions[0:self._dim])
    
    
    @property
    def time(self):
        """
        Return the simulation time at current time step.
        """
        
        return self._serial_reader.time
    
    
    def setStep(self, step):
        """
        Update the metadata from the summary file in the data directory at a new time step.
        """
        
        self._serial_reader.step = step
    
    
    def getStep(self, step):
        """
        Return the time step that is currently set.
        """
        
        return self._serial_reader.step
    
    
    step = property(getStep, setStep)
    
    
    @property
    def sub_domain(self):
        """
        Return two tuples containing the sub-domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """
        
        return tuple(self._subdomain_lo[0:self._dim]), tuple(self._subdomain_hi[0:self._dim])
    
    
    @property
    def interior_chunk(self):
        """
        Return two tuples containing the interior chunk of sub-domain used in this reader as a lower bound (lo)
        and upper bound (hi).
        """
        
        return tuple(self._interior_chunk_lo[0:self._dim]), tuple(self._interior_chunk_hi[0:self._dim])
    
    
    @property
    def interior_chunk_size(self):
        """
        Return a tuple containing the size of the interior chunk of sub-domain used in this reader.
        """
        
        return tuple(self._interior_chunk_size[0:self._dim])
    
    
    @property
    def full_chunk(self):
        """
        Return two tuples containing the full chunk of sub-domain used in this reader as a lower bound (lo) and
        upper bound (hi).
        """
        
        return tuple(self._full_chunk_lo[0:self._dim]), tuple(self._full_chunk_hi[0:self._dim])
    
    
    @property
    def full_chunk_size(self):
        """
        Return a tuple containing the size of the full chunk of sub-domain used in this reader.
        """
        
        return tuple(self._full_chunk_size[0:self._dim])
   
    
    @property
    def interior(self):
        """
        Return a boolean numpy array which is True only in the interior of the domain and False for ghost cells.
        """
        
        return self._interior[0:self._dim]
    
    
    def readCoordinates(self, communicate=False):
        """
        Get the coordinates of the stored sub-domain.
        Default to the full domain when the sub-domain is not set.
        (Not yet well implemented with communication!)
        """
        
        x_c = numpy.zeros( tuple(self._full_chunk_size), dtype=numpy.float64, order='F' )
        y_c = numpy.zeros( tuple(self._full_chunk_size), dtype=numpy.float64, order='F' )
        z_c = numpy.zeros( tuple(self._full_chunk_size), dtype=numpy.float64, order='F' )
        
        if self._dim == 1:
            x_coords, y_coords, z_coords = self._serial_reader.readCoordinates()
            nx = x_coords.shape[0]
            x_coords = x_coords.reshape((nx, 1, 1), order='F')
            x_c[self._interior] = x_coords
        
        elif self._dim == 2:
            x_coords, y_coords, z_coords = self._serial_reader.readCoordinates()
            
            nx = x_coords.shape[0]
            ny = x_coords.shape[1]
            
            x_coords = x_coords.reshape((nx, ny, 1), order='F')
            y_coords = y_coords.reshape((nx, ny, 1), order='F')
            
            x_c[self._interior] = x_coords
            y_c[self._interior] = y_coords
        
        else:
            x_c[self._interior], y_c[self._interior], z_c[self._interior] = self._serial_reader.readCoordinates()
        
        # Communicate to get the coordinates in the ghost cell regions.
        if communicate:
            self._grid_partition.fill_halo_x(x_c)
            
            if self._dim > 1:
                self._grid_partition.fill_halo_y(x_c)
                
                self._grid_partition.fill_halo_x(y_c)
                self._grid_partition.fill_halo_y(y_c)
            
            if self._dim > 2:
                self._grid_partition.fill_halo_z(x_c)
                self._grid_partition.fill_halo_z(y_c)
                
                self._grid_partition.fill_halo_x(z_c)
                self._grid_partition.fill_halo_y(z_c)
                self._grid_partition.fill_halo_z(z_c)
        
        if self._dim == 1:
            return numpy.squeeze(x_c, (1, 2)), [], []
        elif self._dim == 2:
            return numpy.squeeze(x_c, 2), numpy.squeeze(y_c, 2), []
        else:
            return x_c, y_c, z_c
    
    
    def readData(self, var_names, communicate=False):
        """
        Read the data of several variables in the assigned chunk of the stored sub-domain.
        Default to the full domain when the sub-domain is not set.
        (Not yet well implemented with communication and vector!)
        """
        
        if isinstance(var_names, basestring):
            var_names = (var_names,)
        
        data_vars = []
        
        for i in range(len(var_names)):
            data_var, = self._serial_reader.readData(var_names[i])
            
            num_components = 1
            if data_var.ndim == self._dim + 1:
                num_components = data_var.shape[self._dim]
            
            if num_components == 1:
                data_vars.append( numpy.zeros( tuple(self._full_chunk_size[0:self._dim]), dtype=numpy.float64, order='F' ) )
                data_vars[i][self._interior[0:self._dim]] = data_var
            
            else:
                data_vars.append( numpy.zeros( tuple(self._full_chunk_size[0:self._dim]) + (num_components, ), \
                                               dtype=numpy.float64, order='F' ) )
                data_vars[i][ self._interior[0:self._dim] + (slice(0, num_components), ) ] = data_var
            
            # Communicate to get the data in the ghost cell regions.
            if communicate:
                shape_modified = data_vars[i].shape[0:self._dim]
                
                if self._dim == 1:
                    shape_modified = numpy.append(shape_modified, [1, 1])
                elif self._dim == 2:
                    shape_modified = numpy.append(shape_modified, 1)
                
                data_to_communicate = []
                for ic in range(0, num_components):
                    if self._dim == 1:
                        if num_components == 1:
                            data_to_communicate = numpy.reshape(data_vars[i][:], shape_modified)
                        else:
                            data_to_communicate = numpy.reshape(data_vars[i][:, ic], shape_modified)
                        
                        self._grid_partition.fill_halo_x(data_to_communicate)
                    
                    elif self._dim == 2:
                        if num_components == 1:
                            data_to_communicate = numpy.reshape(data_vars[i][:, :], shape_modified)
                        else:
                            data_to_communicate = numpy.reshape(data_vars[i][:, :, ic], shape_modified)
                        
                        self._grid_partition.fill_halo_x(data_to_communicate)
                        self._grid_partition.fill_halo_y(data_to_communicate)
                    
                    else:
                        if num_components == 1:
                            data_to_communicate = data_vars[i]
                        else:
                            data_to_communicate = data_vars[i][:, :, :, ic]
                        
                        self._grid_partition.fill_halo_x(data_to_communicate)
                        self._grid_partition.fill_halo_y(data_to_communicate)
                        self._grid_partition.fill_halo_z(data_to_communicate)
        
        return tuple(data_vars)
