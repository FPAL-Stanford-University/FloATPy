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
        serial_reader : a concrete object that extends BaseReader
        sub_domain : Tuple of size 2 with the first entry being lo and second entry being hi
        num_ghosts : numpy integer array of size 3 with the no. of ghost values in the x, y and z directions respectively
        """
        
        if not isinstance(serial_reader, BaseReader):
            raise RuntimeError("The given serial data reader is not instance of the base reader!")
       
        # Set the communicator and it's Fortran value
        self._comm  = comm
        self._fcomm = comm.py2f()

        # Set the serial data reader to use
        self._serial_reader = serial_reader
       
        # Dimensionality of the data set (1D, 2D or 3D)
        self._dim = serial_reader.dimension
        
        self._periodic_dimensions = None
        self._domain_size = None
        
        if num_ghosts is None:
            self._num_ghosts = numpy.array([0, 0, 0], dtype=numpy.int32)
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
        
        elif self._dim == 3:
            self._periodic_dimensions = numpy.asarray(serial_reader.periodic_dimensions)
            self._domain_size = numpy.asarray(serial_reader.domain_size)
       
        if sub_domain == None:
            self._subdomain_lo = numpy.array([0, 0, 0], dtype=self._domain_size.dtype)
            self._subdomain_hi = self._domain_size
            self._subdomain_size = self._domain_size
        else:
            # Need to change periodic_dimensions if only reading in a sub domain
            raise NotImplementedError('Reading in only a partial sub domain is not yet implemented. Sorry!')
            self._subdomain_lo = numpy.asarray( sub_domain[0] )
            self._subdomain_hi = numpy.asarray( sub_domain[1] )
            self._subdomain_size = self._subdomain_hi - self._subdomain_lo
      
        # Create the parallel grid partition object that handles all the communication stuff.
        self.grid_partition = pyt3d.t3dmod.t3d(self._fcomm, self._subdomain_size[0], self._subdomain_size[1], self._subdomain_size[2],
                                               self._periodic_dimensions, nghosts=self._num_ghosts )

        self._interior_chunk_sz = numpy.zeros(3, dtype=numpy.int32, order='F') # Size of the interior chunk of this process
        self._interior_chunk_lo = numpy.zeros(3, dtype=numpy.int32, order='F') # Index of the start of the interior chunk of this process
        self._interior_chunk_hi = numpy.zeros(3, dtype=numpy.int32, order='F') # Index of the end of the interior chunk of this process
        self.grid_partition.get_sz3d(self._interior_chunk_sz)
        self.grid_partition.get_st3d(self._interior_chunk_lo)
        self.grid_partition.get_en3d(self._interior_chunk_hi)
        self._interior_chunk_lo = self._interior_chunk_lo - 1 # Convert to 0 based indexing

        self._full_chunk_sz = numpy.zeros(3, dtype=numpy.int32, order='F') # Size of the full chunk of this process
        self._full_chunk_lo = numpy.zeros(3, dtype=numpy.int32, order='F') # Index of the start of the full chunk of this process
        self._full_chunk_hi = numpy.zeros(3, dtype=numpy.int32, order='F') # Index of the end of the full chunk of this process
        self.grid_partition.get_sz3d(self._full_chunk_sz)
        self.grid_partition.get_st3d(self._full_chunk_lo)
        self.grid_partition.get_en3d(self._full_chunk_hi)
        self._full_chunk_lo = self._full_chunk_lo - 1 # Convert to 0 based indexing

        # Set the sub domain to read in using the serial data reader.
        self.serial_reader.sub_domain = ( tuple(self._interior_chunk_lo), tuple(self._interior_chunk_hi) )

    
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
    
    
    def getSubDomain(self):
        """
        Return two tuples containing the sub-domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """
        
        return self._subdomain_lo, self._subdomain_hi
    
    
    sub_domain = property(getSubDomain, setSubDomain)
    
    
    @property
    def interior_chunk(self):
        """
        Return two tuples containing the interior chunk of sub-domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """
        
        return tuple(self._interior_chunk_lo), tuple(self._interior_chunk_hi)
    
    
    @property
    def interior_chunk_size(self):
        """
        Return a tuple containing the size of the interior chunk of sub-domain used in this reader.
        """
        
        return tuple(self._interior_chunk_sz)
    
    
    @property
    def full_chunk(self):
        """
        Return two tuples containing the full chunk of sub-domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """
        
        return tuple(self._full_chunk_lo), tuple(self._full_chunk_hi)
    
    
    @property
    def full_chunk_size(self):
        """
        Return a tuple containing the size of the full chunk of sub-domain used in this reader.
        """
        
        return tuple(self._full_chunk_sz)
    
    
    def readCoordinates(self):
        """
        Get the coordinates of the stored sub-domain.
        Default to the full domain when the sub-domain is not set.
        (Not yet implemented!)
        """
        
        return
    
    
    def readData(self, var_names):
        """
        Read the data of several variables in the assigned chunk of the stored sub-domain.
        Default to the full domain when the sub-domain is not set.
        (Not yet implemented!)
        """
        
        return
