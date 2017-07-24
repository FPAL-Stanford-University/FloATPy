from base_reader import BaseReader
from mpi4py import MPI
import numpy
import os
import sys

cwd = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cwd + '/../parallel/pyt3d/')
import pyt3d

comm = MPI.COMM_WORLD
fcomm = MPI.COMM_WORLD.py2f()

class ParallelDataReader():
    """
    Class to read data and exchange data across nodes with MPI.
    """
    
    def __init__(self, serial_reader, num_ghosts = None):
        """
        Constructor of the class.
        """
        
        if not isinstance(serial_reader, BaseReader):
            raise RuntimeError("The given serial data reader is not instance of the base reader!")
        
        self._serial_reader = serial_reader
        
        self._dim = serial_reader.dimension
        
        self._periodic_dimensions = None
        self._domain_size = None
        
        if self._dim == 1:
            self._periodic_dimensions = numpy.array([serial_reader.dimension[0], False, False])
            self._domain_size = numpy.array([serial_reader.domain_size[0], 1, 1])
        
        elif self._dim == 2:
            self._periodic_dimensions = numpy.array([serial_reader.dimension[0], serial_reader.dimension[1] , False])
            self._domain_size = numpy.array([serial_reader.domain_size[0], serial_reader.domain_size[1], 1])
        
        elif self._dim == 3:
            self._periodic_dimensions = numpy.asarray(serial_reader.dimension)
            self._domain_size = numpy.asarray(serial_reader.domain_size)
        
        if num_ghosts is None:
            self._num_ghosts = numpy.array([0, 0, 0], dtype=numpy.int32)
        else:
            if len(num_ghosts) != 3:
                raise RuntimeError('Dimension of num_ghosts should be 3!')
            self._num_ghosts = numpy.asarray(num_ghosts, dtype=numpy.int32)
        
        self._steps = serial_reader.steps
        self._step = serial_reader.step
        
        self._lo_subdomain = numpy.array([0, 0, 0], dtype=self._domain_size.dtype)
        self._hi_subdomain = self._domain_size - numpy.array([1, 1, 1], dtype=self._domain_size.dtype)
        
        # Get the local chunk.
        
        self._lo_chunk = self._lo_subdomain
        self._hi_chunk = self._hi_subdomain
    
    
    @property
    def serial_reader(self):
        """
        Return the serial data reader.
        """
        
        return self._serial_reader
    
    
    @abc.abstractproperty
    def domain_size(self):
        """
        Return a tuple containing the full domain size of this dataset.
        """
        return self._domain_size[0:self._dim]
    
    
    @property
    def periodic_dimensions(self):
        """
        Return a tuple indicating if data is periodic in each dimension.
        """
        
        return self._periodic_dimensions[0:self._dim]
    
    
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
        self._step = step
    
    
    def getStep(self, step):
        """
        Return the time step that is currently set.
        """
        
        return self._step
    
    
    step = property(getStep, setStep)
    
    
    def setSubDomain(self, lo_and_hi):
        """
        Set the sub-domain for reading coordinates and data.
        (Not yet implemented!)
        """
        
        return
    
    
    def getSubDomain(self):
        """
        Return two tuples containing the sub-domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """
        
        return self._lo_subdomain, self._hi_subdomain
    
    
    sub_domain = property(getSubDomain, setSubDomain)
    
    
    @property
    def chunk(self):
        """
        Return two tuples containing the chunk of sub-domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """
        
        return self._lo_chunk, self._hi_chunk
    
    
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
