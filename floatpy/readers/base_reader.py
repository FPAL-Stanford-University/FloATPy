import abc

class BaseReader(object):
    """
    Abstract base class to read data.
    
    Steps in using the data reader to post-process data in many time steps:
    1. Use the constructor BaseReader(data_directory_path) to create the object and initialize it.
    2. Call getDomainSize() to get the full domain size.
    3. Call setSubDomain(lo, hi) to set the sub-domain to read in.
    4. Call readCoordinates() to get coordinates of the sub-domain.
    For each time step:
        a. Call updateSummary(step) for each timestep.
        b. Call readData().
        c. Do your post-processing...
    
    To write a concrete class (called MyReaderImplementation, say) that derives from this, implement the following
    abstract methods and in the end of the file add the following code to register the concrete class
    BaseReader.register(MyReaderImplementation).
    """
    
    __metaclass__ = abc.ABCMeta
    
    
    @abc.abstractmethod
    def updateSummary(self, step):
        """
        Update the metadata from the summary file in the data directory at new time step.
        """
        return
    
    
    @abc.abstractmethod
    def setSubDomain(self, lo, hi):
        """
        Set the sub-domain for reading coordinates and data.
        """
        return
    
    
    @abc.abstractproperty
    def getDomainSize(self):
        """
        Return a tuple containing the full domain size of this dataset.
        """
        return
    
    
    @abc.abstractproperty
    def getSubDomain(self):
        """
        Return two tuples containing the sub-domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """
        return
    
    
    @abc.abstractproperty
    def getPeriodicDimensions(self):
        """
        Return a tuple indicating if data is periodic in each dimension.
        """
        return

    @abc.abstractproperty
    def getTime(self):
        """
        Return the current simulation time.
        """
        return
    
    
    @abc.abstractmethod
    def readCoordinates(self):
        """
        Get the coordinates of the stored sub-domain.
        Default to the full domain when the sub-domain is not set.
        """
        return
    
    
    @abc.abstractmethod
    def readData(self, var_names):
        """
        Read the data of several variables in the stored sub-domain.
        Default to the full domain when the sub-domain is not set.
        """
        return

