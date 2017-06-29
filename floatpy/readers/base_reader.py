import abc

class BaseReader(object):
    """
    Abstract base class to read data.
    
    Steps in using the data reader to post-process data in many time steps:
    1. Use the constructor BaseReader(data_directory_path) to create the object and initialize it.
    4. Call setSubDomain(lo,hi) to set the sub domain to read in (may only call once depending on the usage)
    2. Call readSummary(step) for each timestep
    5. Call readCoordinates() (may only call once depending on the usage)
    5. Call readData()
    8. Do your post-processing...

    To write a concrete class (called MyReaderImplementation, say) that derives from this, implement the following
    abstract methods and in the end of the file add the following code to register the concrete class
    BaseReader.register(MyReaderImplementation)
    """

    __metaclass__ = abc.ABCMeta
    
    
    @abc.abstractmethod
    def readSummary(self, step):
        """
        Get the meta data from the summary file in the data directory.
        Return error when data directory is not set.
        """
        return
        
    
    @abc.abstractmethod
    def setSubDomain(self, lo, hi):
        """
        Set the sub-domain for reading coordinates and data.
        Return error when summary is not read because need to check whether the subdomain is
        in the full domain.
        Return error when data is already loaded.
        Should clear data first before re-setting the sub-domain.
        """
        return
        
    
    @abc.abstractproperty
    def domain_size(self):
        """
        Return a tuple containing the full domain size of this dataset.
        """
        return

    @abc.abstractproperty
    def sub_domain(self):
        """
        Return two tuples containing the sub-domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """
        return

    @abc.abstractproperty
    def periodic_dimensions(self):
        """
        Return a tuple indicating if data is periodic in each dimension.
        """
        return

    @abc.abstractproperty
    def time(self):
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
    def readData(self, var_names, step):
        """
        Read the data of several variables in the stored sub-domain.
        Default to the full domain when the sub-domain is not set.
        """
        return
        
    
