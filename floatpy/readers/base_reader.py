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
    
    # def __init__(self, data_directory_path):
    #     """
    #     Set the absolute path of the data directory.
    #     Return error when data is already loaded.
    #     Should clear data first before re-setting the path.
    #     """
    #     
    
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
    def domainSize(self):
        """
        Return the full domain size of this dataset.
        """
        return

    @abc.abstractproperty
    def getSubDomain(self):
        """
        Return the sub-domain used in this reader.
        """
        return

    @abc.abstractproperty
    def getPeriodicDimensions(self):
        """
        Return an array indicating if data is periodic in each dimension.
        """
        return

    @abc.abstractmethod
    def readCoordinates(self):
        """
        Get the coordinates of the stored sub-domain.
        Return error when the sub-domain is not set.
        """
        return
        
    
    @abc.abstractmethod
    def readData(self, var_names, step):
        """
        Read the data of several variables in the stored sub-domain.
        Return error when the sub-domain is not set.
        """
        return
        
    
    
    # def __del__(self):
    #     """
    #     Clear all data in the class.
    #     """
    #     
