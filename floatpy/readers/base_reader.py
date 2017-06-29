class BaseReader:
    """
    Abstract class to reader data.
    
    Steps in using the data reader to post-process data in many time steps:
    1. Call setDataDirectoryPath()
    2. Call readSummary() (may only call once depending on the data format)
    3. Call getFullDomainSize() (may only call once depending on the usage)
    4. Call setSubDomain() (may only call once depending on the usage)
    5. Call getCoordinates() (may only call once depending on the usage)
    5. Call readData()
    6. Call getData()
    8. Do some post-processing...
    7. Call clearData()/clear() (call clear if user wants to read summary from another directory)
    """
    
    def setDataDirectoryPath(self, data_directory_path):
        """
        Set the absolute path of the data directory.
        Return error when data is already loaded.
        Should clear data first before re-setting the path.
        """
        
        raise RuntimeError("Not yet implemented")
    
    
    def getDataDirectoryPath(self):
        """
        Get the stored absolute path of the data directory.
        """
        
        raise RuntimeError("Not yet implemented")
    
    
    def readSummary(self):
        """
        Get the meta data from the summary file in the data directory.
        Return error when data directory is not set.
        Return error when data is already loaded.
        """
        
        raise RuntimeError("Not yet implemented")
    
    
    def readSummary(self, summary_file_path):
        """
        Get the meta data from the summary file in the given path.
        Return error when data is already loaded.
        """
        
        raise RuntimeError("Not yet implemented")
    
    
    def getFullDomainSize(self):
        """
        Get the full domain size.
        Return error when summary is not read.
        """
        
        raise RuntimeError("Not yet implemented")
    
    
    def setSubDomain(self, lo, hi):
        """
        Set the sub-domain for reading coordinates and data.
        Return error when summary is not read because need to check whether the subdomain is
        in the full domain.
        Return error when data is already loaded.
        Should clear data first before re-setting the sub-domain.
        """
        
        raise RuntimeError("Not yet implemented")
    
    
    def getCoordinates(self):
        """
        Get the coordinates of the stored sub-domain.
        Return error when the sub-domain is not set.
        """
        
        raise RuntimeError("Not yet implemented")
    
    
    def readData(self, var_names, periodic_dimension):
        """
        Read the data of several variables in the stored sub-domain.
        Return error when the sub-domain is not set.
        Return error when data is already loaded.
        """
        
        raise RuntimeError("Not yet implemented")
    
    
    def getData(self, var_name):
        """
        Get the loaded data of a variable.
        Return error when data is not loaded yet.
        """
        
        raise RuntimeError("Not yet implemented")
    
    
    def clearData(self):
        """
        Clear any loaded data.
        """
        
        raise RuntimeError("Not yet implemented")
    
    
    def clear(self):
        """
        Clear all data in the class.
        """
        
        raise RuntimeError("Not yet implemented")
