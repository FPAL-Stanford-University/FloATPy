from base_reader import BaseReader
from mpi4py import MPI
import numpy
import os
import sys

cwd = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cwd + '/../parallel/pyt3d/')
import pyt3d

class ParallelDataReader():
    """
    Class to read data and exchange data across nodes with MPI.
    """
    
    def __init__(self, serial_reader):
        """
        Constructor of the class.
        """
        
        if not isinstance(serial_reader, BaseReader):
            raise RuntimeError("The given serial data reader is not instance of the base reader!")
        
        self._serial_reader = serial_reader
        
   