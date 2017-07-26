from mpi4py import MPI
import numpy
import os
import unittest

import floatpy.readers.wchr_ascii_reader as war
import floatpy.readers.parallel_reader as pdr

class TestReaderWchrAscii(unittest.TestCase):
    
    def setUp(self):
        self.filename_prefix = os.path.join(os.path.dirname(__file__), 'test_data_wchr_ascii/WCHR_')
        self.serial_reader = war.WchrAsciiReader(self.filename_prefix)
        self.serial_reader.step = 0
        
        self.comm  = MPI.COMM_WORLD
        self.reader = pdr.ParallelDataReader( MPI.COMM_WORLD, war.WchrAsciiReader(self.filename_prefix), num_ghosts=(1,1,1) )
        self.reader.step = 0

        self.lo, self.hi = self.reader.interior_chunk
    
    
    def testReadCoordinatesChunk(self):
        
        # Read full coordinates.
        
        self.serial_reader.sub_domain = (0,0,0), (self.serial_reader.domain_size[0]-1, self.serial_reader.domain_size[1]-1, self.serial_reader.domain_size[2]-1)
        x, y, z = self.serial_reader.readCoordinates()

        # Read chunked coordinates.
        
        x_c, y_c, z_c = self.reader.readCoordinates()

        # Check that the chunked coordinates are equal to the corresponding full coords.
        
        xerr = numpy.absolute(x[ self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1, self.lo[2]:self.hi[2]+1 ] - x_c[self.reader.interior]).max()
        yerr = numpy.absolute(y[ self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1, self.lo[2]:self.hi[2]+1 ] - y_c[self.reader.interior]).max()
        zerr = numpy.absolute(z[ self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1, self.lo[2]:self.hi[2]+1 ] - z_c[self.reader.interior]).max()
        
        self.assertEqual(xerr, 0., "Incorrect chunked coordinate data reader in X")
        self.assertEqual(yerr, 0., "Incorrect chunked coordinate data reader in Y")
        self.assertEqual(zerr, 0., "Incorrect chunked coordinate data reader in Z")
    
    
    def testReadDataChunk(self):
        
        # Read full data.
        
        self.serial_reader.sub_domain = (0,0,0), (self.serial_reader.domain_size[0]-1, self.serial_reader.domain_size[1]-1, self.serial_reader.domain_size[2]-1)
        rho,    = self.serial_reader.readData('rho')
        u, v, w = self.serial_reader.readData(('u','v','w'))
        p,      = self.serial_reader.readData('p')
        
        # Read in chunked data.
        
        self.reader.sub_domain = self.lo, self.hi
        rho_c,        = self.reader.readData('rho')
        u_c, v_c, w_c = self.reader.readData(('u','v','w'))
        p_c,          = self.reader.readData('p')
        
        rerr = numpy.absolute(rho[ self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1, self.lo[2]:self.hi[2]+1 ] - rho_c[self.reader.interior]).max()
        uerr = numpy.absolute(u  [ self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1, self.lo[2]:self.hi[2]+1 ] - u_c  [self.reader.interior]).max()
        verr = numpy.absolute(v  [ self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1, self.lo[2]:self.hi[2]+1 ] - v_c  [self.reader.interior]).max()
        werr = numpy.absolute(w  [ self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1, self.lo[2]:self.hi[2]+1 ] - w_c  [self.reader.interior]).max()
        perr = numpy.absolute(p  [ self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1, self.lo[2]:self.hi[2]+1 ] - p_c  [self.reader.interior]).max()
        
        self.assertEqual(rerr, 0., "Incorrect chunked variable data reader for rho")
        self.assertEqual(uerr, 0., "Incorrect chunked variable data reader for u  ")
        self.assertEqual(verr, 0., "Incorrect chunked variable data reader for v  ")
        self.assertEqual(werr, 0., "Incorrect chunked variable data reader for w  ")
        self.assertEqual(perr, 0., "Incorrect chunked variable data reader for p  ")


if __name__ == '__main__':
    unittest.main()
