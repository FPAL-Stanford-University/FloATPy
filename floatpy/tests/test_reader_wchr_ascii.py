import numpy
import os
import unittest

import floatpy.readers.wchr_ascii_reader as war

class TestReaderWchrAscii(unittest.TestCase):
    
    def setUp(self):
        self.filename_prefix = os.path.join(os.path.dirname(__file__), 'test_data_wchr_ascii/WCHR_')
        self.reader = war.WchrAsciiReader(self.filename_prefix)
        
        self.lo = (2 , 4, 1)
        self.hi = (5 , 8, 7)
        self.reader.sub_domain = (self.lo, self.hi)
        self.reader.step = 0
    
    
    def testReadCoordinatesChunk(self):
        
        # Read full coordinates.
        
        self.reader.sub_domain = (0,0,0), self.reader.domain_size
        x, y, z = self.reader.readCoordinates()
        
        # Read chunked coordinates.
        
        self.reader.sub_domain = self.lo, self.hi
        x_c, y_c, z_c = self.reader.readCoordinates()

        # Check that the chunked coordinates are equal to the corresponding full coords.
        
        xerr = numpy.absolute(x[ self.lo[0]:self.hi[0], self.lo[1]:self.hi[1], self.lo[2]:self.hi[2] ] - x_c).max()
        yerr = numpy.absolute(y[ self.lo[0]:self.hi[0], self.lo[1]:self.hi[1], self.lo[2]:self.hi[2] ] - y_c).max()
        zerr = numpy.absolute(z[ self.lo[0]:self.hi[0], self.lo[1]:self.hi[1], self.lo[2]:self.hi[2] ] - z_c).max()
        
        self.assertEqual(xerr, 0., "Incorrect chunked coordinate data reader in X")
        self.assertEqual(yerr, 0., "Incorrect chunked coordinate data reader in Y")
        self.assertEqual(zerr, 0., "Incorrect chunked coordinate data reader in Z")
    
    
    def testReadDataChunk(self):
        
        # Read full data.
        
        self.reader.sub_domain = (0,0,0), self.reader.domain_size
        rho,    = self.reader.readData('rho')
        u, v, w = self.reader.readData(('u','v','w'))
        p,      = self.reader.readData('p')
        
        # Read in chunked data.
        
        self.reader.sub_domain = self.lo, self.hi
        rho_c,        = self.reader.readData('rho')
        u_c, v_c, w_c = self.reader.readData(('u','v','w'))
        p_c,          = self.reader.readData('p')
        
        rerr = numpy.absolute(rho[ self.lo[0]:self.hi[0], self.lo[1]:self.hi[1], self.lo[2]:self.hi[2] ] - rho_c).max()
        uerr = numpy.absolute(u  [ self.lo[0]:self.hi[0], self.lo[1]:self.hi[1], self.lo[2]:self.hi[2] ] - u_c  ).max()
        verr = numpy.absolute(v  [ self.lo[0]:self.hi[0], self.lo[1]:self.hi[1], self.lo[2]:self.hi[2] ] - v_c  ).max()
        werr = numpy.absolute(w  [ self.lo[0]:self.hi[0], self.lo[1]:self.hi[1], self.lo[2]:self.hi[2] ] - w_c  ).max()
        perr = numpy.absolute(p  [ self.lo[0]:self.hi[0], self.lo[1]:self.hi[1], self.lo[2]:self.hi[2] ] - p_c  ).max()
        
        self.assertEqual(rerr, 0., "Incorrect chunked variable data reader for rho")
        self.assertEqual(uerr, 0., "Incorrect chunked variable data reader for u  ")
        self.assertEqual(verr, 0., "Incorrect chunked variable data reader for v  ")
        self.assertEqual(werr, 0., "Incorrect chunked variable data reader for w  ")
        self.assertEqual(perr, 0., "Incorrect chunked variable data reader for p  ")


if __name__ == '__main__':
    unittest.main()
