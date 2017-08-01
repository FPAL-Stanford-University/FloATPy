import numpy
import os
import unittest

from floatpy.readers import samrai_reader

class TestSamraiDataReader2D(unittest.TestCase):

    def setUp(self):
        self.directory_name = os.path.join(os.path.dirname(__file__), 'test_data_samrai_2D')
        self.reader = samrai_reader.SamraiDataReader(self.directory_name)
        
        self.lo = (12, 36)
        self.hi = (23, 71)
        self.reader.sub_domain = (self.lo, self.hi)
        self.reader.step = 0
    
    
    def testReadCoordinatesSubdomain(self):
        
        # Read full coordinates.
        
        self.reader.sub_domain = (0,0), (self.reader.domain_size[0]-1, self.reader.domain_size[1]-1)
        
        x, y = self.reader.readCoordinates()
        
        # Read coordinates in sub-domain.
        
        self.reader.sub_domain = self.lo, self.hi
        x_s, y_s = self.reader.readCoordinates()
        
        # Check that the coordinates in sub-domain are equal to the corresponding full coordinates.
        
        x_err = numpy.absolute(x[ self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1 ] - x_s).max()
        y_err = numpy.absolute(y[ self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1 ] - y_s).max()
        
        self.assertEqual(x_err, 0., "Incorrect sub-domain coordinate data reader in x direction!")
        self.assertEqual(y_err, 0., "Incorrect sub-domain coordinate data reader in y direction!")
    
    
    def testReadDataSubdomain(self):
        
        # Read full data.
        
        self.reader.sub_domain = (0,0), (self.reader.domain_size[0]-1, self.reader.domain_size[1]-1)
        
        rho, vel, p = self.reader.readData(('density', 'velocity', 'pressure'))
        
        # Read data in sub-domain.
        
        self.reader.sub_domain = self.lo, self.hi
        rho_s, vel_s, p_s = self.reader.readData(('density', 'velocity', 'pressure'))
        
        rho_err = numpy.absolute(rho[self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1   ] - rho_s).max()
        vel_err = numpy.absolute(vel[self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1, :] - vel_s).max()
        p_err   = numpy.absolute(p  [self.lo[0]:self.hi[0]+1, self.lo[1]:self.hi[1]+1   ] - p_s  ).max()
        
        self.assertEqual(rho_err, 0., "Incorrect sub-domain variable data reader for density!")
        self.assertEqual(vel_err, 0., "Incorrect sub-domain variable data reader for velocity!")
        self.assertEqual(p_err,   0., "Incorrect sub-domain variable data reader for pressure!")


if __name__ == '__main__':
    unittest.main()
