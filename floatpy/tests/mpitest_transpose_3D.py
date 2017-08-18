from mpi4py import MPI
import numpy
import os
import unittest

from floatpy.parallel import transpose_wrapper
from floatpy.readers import samrai_reader, parallel_reader

class TestTranspose2D(unittest.TestCase):
    
    def setUp(self):
        self.directory_name = os.path.join(os.path.dirname(__file__), 'test_data_samrai_3D')
        self.serial_reader = samrai_reader.SamraiDataReader(self.directory_name)
        self.serial_reader.step = 0
        
        self.comm = MPI.COMM_WORLD
        self.reader = parallel_reader.ParallelDataReader(MPI.COMM_WORLD, samrai_reader.SamraiDataReader(self.directory_name))
        self.reader.step = 0
    
    
    def testTransposeInX(self):
        
        # Read full data.
        
        self.serial_reader.sub_domain = (0, 0, 0), \
            (self.serial_reader.domain_size[0]-1, self.serial_reader.domain_size[1]-1, self.serial_reader.domain_size[2]-1)
        
        rho, vel = self.serial_reader.readData(('density', 'velocity'))
        
        # Read data in parallel region.
        
        rho_p, vel_p = self.reader.readData(('density', 'velocity'))
        
        tw = transpose_wrapper.TransposeWrapper(self.reader.grid_partition, direction=0, dim=3)
        lo_t, hi_t = tw.full_chunk
        
        rho_t = tw.transposeToPencil(rho_p)
        rho_err = numpy.absolute(rho[lo_t[0]:hi_t[0]+1, lo_t[1]:hi_t[1]+1, lo_t[2]:hi_t[2]+1] - rho_t).max()
        self.assertEqual(rho_err, 0.0, "Incorrect transposed data in x-direction for scalar!")
        
        vel_t = tw.transposeToPencil(vel_p)
        vel_err = numpy.absolute(vel[lo_t[0]:hi_t[0]+1, lo_t[1]:hi_t[1]+1, lo_t[2]:hi_t[2]+1, :] - vel_t).max()
        self.assertEqual(vel_err, 0.0, "Incorrect transposed data in x-direction for vector!")
    
    
    def testTransposeInY(self):
        
        # Read full data.
        
        self.serial_reader.sub_domain = (0, 0, 0), \
            (self.serial_reader.domain_size[0]-1, self.serial_reader.domain_size[1]-1, self.serial_reader.domain_size[2]-1)
        
        rho, vel = self.serial_reader.readData(('density', 'velocity'))
        
        # Read data in parallel region.
        
        rho_p, vel_p = self.reader.readData(('density', 'velocity'))
        
        tw = transpose_wrapper.TransposeWrapper(self.reader.grid_partition, direction=1, dim=3)
        lo_t, hi_t = tw.full_chunk
        
        rho_t = tw.transposeToPencil(rho_p)
        rho_err = numpy.absolute(rho[lo_t[0]:hi_t[0]+1, lo_t[1]:hi_t[1]+1, lo_t[2]:hi_t[2]+1] - rho_t).max()
        self.assertEqual(rho_err, 0.0, "Incorrect transposed data in y-direction for scalar!")
        
        vel_t = tw.transposeToPencil(vel_p)
        vel_err = numpy.absolute(vel[lo_t[0]:hi_t[0]+1, lo_t[1]:hi_t[1]+1, lo_t[2]:hi_t[2]+1, :] - vel_t).max()
        self.assertEqual(vel_err, 0.0, "Incorrect transposed data in y-direction for vector!")
    
    
    def testTransposeInZ(self):
        
        # Read full data.
        
        self.serial_reader.sub_domain = (0, 0, 0), \
            (self.serial_reader.domain_size[0]-1, self.serial_reader.domain_size[1]-1, self.serial_reader.domain_size[2]-1)
        
        rho, vel = self.serial_reader.readData(('density', 'velocity'))
        
        # Read data in parallel region.
        
        rho_p, vel_p = self.reader.readData(('density', 'velocity'))
        
        tw = transpose_wrapper.TransposeWrapper(self.reader.grid_partition, direction=2, dim=3)
        lo_t, hi_t = tw.full_chunk
        
        rho_t = tw.transposeToPencil(rho_p)
        rho_err = numpy.absolute(rho[lo_t[0]:hi_t[0]+1, lo_t[1]:hi_t[1]+1, lo_t[2]:hi_t[2]+1] - rho_t).max()
        self.assertEqual(rho_err, 0.0, "Incorrect transposed data in z-direction for scalar!")
        
        vel_t = tw.transposeToPencil(vel_p)
        vel_err = numpy.absolute(vel[lo_t[0]:hi_t[0]+1, lo_t[1]:hi_t[1]+1, lo_t[2]:hi_t[2]+1, :] - vel_t).max()
        self.assertEqual(vel_err, 0.0, "Incorrect transposed data in z-direction for vector!")


if __name__ == '__main__':
    unittest.main()
