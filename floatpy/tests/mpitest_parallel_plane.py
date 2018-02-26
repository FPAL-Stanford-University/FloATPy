from mpi4py import MPI
import numpy
import os
import unittest

import floatpy.readers.miranda_reader as mir
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.parallel_plane as pp

class TestReaderParallel(unittest.TestCase):
    
    def setUp(self):
        self.filename_prefix = os.path.join(os.path.dirname(__file__), 'test_data_miranda/plot.mir')
        self.serial_reader = mir.MirandaReader(self.filename_prefix, periodic_dimensions=(False,True,True))
        self.serial_reader.step = 0
        
        self.comm = MPI.COMM_WORLD
        self.num_ghosts = (1, 1, 1)
        self.reader = pdr.ParallelDataReader( MPI.COMM_WORLD, mir.MirandaReader(self.filename_prefix, periodic_dimensions=(False,True,True)), num_ghosts=self.num_ghosts )
        self.reader.step = 0
        
        self.lo, self.hi = self.reader.interior_chunk

        self.direction = 2
        self.index = 0
        self.pp = pp.ParallelPlane(self.reader.grid_partition, 2, 0)
    
    
    def testReadCoordinatesPlane(self):
        
        # Read full coordinates.
        
        self.serial_reader.sub_domain = (0,0,0), (self.serial_reader.domain_size[0]-1, self.serial_reader.domain_size[1]-1, self.serial_reader.domain_size[2]-1)
        x, y, z = self.serial_reader.readCoordinates()
        
        # Read chunked coordinates.
        
        x_c, y_c, z_c = self.reader.readCoordinates()

        # Get a plane of the coordinates
        has_plane_x, x_p = self.pp.get_plane(x_c[self.reader.interior])
        has_plane_y, y_p = self.pp.get_plane(y_c[self.reader.interior])
        has_plane_z, z_p = self.pp.get_plane(z_c[self.reader.interior])
        
        # Check that the plane is gathered correctly
        if has_plane_x:
            xerr = numpy.absolute(x[ :, :, self.index ] - x_p).max()
            self.assertEqual(xerr, 0., "Incorrect chunked coordinate data reader in X")
        if has_plane_y:
            yerr = numpy.absolute(y[ :, :, self.index ] - y_p).max()
            self.assertEqual(yerr, 0., "Incorrect chunked coordinate data reader in Y")
        if has_plane_z:
            zerr = numpy.absolute(z[ :, :, self.index ] - z_p).max()
            self.assertEqual(zerr, 0., "Incorrect chunked coordinate data reader in Z")
    
    
    def testReadDataPlane(self):
        
        # Read full data.
        
        self.serial_reader.sub_domain = (0,0,0), (self.serial_reader.domain_size[0]-1, self.serial_reader.domain_size[1]-1, self.serial_reader.domain_size[2]-1)
        rho,    = self.serial_reader.readData('density')
        
        # Read in chunked data.
        
        rho_c,        = self.reader.readData('density')
        
        # Get a plane of data
        has_plane, rho_p = self.pp.get_plane(rho_c[self.reader.interior])
        
        if has_plane:
            rerr = numpy.absolute(rho[ :, :, self.index ] - rho_p).max()
            self.assertEqual(rerr, 0., "Incorrect chunked variable data reader for rho")


if __name__ == '__main__':
    unittest.main()
