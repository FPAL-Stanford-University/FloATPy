from mpi4py import MPI
import numpy
import os
import unittest

import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr

class TestReaderParallel(unittest.TestCase):
    
    def setUp(self):
        self.filename = os.path.join(os.path.dirname(__file__), 'test_data_padeops/taylorgreen.h5')
        self.serial_reader = por.PadeopsReader(self.filename, periodic_dimensions=(True,True,True))
        self.serial_reader.step = 0
        
        self.comm = MPI.COMM_WORLD
        self.num_ghosts = (1, 1, 1)
        self.reader = pdr.ParallelDataReader( MPI.COMM_WORLD, por.PadeopsReader(self.filename, periodic_dimensions=(True,True,True)), num_ghosts=self.num_ghosts )
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
    
    
    def testReadCoordinatesChunkWithCommunication(self):
        
        # Read full coordinates.
        
        self.serial_reader.sub_domain = (0,0,0), (self.serial_reader.domain_size[0]-1, self.serial_reader.domain_size[1]-1, self.serial_reader.domain_size[2]-1)
        x, y, z = self.serial_reader.readCoordinates()
        
        # Read chunked coordinates.
        
        x_c, y_c, z_c = self.reader.readCoordinates(communicate=True)
        
        # Get the indices with ghost cells.
        
        indices_ghost_x = range(self.lo[0]-self.num_ghosts[0], self.hi[0]+self.num_ghosts[0]+1)
        indices_ghost_y = range(self.lo[1]-self.num_ghosts[1], self.hi[1]+self.num_ghosts[1]+1)
        indices_ghost_z = range(self.lo[2]-self.num_ghosts[2], self.hi[2]+self.num_ghosts[2]+1)
        
        # Create the correct coordinates in domain with ghost cells.
        
        x_ghost = x.take(indices_ghost_x, axis=0, mode='wrap')
        x_ghost = x_ghost.take(indices_ghost_y, axis=1, mode='wrap')
        x_ghost = x_ghost.take(indices_ghost_z, axis=2, mode='wrap')
        
        y_ghost = y.take(indices_ghost_x, axis=0, mode='wrap')
        y_ghost = y_ghost.take(indices_ghost_y, axis=1, mode='wrap')
        y_ghost = y_ghost.take(indices_ghost_z, axis=2, mode='wrap')
        
        z_ghost = z.take(indices_ghost_x, axis=0, mode='wrap')
        z_ghost = z_ghost.take(indices_ghost_y, axis=1, mode='wrap')
        z_ghost = z_ghost.take(indices_ghost_z, axis=2, mode='wrap')
        
        # Check that the chunked coordinates are equal to the corresponding full coords.
        
        xerr = numpy.absolute(x_ghost - x_c).max()
        yerr = numpy.absolute(y_ghost - y_c).max()
        zerr = numpy.absolute(z_ghost - z_c).max()
        
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


    def testReadDataChunkWithCommunication(self):
        
        # Read full data.
        
        self.serial_reader.sub_domain = (0,0,0), (self.serial_reader.domain_size[0]-1, self.serial_reader.domain_size[1]-1, self.serial_reader.domain_size[2]-1)
        rho,    = self.serial_reader.readData('rho')
        u, v, w = self.serial_reader.readData(('u','v','w'))
        p,      = self.serial_reader.readData('p')
        
        # Read in chunked data.
        
        rho_c,        = self.reader.readData('rho', communicate=True)
        u_c, v_c, w_c = self.reader.readData(('u','v','w'), communicate=True)
        p_c,          = self.reader.readData('p', communicate=True)
        
        # Get the indices with ghost cells.
        
        indices_ghost_x = range(self.lo[0]-self.num_ghosts[0], self.hi[0]+self.num_ghosts[0]+1)
        indices_ghost_y = range(self.lo[1]-self.num_ghosts[1], self.hi[1]+self.num_ghosts[1]+1)
        indices_ghost_z = range(self.lo[2]-self.num_ghosts[2], self.hi[2]+self.num_ghosts[2]+1)
        
        # Create the correct data in domain with ghost cells.
        
        rho_ghost = rho.take(indices_ghost_x, axis=0, mode='wrap')
        rho_ghost = rho_ghost.take(indices_ghost_y, axis=1, mode='wrap')
        rho_ghost = rho_ghost.take(indices_ghost_z, axis=2, mode='wrap')
        
        u_ghost = u.take(indices_ghost_x, axis=0, mode='wrap')
        u_ghost = u_ghost.take(indices_ghost_y, axis=1, mode='wrap')
        u_ghost = u_ghost.take(indices_ghost_z, axis=2, mode='wrap')
        
        v_ghost = v.take(indices_ghost_x, axis=0, mode='wrap')
        v_ghost = v_ghost.take(indices_ghost_y, axis=1, mode='wrap')
        v_ghost = v_ghost.take(indices_ghost_z, axis=2, mode='wrap')
        
        w_ghost = w.take(indices_ghost_x, axis=0, mode='wrap')
        w_ghost = w_ghost.take(indices_ghost_y, axis=1, mode='wrap')
        w_ghost = w_ghost.take(indices_ghost_z, axis=2, mode='wrap')
        
        p_ghost = p.take(indices_ghost_x, axis=0, mode='wrap')
        p_ghost = p_ghost.take(indices_ghost_y, axis=1, mode='wrap')
        p_ghost = p_ghost.take(indices_ghost_z, axis=2, mode='wrap')
        
        rerr = numpy.absolute(rho_ghost - rho_c).max()
        uerr = numpy.absolute(u_ghost   - u_c  ).max()
        verr = numpy.absolute(v_ghost   - v_c  ).max()
        werr = numpy.absolute(w_ghost   - w_c  ).max()
        perr = numpy.absolute(p_ghost   - p_c  ).max()
        
        self.assertEqual(rerr, 0., "Incorrect chunked variable data reader for rho")
        self.assertEqual(uerr, 0., "Incorrect chunked variable data reader for u  ")
        self.assertEqual(verr, 0., "Incorrect chunked variable data reader for v  ")
        self.assertEqual(werr, 0., "Incorrect chunked variable data reader for w  ")
        self.assertEqual(perr, 0., "Incorrect chunked variable data reader for p  ")


if __name__ == '__main__':
    unittest.main()
