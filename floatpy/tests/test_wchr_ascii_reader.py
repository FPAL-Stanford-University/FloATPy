import numpy
import os
import unittest

import floatpy.readers.wchr_ascii_reader as war

class TestWchrAsciiReader(unittest.TestCase):

    def setUp(self):
        self.filename_prefix = os.path.join(os.path.dirname(__file__), 'test_wchr_ascii_reader_data/WCHR_')
        self.reader = war.WchrAsciiReader(self.filename_prefix)

        # Read in the full domain coordinates
        self.reader.readXCoord()
        self.reader.readYCoord()
        self.reader.readZCoord()


    def test_readCoordinates_chunk(self):

        # Define a chunk
        chunk = ((2,5),(4,8),(1,7))

        # Read chunked coordinates
        x, y, z = self.reader.readCoordinates(chunk)

        # Check that the chunked coordinates are equal to the corresponding full coords
        xerr = numpy.absolute(self.reader.x_c[ chunk[0][0]:chunk[0][1], chunk[1][0]:chunk[1][1], chunk[2][0]:chunk[2][1] ] - x).max()
        yerr = numpy.absolute(self.reader.y_c[ chunk[0][0]:chunk[0][1], chunk[1][0]:chunk[1][1], chunk[2][0]:chunk[2][1] ] - y).max()
        zerr = numpy.absolute(self.reader.z_c[ chunk[0][0]:chunk[0][1], chunk[1][0]:chunk[1][1], chunk[2][0]:chunk[2][1] ] - z).max()

        self.assertEqual(xerr, 0., "Incorrect chunked coordinate data reader in X")
        self.assertEqual(yerr, 0., "Incorrect chunked coordinate data reader in Y")
        self.assertEqual(zerr, 0., "Incorrect chunked coordinate data reader in Z")


    def test_readVariables_chunk(self):

        # Define a chunk
        chunk = ((2,5),(4,8),(1,7))

        # Read full data
        rho,    = self.reader.readVariables('rho', 0, chunk=None)
        u, v, w = self.reader.readVariables(('u','v','w'), 0, chunk=None)
        p,      = self.reader.readVariables('p', 0, chunk=None)

        # Read in chunked data
        rho_c,        = self.reader.readVariables('rho', 0, chunk=chunk)
        u_c, v_c, w_c = self.reader.readVariables(('u','v','w'), 0, chunk=chunk)
        p_c,          = self.reader.readVariables('p', 0, chunk=chunk)

        rerr = numpy.absolute(rho[ chunk[0][0]:chunk[0][1], chunk[1][0]:chunk[1][1], chunk[2][0]:chunk[2][1] ] - rho_c).max()
        uerr = numpy.absolute(u  [ chunk[0][0]:chunk[0][1], chunk[1][0]:chunk[1][1], chunk[2][0]:chunk[2][1] ] - u_c  ).max()
        verr = numpy.absolute(v  [ chunk[0][0]:chunk[0][1], chunk[1][0]:chunk[1][1], chunk[2][0]:chunk[2][1] ] - v_c  ).max()
        werr = numpy.absolute(w  [ chunk[0][0]:chunk[0][1], chunk[1][0]:chunk[1][1], chunk[2][0]:chunk[2][1] ] - w_c  ).max()
        perr = numpy.absolute(p  [ chunk[0][0]:chunk[0][1], chunk[1][0]:chunk[1][1], chunk[2][0]:chunk[2][1] ] - p_c  ).max()

        self.assertEqual(rerr, 0., "Incorrect chunked variable data reader for rho")
        self.assertEqual(uerr, 0., "Incorrect chunked variable data reader for u  ")
        self.assertEqual(verr, 0., "Incorrect chunked variable data reader for v  ")
        self.assertEqual(werr, 0., "Incorrect chunked variable data reader for w  ")
        self.assertEqual(perr, 0., "Incorrect chunked variable data reader for p  ")


if __name__ == '__main__':
    unittest.main()
