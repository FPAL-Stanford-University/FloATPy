from mpi4py import MPI
import numpy
import unittest

from floatpy.parallel import _t3dmod
from floatpy.derivatives.compact import CompactDerivative


class TestDerivativesCompact(unittest.TestCase):
    
    def setUp(self):
        self.nx, self.ny, self.nz = 64, 64, 64
        self.omega = 1.

        self.comm = MPI.COMM_WORLD
        self.fcomm = self.comm.py2f()
        self.periodic = numpy.array([True, True, True])
        self.order = (10,10,10)

        self.dx, self.dy, self.dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz

        self.grid_partition = _t3dmod.t3d(self.fcomm, self.nx, self.ny, self.nz, self.periodic )

        self.chunk_3d_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.grid_partition.get_sz3d(self.chunk_3d_size)
        self.grid_partition.get_st3d(self.chunk_3d_lo)
        self.grid_partition.get_en3d(self.chunk_3d_hi)
        self.chunk_3d_lo = self.chunk_3d_lo - 1 # Convert to 0 based indexing
        self.chunk_3d_hi = self.chunk_3d_hi - 1 # Convert to 0 based indexing

        self.der = CompactDerivative(self.grid_partition, (self.dx, self.dy, self.dz), self.order, self.periodic)
        
        self.x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[self.chunk_3d_lo[0]:self.chunk_3d_hi[0]+1]
        self.y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[self.chunk_3d_lo[1]:self.chunk_3d_hi[1]+1]
        self.z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[self.chunk_3d_lo[2]:self.chunk_3d_hi[2]+1]

        self.x, self.y, self.z = numpy.meshgrid(self.x, self.y, self.z, indexing='ij')
        self.x = numpy.asfortranarray(self.x)
        self.y = numpy.asfortranarray(self.y)
        self.z = numpy.asfortranarray(self.z)

        self.f = numpy.sin(self.omega*self.x) * numpy.cos(self.omega*self.y) * numpy.cos(self.omega*self.z)

        self.dfdx_exact =  self.omega * numpy.cos(self.omega*self.x) * numpy.cos(self.omega*self.y) * numpy.cos(self.omega*self.z)
        self.dfdy_exact = -self.omega * numpy.sin(self.omega*self.x) * numpy.sin(self.omega*self.y) * numpy.cos(self.omega*self.z)
        self.dfdz_exact = -self.omega * numpy.sin(self.omega*self.x) * numpy.cos(self.omega*self.y) * numpy.sin(self.omega*self.z)

        self.d2fdx2_exact = -self.omega**2 * numpy.sin(self.omega*self.x) * numpy.cos(self.omega*self.y) * numpy.cos(self.omega*self.z)
        self.d2fdy2_exact = -self.omega**2 * numpy.sin(self.omega*self.x) * numpy.cos(self.omega*self.y) * numpy.cos(self.omega*self.z)
        self.d2fdz2_exact = -self.omega**2 * numpy.sin(self.omega*self.x) * numpy.cos(self.omega*self.y) * numpy.cos(self.omega*self.z)


    def testFirstDerivativeX(self):
        """
        Test the first derivative in the X direction.
        """

        dfdx = numpy.empty( self.chunk_3d_size, dtype=numpy.float64, order='F' )
        self.der.ddx(self.f, dfdx)

        myerror = numpy.array([ numpy.absolute(self.dfdx_exact - dfdx).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 5.0e-14, "Incorrect periodic first derivative in X!")


    def testFirstDerivativeY(self):
        """
        Test the first derivative in the Y direction.
        """

        dfdy = numpy.empty( self.chunk_3d_size, dtype=numpy.float64, order='F' )
        self.der.ddy(self.f, dfdy)

        myerror = numpy.array([ numpy.absolute(self.dfdy_exact - dfdy).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 5.0e-14, "Incorrect periodic first derivative in Y!")


    def testFirstDerivativeZ(self):
        """
        Test the first derivative in the Z direction.
        """

        dfdz = numpy.empty( self.chunk_3d_size, dtype=numpy.float64, order='F' )
        self.der.ddz(self.f, dfdz)

        myerror = numpy.array([ numpy.absolute(self.dfdz_exact - dfdz).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 5.0e-14, "Incorrect periodic first derivative in Z!")
    

    def testSecondDerivativeX(self):
        """
        Test the second derivative in the X direction.
        """

        d2fdx2 = numpy.empty( self.chunk_3d_size, dtype=numpy.float64, order='F' )
        self.der.d2dx2(self.f, d2fdx2)

        myerror = numpy.array([ numpy.absolute(self.d2fdx2_exact - d2fdx2).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 5.0e-12, "Incorrect periodic second derivative in X!")


    def testSecondDerivativeY(self):
        """
        Test the second derivative in the Y direction.
        """

        d2fdy2 = numpy.empty( self.chunk_3d_size, dtype=numpy.float64, order='F' )
        self.der.d2dy2(self.f, d2fdy2)

        myerror = numpy.array([ numpy.absolute(self.d2fdy2_exact - d2fdy2).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 5.0e-12, "Incorrect periodic second derivative in X!")


    def testSecondDerivativeZ(self):
        """
        Test the second derivative in the Z direction.
        """

        d2fdz2 = numpy.empty( self.chunk_3d_size, dtype=numpy.float64, order='F' )
        self.der.d2dz2(self.f, d2fdz2)

        myerror = numpy.array([ numpy.absolute(self.d2fdz2_exact - d2fdz2).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 5.0e-12, "Incorrect periodic second derivative in X!")


    def testGradient(self):
        """
        Test the gradient function
        """

        dfdx, dfdy, dfdz = self.der.gradient(self.f)

        myerror = numpy.zeros(3)
        myerror[0] = numpy.absolute(self.dfdx_exact - dfdx).max()
        myerror[1] = numpy.absolute(self.dfdy_exact - dfdy).max()
        myerror[2] = numpy.absolute(self.dfdz_exact - dfdz).max()

        error = numpy.zeros(3)
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        for i in range(3):
            self.assertLess(error[i], 5.0e-14, "Incorrect gradient!")


    def testDivergence(self):
        """
        Test the divergence function
        """

        laplacian = self.der.divergence(self.dfdx_exact, self.dfdy_exact, self.dfdz_exact)

        myerror = numpy.zeros(1)
        myerror[0] = numpy.absolute(self.d2fdx2_exact + self.d2fdy2_exact + self.d2fdx2_exact - laplacian).max()

        error = numpy.zeros(1)
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 5.0e-14, "Incorrect divergence!")


    def testCurl(self):
        """
        Test the curl function
        """

        curl = self.der.curl(self.dfdx_exact, self.dfdy_exact, self.dfdz_exact)

        myerror = numpy.zeros(3)
        myerror[0] = numpy.absolute(curl[:,:,:,0]).max()
        myerror[1] = numpy.absolute(curl[:,:,:,1]).max()
        myerror[2] = numpy.absolute(curl[:,:,:,2]).max()

        error = numpy.zeros(3)
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        for i in range(3):
            self.assertLess(error[i], 5.0e-14, "Incorrect curl!")


    def testLaplacian(self):
        """
        Test the laplacian function
        """

        laplacian = self.der.laplacian(self.f)

        myerror = numpy.zeros(1)
        myerror[0] = numpy.absolute(self.d2fdx2_exact + self.d2fdy2_exact + self.d2fdx2_exact - laplacian).max()

        error = numpy.zeros(1)
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 5.0e-12, "Incorrect laplacian!")


if __name__ == '__main__':
    unittest.main()
