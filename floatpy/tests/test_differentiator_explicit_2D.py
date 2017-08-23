import numpy
import unittest

from floatpy.derivatives import ExplicitDifferentiator

class TestDifferentiatorExplicit(unittest.TestCase):
    
    def setUp(self):
        self.nx, self.ny = 64, 64
        self.omega = 1.

        self.order = (6, 6)

        self.dx, self.dy = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny

        self.der = ExplicitDifferentiator((self.dx, self.dy), self.order, 2, 'F')
        
        self.x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)
        self.y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)

        self.x, self.y = numpy.meshgrid(self.x, self.y, indexing='ij')
        self.x = numpy.asfortranarray(self.x)
        self.y = numpy.asfortranarray(self.y)

        self.f = numpy.sin(self.omega*self.x) * numpy.cos(self.omega*self.y)

        self.dfdx_exact =  self.omega * numpy.cos(self.omega*self.x) * numpy.cos(self.omega*self.y)
        self.dfdy_exact = -self.omega * numpy.sin(self.omega*self.x) * numpy.sin(self.omega*self.y)

        self.d2fdx2_exact = -self.omega**2 * numpy.sin(self.omega*self.x) * numpy.cos(self.omega*self.y)
        self.d2fdy2_exact = -self.omega**2 * numpy.sin(self.omega*self.x) * numpy.cos(self.omega*self.y)


    def testFirstDerivativeX(self):
        """
        Test the first derivative in the X direction.
        """

        dfdx = numpy.empty( (self.nx+1, self.ny+1), dtype=numpy.float64, order='F' )
        self.der.ddx(self.f, dfdx, use_one_sided=True)

        error = numpy.array([ numpy.absolute(self.dfdx_exact - dfdx).max() ])
        self.assertLess(error[0], 1.0e-6, "Incorrect first derivative in X direction!")


    def testFirstDerivativeY(self):
        """
        Test the first derivative in the Y direction.
        """

        dfdy = numpy.empty( (self.nx+1, self.ny+1), dtype=numpy.float64, order='F' )
        self.der.ddy(self.f, dfdy, use_one_sided=True)

        error = numpy.array([ numpy.absolute(self.dfdy_exact - dfdy).max() ])
        self.assertLess(error[0], 1.0e-6, "Incorrect first derivative in Y direction!")


    def testSecondDerivativeX(self):
        """
        Test the second derivative in the X direction.
        """

        d2fdx2 = numpy.empty( (self.nx+1, self.ny+1), dtype=numpy.float64, order='F' )
        self.der.d2dx2(self.f, d2fdx2, use_one_sided=True)

        error = numpy.array([ numpy.absolute(self.d2fdx2_exact - d2fdx2).max() ])
        self.assertLess(error[0], 1.0e-6, "Incorrect second derivative in the X direction!")


    def testSecondDerivativeY(self):
        """
        Test the second derivative in the Y direction.
        """

        d2fdy2 = numpy.empty( (self.nx+1, self.ny+1), dtype=numpy.float64, order='F' )
        self.der.d2dy2(self.f, d2fdy2, use_one_sided=True)

        error = numpy.array([ numpy.absolute(self.d2fdy2_exact - d2fdy2).max() ])
        self.assertLess(error[0], 1.0e-6, "Incorrect second derivative in the Y direction!")


    def testGradient(self):
        """
        Test the gradient function.
        """

        dfdx, dfdy = self.der.gradient(self.f, use_one_sided=True)

        error = numpy.zeros(2)
        error[0] = numpy.absolute(self.dfdx_exact - dfdx).max()
        error[1] = numpy.absolute(self.dfdy_exact - dfdy).max()

        for i in range(2):
            self.assertLess(error[i], 1.0e-6, "Incorrect gradient!")


    def testDivergence(self):
        """
        Test the divergence function.
        """

        df = numpy.concatenate((self.dfdx_exact[..., numpy.newaxis], \
            self.dfdy_exact[..., numpy.newaxis]), axis=2)

        divergence = self.der.divergence(df, use_one_sided=True)

        error = numpy.zeros(1)
        error[0] = numpy.absolute(self.d2fdx2_exact + self.d2fdy2_exact - divergence).max()

        self.assertLess(error[0], 1.0e-6, "Incorrect divergence!")


    def testCurl(self):
        """
        Test the curl function.
        """

        df = numpy.concatenate((self.dfdx_exact[..., numpy.newaxis], \
            self.dfdy_exact[..., numpy.newaxis]), axis=2)

        curl = self.der.curl(df, use_one_sided=True)

        error = numpy.zeros(1)
        error[0] = numpy.absolute(curl).max()

        self.assertLess(error[0], 1.0e-6, "Incorrect curl!")


    def testLaplacian(self):
        """
        Test the laplacian function.
        """

        laplacian = self.der.laplacian(self.f, use_one_sided=True)

        error = numpy.zeros(1)
        error[0] = numpy.absolute(self.d2fdx2_exact + self.d2fdy2_exact - laplacian).max()

        self.assertLess(error[0], 5.0e-5, "Incorrect laplacian!")


if __name__ == '__main__':
    unittest.main()
