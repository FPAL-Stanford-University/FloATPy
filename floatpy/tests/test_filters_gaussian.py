import numpy
import unittest

import floatpy.filters.pygaussian as pygaussian

def getTransferFunction(k):
    agf = 3565. /  10368.
    bgf = 3091. /  12960.
    cgf = 1997. /  25920.
    dgf =  149. /  12960.
    egf =  107. / 103680.

    T = (agf + 2.* bgf*numpy.cos(k) + 2.*cgf*numpy.cos(2.*k) +2.*dgf*numpy.cos(3.*k) + 2.*egf*numpy.cos(4.*k) )
    return T

class TestFiltersGaussian(unittest.TestCase):
    
    def setUp(self):
        self.nx, self.ny, self.nz = 32, 32, 32
        self.omega = 12.

    def testFilterPeriodicX(self):
        """
        Test the periodic first derivatives in the X direction.
        """

        dx, dy, dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz
        k_norm = self.omega * dx               
        TF = getTransferFunction(k_norm)

        x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[:-1]
        y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[:-1]
        z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[:-1]

        x, y, z = numpy.meshgrid(x, y, z, indexing='ij')
        x = numpy.asfortranarray(x)
        y = numpy.asfortranarray(y)
        z = numpy.asfortranarray(z)

        f = numpy.sin(self.omega*x) * numpy.cos(self.omega*y) * numpy.cos(self.omega*z)
        f_tilde_exact = TF * f

        fil = pygaussian.gaussianstuff.gaussian(self.nx, True)     

        f_tilde = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        fil.filter1(f, f_tilde, self.ny, self.nz)

        error = numpy.absolute(f_tilde_exact - f_tilde).max()
        self.assertLess(error, 5.0e-15, "Incorrect gaussian filter in X direction!")


    def testFilterPeriodicY(self):
        """
        Test the periodic first derivatives in the Y direction.
        """

        dx, dy, dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz
        k_norm = self.omega * dy               
        TF = getTransferFunction(k_norm)

        x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[:-1]
        y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[:-1]
        z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[:-1]

        x, y, z = numpy.meshgrid(x, y, z, indexing='ij')
        x = numpy.asfortranarray(x)
        y = numpy.asfortranarray(y)
        z = numpy.asfortranarray(z)

        f = numpy.sin(self.omega*x) * numpy.cos(self.omega*y) * numpy.cos(self.omega*z)
        f_tilde_exact = TF * f

        fil = pygaussian.gaussianstuff.gaussian(self.ny, True)     

        f_tilde = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        fil.filter2(f, f_tilde, self.nx, self.nz)

        error = numpy.absolute(f_tilde_exact - f_tilde).max()
        self.assertLess(error, 5.0e-15, "Incorrect gaussian filter in Y direction!")
    
    
    def testFilterPeriodicZ(self):
        """
        Test the periodic first derivatives in the Z direction.
        """

        dx, dy, dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz
        k_norm = self.omega * dz               
        TF = getTransferFunction(k_norm)

        x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[:-1]
        y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[:-1]
        z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[:-1]

        x, y, z = numpy.meshgrid(x, y, z, indexing='ij')
        x = numpy.asfortranarray(x)
        y = numpy.asfortranarray(y)
        z = numpy.asfortranarray(z)

        f = numpy.sin(self.omega*x) * numpy.cos(self.omega*y) * numpy.cos(self.omega*z)
        f_tilde_exact = TF * f

        fil = pygaussian.gaussianstuff.gaussian(self.nz, True)     

        f_tilde = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        fil.filter3(f, f_tilde, self.nx, self.ny)

        error = numpy.absolute(f_tilde_exact - f_tilde).max()
        self.assertLess(error, 5.0e-15, "Incorrect gaussian filter in Z direction!")
    
    
if __name__ == '__main__':
    unittest.main()
