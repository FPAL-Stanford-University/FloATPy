import numpy
import unittest

import floatpy.filters.pycf90 as pycf90

def getTransferFunction(k):
    alpha90 = 6.6624e-1
    beta90  = 1.6688e-1
    a90     = 9.9965e-1
    b90     = 6.6652e-1 # Alreday divided by factor of 2
    c90     = 1.6674e-1 # Alreday divided by factor of 2
    d90     = 4.0e-5    # Alreday divided by factor of 2
    e90     = -5.0e-6   # Alreday divided by factor of 2
    T = (a90 + 2.* b90*numpy.cos(k) + 2.*c90*numpy.cos(2.*k) + 2.*d90*numpy.cos(3.*k) + 2.*e90*numpy.cos(4.*k) ) \
      / (1. + 2.*alpha90*numpy.cos(k) + 2.*beta90*numpy.cos(2.*k) )
    return T

class TestFiltersCompact(unittest.TestCase):
    
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

        fil = pycf90.cf90stuff.cf90(self.nx, True)     

        f_tilde = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        fil.filter1(f, f_tilde, self.ny, self.nz)

        error = numpy.absolute(f_tilde_exact - f_tilde).max()
        self.assertLess(error, 5.0e-14, "Incorrect compact filter in X direction!")


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

        fil = pycf90.cf90stuff.cf90(self.ny, True)     

        f_tilde = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        fil.filter2(f, f_tilde, self.nx, self.nz)

        error = numpy.absolute(f_tilde_exact - f_tilde).max()
        self.assertLess(error, 5.0e-14, "Incorrect compact filter in Y direction!")
    
    
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

        fil = pycf90.cf90stuff.cf90(self.nz, True)     

        f_tilde = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        fil.filter3(f, f_tilde, self.nx, self.ny)

        error = numpy.absolute(f_tilde_exact - f_tilde).max()
        self.assertLess(error, 5.0e-14, "Incorrect compact filter in Z direction!")
    
    
if __name__ == '__main__':
    unittest.main()
