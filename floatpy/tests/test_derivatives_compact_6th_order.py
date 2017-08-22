import numpy
import unittest

import floatpy.derivatives.compact.pycd06 as pycd06

class TestDerivativesCompact6th(unittest.TestCase):
    
    def setUp(self):
        self.nx, self.ny, self.nz = 64, 64, 64
        self.omega = 1.

    def testDerivativePeriodicX(self):
        """
        Test the periodic first derivatives in the X direction.
        """

        dx, dy, dz = 0., 0., 0.
        periodic = True

        if periodic:
            dx, dy, dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[:-1]
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[:-1]
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[:-1]
        else:
            dx, dy, dz = 2.*numpy.pi / (self.nx-1), 2.*numpy.pi / (self.ny-1), 2.*numpy.pi / (self.nz-1)
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx)
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny)
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz)

        x, y, z = numpy.meshgrid(x, y, z, indexing='ij')
        x = numpy.asfortranarray(x)
        y = numpy.asfortranarray(y)
        z = numpy.asfortranarray(z)

        f = numpy.sin(self.omega*x) * numpy.cos(self.omega*y) * numpy.cos(self.omega*z)
        dfdx_exact = self.omega * numpy.cos(self.omega*x) * numpy.cos(self.omega*y) * numpy.cos(self.omega*z)

        der = pycd06.cd06stuff.cd06(self.nx, dx, periodic, 0, 0)     

        dfdx = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        der.dd1(f, dfdx, self.ny, self.nz)

        error = numpy.absolute(dfdx_exact - dfdx).max()
        self.assertLess(error, 1.0e-9, "Incorrect derivative in first direction for second order finite difference!")
    
    def testDerivativeNonperiodicX(self):
        """
        Test the non-periodic first derivatives in the X direction.
        """

        dx, dy, dz = 0., 0., 0.
        periodic = False

        if periodic:
            dx, dy, dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[:-1]
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[:-1]
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[:-1]
        else:
            dx, dy, dz = 2.*numpy.pi / (self.nx-1), 2.*numpy.pi / (self.ny-1), 2.*numpy.pi / (self.nz-1)
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx)
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny)
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz)

        x, y, z = numpy.meshgrid(x, y, z, indexing='ij')
        x = numpy.asfortranarray(x)
        y = numpy.asfortranarray(y)
        z = numpy.asfortranarray(z)

        f = numpy.sin(self.omega*x) * numpy.cos(self.omega*y) * numpy.cos(self.omega*z)
        dfdx_exact = self.omega * numpy.cos(self.omega*x) * numpy.cos(self.omega*y) * numpy.cos(self.omega*z)

        der = pycd06.cd06stuff.cd06(self.nx, dx, periodic, 0, 0)     

        dfdx = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        der.dd1(f, dfdx, self.ny, self.nz)

        error = numpy.absolute(dfdx_exact - dfdx).max()
        self.assertLess(error, 5.0e-5, "Incorrect derivative in first direction for second order finite difference!")
    
    def testDerivativePeriodicY(self):
        """
        Test the periodic first derivatives in the Y direction.
        """

        dx, dy, dz = 0., 0., 0.
        periodic = True

        if periodic:
            dx, dy, dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[:-1]
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[:-1]
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[:-1]
        else:
            dx, dy, dz = 2.*numpy.pi / (self.nx-1), 2.*numpy.pi / (self.ny-1), 2.*numpy.pi / (self.nz-1)
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx)
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny)
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz)

        x, y, z = numpy.meshgrid(x, y, z, indexing='ij')
        x = numpy.asfortranarray(x)
        y = numpy.asfortranarray(y)
        z = numpy.asfortranarray(z)

        f = numpy.sin(self.omega*x) * numpy.cos(self.omega*y) * numpy.cos(self.omega*z)
        dfdy_exact = -self.omega * numpy.sin(self.omega*x) * numpy.sin(self.omega*y) * numpy.cos(self.omega*z)

        der = pycd06.cd06stuff.cd06(self.ny, dy, periodic, 0, 0)     

        dfdy = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        der.dd2(f, dfdy, self.nx, self.nz)

        error = numpy.absolute(dfdy_exact - dfdy).max()
        self.assertLess(error, 1.0e-9, "Incorrect derivative in first direction for second order finite difference!")
    
    def testDerivativeNonperiodicY(self):
        """
        Test the non-periodic first derivatives in the Y direction.
        """

        dx, dy, dz = 0., 0., 0.
        periodic = False

        if periodic:
            dx, dy, dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[:-1]
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[:-1]
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[:-1]
        else:
            dx, dy, dz = 2.*numpy.pi / (self.nx-1), 2.*numpy.pi / (self.ny-1), 2.*numpy.pi / (self.nz-1)
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx)
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny)
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz)

        x, y, z = numpy.meshgrid(x, y, z, indexing='ij')
        x = numpy.asfortranarray(x)
        y = numpy.asfortranarray(y)
        z = numpy.asfortranarray(z)

        f = numpy.sin(self.omega*x) * numpy.cos(self.omega*y) * numpy.cos(self.omega*z)
        dfdy_exact = -self.omega * numpy.sin(self.omega*x) * numpy.sin(self.omega*y) * numpy.cos(self.omega*z)

        der = pycd06.cd06stuff.cd06(self.ny, dy, periodic, 0, 0)     

        dfdy = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        der.dd2(f, dfdy, self.nx, self.nz)

        error = numpy.absolute(dfdy_exact - dfdy).max()
        self.assertLess(error, 5.0e-5, "Incorrect derivative in first direction for second order finite difference!")
    
    def testDerivativePeriodicZ(self):
        """
        Test the periodic first derivatives in the Z direction.
        """

        dx, dy, dz = 0., 0., 0.
        periodic = True

        if periodic:
            dx, dy, dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[:-1]
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[:-1]
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[:-1]
        else:
            dx, dy, dz = 2.*numpy.pi / (self.nx-1), 2.*numpy.pi / (self.ny-1), 2.*numpy.pi / (self.nz-1)
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx)
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny)
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz)

        x, y, z = numpy.meshgrid(x, y, z, indexing='ij')
        x = numpy.asfortranarray(x)
        y = numpy.asfortranarray(y)
        z = numpy.asfortranarray(z)

        f = numpy.sin(self.omega*x) * numpy.cos(self.omega*y) * numpy.cos(self.omega*z)
        dfdz_exact = -self.omega * numpy.sin(self.omega*x) * numpy.cos(self.omega*y) * numpy.sin(self.omega*z)

        der = pycd06.cd06stuff.cd06(self.nz, dz, periodic, 0, 0)

        dfdz = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        der.dd3(f, dfdz, self.nx, self.ny)

        error = numpy.absolute(dfdz_exact - dfdz).max()
        self.assertLess(error, 1.0e-9, "Incorrect derivative in first direction for second order finite difference!")
    
    def testDerivativeNonperiodicZ(self):
        """
        Test the non-periodic first derivatives in the Z direction.
        """

        dx, dy, dz = 0., 0., 0.
        periodic = False

        if periodic:
            dx, dy, dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[:-1]
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[:-1]
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[:-1]
        else:
            dx, dy, dz = 2.*numpy.pi / (self.nx-1), 2.*numpy.pi / (self.ny-1), 2.*numpy.pi / (self.nz-1)
            x = numpy.linspace(0., 2.*numpy.pi, num=self.nx)
            y = numpy.linspace(0., 2.*numpy.pi, num=self.ny)
            z = numpy.linspace(0., 2.*numpy.pi, num=self.nz)

        x, y, z = numpy.meshgrid(x, y, z, indexing='ij')
        x = numpy.asfortranarray(x)
        y = numpy.asfortranarray(y)
        z = numpy.asfortranarray(z)

        f = numpy.sin(self.omega*x) * numpy.cos(self.omega*y) * numpy.cos(self.omega*z)
        dfdz_exact = -self.omega * numpy.sin(self.omega*x) * numpy.cos(self.omega*y) * numpy.sin(self.omega*z)

        der = pycd06.cd06stuff.cd06(self.nz, dz, periodic, 0, 0)     

        dfdz = numpy.empty( (self.nx, self.ny, self.nz), dtype=numpy.float64, order='F' )
        der.dd3(f, dfdz, self.nx, self.ny)

        error = numpy.absolute(dfdz_exact - dfdz).max()
        self.assertLess(error, 5.0e-5, "Incorrect derivative in first direction for second order finite difference!")
    
if __name__ == '__main__':
    unittest.main()
