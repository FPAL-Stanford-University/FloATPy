import numpy
import unittest

import floatpy.upsampling.Lagrange_upsampler

class TestLagrangeUnsampler(unittest.TestCase):
    
    def testUpsamling1D(self):
        """
        Test the 1D upsampling.
        """
        
        r = numpy.array([10])
        
        x_constant = numpy.linspace(0, 1.0, 2)
        y_constant = numpy.ones([2])
        
        upsampler_constant = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='constant')
        y_upsampled_constant = upsampler_constant.upsample(y_constant, r)
        
        error_constant = numpy.absolute(y_upsampled_constant - 1.0).max()
        self.assertLess(error_constant, 1.0e-10, "Incorrect 1D constant upsampling!")
        
        
        x_second = numpy.linspace(0.0, 2.0, 3)
        y_second = x_second
        
        upsampler_second = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='second_order')
        y_upsampled_second = upsampler_second.upsample(y_second, r)
        
        n_ghosts_second = upsampler_second.getNumberOfGhostCells()
        x_second_upsampled = numpy.linspace(0.0 + 0.5/r[0], 2.0 - 0.5/r[0], 2*r[0])
        y_upsampled_second_correct = x_second_upsampled
        error_second = numpy.absolute(
            y_upsampled_second[(n_ghosts_second*r[0]-r[0]/2):(3*r[0]-(n_ghosts_second*r[0]-r[0]/2))] - \
            y_upsampled_second_correct).max()
        self.assertLess(error_second, 1.0e-10, "Incorrect 1D second order Lagrange upsampling!")
        
        x_fourth = numpy.linspace(0.0, 4.0, 5)
        y_fourth = x_fourth*x_fourth*x_fourth
        
        upsampler_fourth = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='fourth_order')
        y_upsampled_fourth = upsampler_fourth.upsample(y_fourth, r)
        
        n_ghosts_fourth = upsampler_fourth.getNumberOfGhostCells()
        x_fourth_upsampled = numpy.linspace(1.0 + 0.5/r[0], 3.0 - 0.5/r[0], 2*r[0])
        y_upsampled_fourth_correct = x_fourth_upsampled*x_fourth_upsampled*x_fourth_upsampled
        error_fourth = numpy.absolute(
            y_upsampled_fourth[(n_ghosts_fourth*r[0]-r[0]/2):(5*r[0]-(n_ghosts_fourth*r[0]-r[0]/2))] - \
            y_upsampled_fourth_correct).max()
        self.assertLess(error_fourth, 1.0e-10, "Incorrect 1D fourth order Lagrange upsampling!")
        
        x_sixth = numpy.linspace(0.0, 6.0, 7)
        y_sixth = x_sixth*x_sixth*x_sixth*x_sixth*x_sixth
        
        upsampler_sixth = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='sixth_order')
        y_upsampled_sixth = upsampler_sixth.upsample(y_sixth, r)
        
        n_ghosts_sixth = upsampler_sixth.getNumberOfGhostCells()
        x_sixth_upsampled = numpy.linspace(2.0 + 0.5/r[0], 4.0 - 0.5/r[0], 2*r[0])
        y_upsampled_sixth_correct = x_sixth_upsampled*x_sixth_upsampled*x_sixth_upsampled*x_sixth_upsampled*x_sixth_upsampled
        error_sixth = numpy.absolute(
            y_upsampled_sixth[(n_ghosts_sixth*r[0]-r[0]/2):(7*r[0]-(n_ghosts_sixth*r[0]-r[0]/2))] - \
            y_upsampled_sixth_correct).max()
        self.assertLess(error_sixth, 1.0e-10, "Incorrect 1D sixth order Lagrange upsampling!")
    
    
    def testUpsamling2D(self):
        """
        Test the 2D upsampling.
        """
        
        r = numpy.array([1, 10])
        
        x_constant = numpy.empty([1, 2])
        x_constant[0, :] = numpy.linspace(0, 1.0, 2)
        x_constant = numpy.repeat(x_constant, 2, axis=0)
        y_constant = numpy.ones(x_constant.shape)
        
        upsampler_constant = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='constant')
        y_upsampled_constant = upsampler_constant.upsample(y_constant, r)
        
        error_constant = numpy.absolute(y_upsampled_constant - 1.0).max()
        self.assertLess(error_constant, 1.0e-10, "Incorrect 2D constant upsampling!")
        
        x_second = numpy.empty([1, 3])
        x_second[0, :] = numpy.linspace(0.0, 2.0, 3)
        x_second = numpy.repeat(x_second, 3, axis=0)
        y_second = x_second
        
        upsampler_second = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='second_order')
        y_upsampled_second = upsampler_second.upsample(y_second, r)
        
        n_ghosts_second = upsampler_second.getNumberOfGhostCells()
        x_second_upsampled = numpy.empty([1, 2*r[1]])
        x_second_upsampled[0, :] = numpy.linspace(0.0 + 0.5/r[1], 2.0 - 0.5/r[1], 2*r[1])
        x_second_upsampled = numpy.repeat(x_second_upsampled, 3, axis=0)
        y_upsampled_second_correct = x_second_upsampled
        error_second = numpy.absolute(\
            y_upsampled_second[:, (n_ghosts_second*r[1]-r[1]/2):(3*r[1]-(n_ghosts_second*r[1]-r[1]/2))] - \
            y_upsampled_second_correct).max()
        self.assertLess(error_second, 1.0e-10, "Incorrect 2D second order Lagrange upsampling!")
        
        x_fourth = numpy.empty([1, 5])
        x_fourth[0, :] = numpy.linspace(0.0, 4.0, 5)
        x_fourth = numpy.repeat(x_fourth, 5, axis=0)
        y_fourth = x_fourth*x_fourth*x_fourth
        
        upsampler_fourth = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='fourth_order')
        y_upsampled_fourth = upsampler_fourth.upsample(y_fourth, r)
        
        n_ghosts_fourth = upsampler_fourth.getNumberOfGhostCells()
        x_fourth_upsampled = numpy.empty([1, 2*r[1]])
        x_fourth_upsampled[0, :]  = numpy.linspace(1.0 + 0.5/r[1], 3.0 - 0.5/r[1], 2*r[1])
        x_fourth_upsampled = numpy.repeat(x_fourth_upsampled, 5, axis=0)
        y_upsampled_fourth_correct = x_fourth_upsampled*x_fourth_upsampled*x_fourth_upsampled
        error_fourth = numpy.absolute(
            y_upsampled_fourth[1:-1, (n_ghosts_fourth*r[1]-r[1]/2):(5*r[1]-(n_ghosts_fourth*r[1]-r[1]/2))] - \
            y_upsampled_fourth_correct[1:-1, :]).max()
        self.assertLess(error_fourth, 1.0e-10, "Incorrect 2D fourth order Lagrange upsampling!")
        
        x_sixth = numpy.empty([1, 7])
        x_sixth[0, :] = numpy.linspace(0.0, 6.0, 7)
        x_sixth = numpy.repeat(x_sixth, 7, axis=0)
        y_sixth = x_sixth*x_sixth*x_sixth*x_sixth*x_sixth
        
        upsampler_sixth = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='sixth_order')
        y_upsampled_sixth = upsampler_sixth.upsample(y_sixth, r)
        
        n_ghosts_sixth = upsampler_sixth.getNumberOfGhostCells()
        x_sixth_upsampled = numpy.empty([1, 2*r[1]])
        x_sixth_upsampled[0, :] = numpy.linspace(2.0 + 0.5/r[1], 4.0 - 0.5/r[1], 2*r[1])
        x_sixth_upsampled = numpy.repeat(x_sixth_upsampled, 7, axis=0)
        y_upsampled_sixth_correct = x_sixth_upsampled*x_sixth_upsampled*x_sixth_upsampled*x_sixth_upsampled*x_sixth_upsampled
        error_sixth = numpy.absolute(
            y_upsampled_sixth[2:-2, (n_ghosts_sixth*r[1]-r[1]/2):(7*r[1]-(n_ghosts_sixth*r[1]-r[1]/2))] - \
            y_upsampled_sixth_correct[2:-2, :]).max()
        self.assertLess(error_sixth, 1.0e-10, "Incorrect 2D sixth order Lagrange upsampling!")
    
    
    def testUpsamling3D(self):
        """
        Test the 3D upsampling.
        """
        
        r = numpy.array([1, 1, 10])
        
        x_constant = numpy.empty([1, 1, 2])
        x_constant[0, 0,:] = numpy.linspace(0, 1.0, 2)
        x_constant = numpy.repeat(x_constant, 2, axis=0)
        x_constant = numpy.repeat(x_constant, 2, axis=1)
        y_constant = numpy.ones(x_constant.shape)
        
        upsampler_constant = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='constant')
        y_upsampled_constant = upsampler_constant.upsample(y_constant, r)
        
        error_constant = numpy.absolute(y_upsampled_constant - 1.0).max()
        self.assertLess(error_constant, 1.0e-10, "Incorrect 3D constant upsampling!")
        
        x_second = numpy.empty([1, 1, 3])
        x_second[0, 0, :] = numpy.linspace(0.0, 2.0, 3)
        x_second = numpy.repeat(x_second, 3, axis=0)
        x_second = numpy.repeat(x_second, 3, axis=1)
        y_second = x_second
        
        upsampler_second = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='second_order')
        y_upsampled_second = upsampler_second.upsample(y_second, r)
        
        n_ghosts_second = upsampler_second.getNumberOfGhostCells()
        x_second_upsampled = numpy.empty([1, 1, 2*r[2]])
        x_second_upsampled[0, 0, :] = numpy.linspace(0.0 + 0.5/r[2], 2.0 - 0.5/r[2], 2*r[2])
        x_second_upsampled = numpy.repeat(x_second_upsampled, 3, axis=0)
        x_second_upsampled = numpy.repeat(x_second_upsampled, 3, axis=1)
        y_upsampled_second_correct = x_second_upsampled
        error_second = numpy.absolute(
            y_upsampled_second[:, :, (n_ghosts_second*r[2]-r[2]/2):(3*r[2]-(n_ghosts_second*r[2]-r[2]/2))] - \
            y_upsampled_second_correct).max()
        self.assertLess(error_second, 1.0e-10, "Incorrect 3D second order Lagrange upsampling!")
        
        x_fourth = numpy.empty([1, 1, 5])
        x_fourth[0, 0, :] = numpy.linspace(0.0, 4.0, 5)
        x_fourth = numpy.repeat(x_fourth, 5, axis=0)
        x_fourth = numpy.repeat(x_fourth, 5, axis=1)
        y_fourth = x_fourth*x_fourth*x_fourth
        
        upsampler_fourth = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='fourth_order')
        y_upsampled_fourth = upsampler_fourth.upsample(y_fourth, r)
        
        n_ghosts_fourth = upsampler_fourth.getNumberOfGhostCells()
        x_fourth_upsampled = numpy.empty([1, 1, 2*r[2]])
        x_fourth_upsampled[0, 0, :]  = numpy.linspace(1.0 + 0.5/r[2], 3.0 - 0.5/r[2], 2*r[2])
        x_fourth_upsampled = numpy.repeat(x_fourth_upsampled, 5, axis=0)
        x_fourth_upsampled = numpy.repeat(x_fourth_upsampled, 5, axis=1)
        y_upsampled_fourth_correct = x_fourth_upsampled*x_fourth_upsampled*x_fourth_upsampled
        error_fourth = numpy.absolute(
            y_upsampled_fourth[1:-1, 1:-1, (n_ghosts_fourth*r[2]-r[2]/2):(5*r[2]-(n_ghosts_fourth*r[2]-r[2]/2))] - \
            y_upsampled_fourth_correct[1:-1, 1:-1, :]).max()
        self.assertLess(error_fourth, 1.0e-10, "Incorrect 3D fourth order Lagrange upsampling!")
        
        x_sixth = numpy.empty([1, 1, 7])
        x_sixth[0, 0, :] = numpy.linspace(0.0, 6.0, 7)
        x_sixth = numpy.repeat(x_sixth, 7, axis=0)
        x_sixth = numpy.repeat(x_sixth, 7, axis=1)
        y_sixth = x_sixth*x_sixth*x_sixth*x_sixth*x_sixth
        
        upsampler_sixth = floatpy.upsampling.Lagrange_upsampler.LagrangeUpsampler(method='sixth_order')
        y_upsampled_sixth = upsampler_sixth.upsample(y_sixth, r)
        
        n_ghosts_sixth = upsampler_sixth.getNumberOfGhostCells()
        x_sixth_upsampled = numpy.empty([1, 1, 2*r[2]])
        x_sixth_upsampled[0, 0, :] = numpy.linspace(2.0 + 0.5/r[2], 4.0 - 0.5/r[2], 2*r[2])
        x_sixth_upsampled = numpy.repeat(x_sixth_upsampled, 7, axis=0)
        x_sixth_upsampled = numpy.repeat(x_sixth_upsampled, 7, axis=1)
        y_upsampled_sixth_correct = x_sixth_upsampled*x_sixth_upsampled*x_sixth_upsampled*x_sixth_upsampled*x_sixth_upsampled
        error_sixth = numpy.absolute(
            y_upsampled_sixth[2:-2, 2:-2, (n_ghosts_sixth*r[2]-r[2]/2):(7*r[2]-(n_ghosts_sixth*r[2]-r[2]/2))] - \
            y_upsampled_sixth_correct[2:-2, 2:-2, :]).max()
        self.assertLess(error_sixth, 1.0e-10, "Incorrect 3D sixth order Lagrange upsampling!")


if __name__ == '__main__':
    unittest.main()
