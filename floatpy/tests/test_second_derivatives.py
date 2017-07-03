import numpy
import unittest

import floatpy.derivatives.second

class TestSecondDerivatives(unittest.TestCase):
    
    def TestFirstDirection(self):
        """
        Test the second derivatives in the first direction.
        """
        
        x = numpy.empty([1, 100])
        x[0, :] = numpy.linspace(0, 2*numpy.pi, 100)
        dx = x[0, 1] - x[0, 0]
        
        y = numpy.sin(x)
        y_prime_exact = -numpy.sin(x)
        
        y_prime_2 = floatpy.derivatives.second.computeSecondDerivative(y, dx, 0, 0, method = 'second_order')
        y_prime_4 = floatpy.derivatives.second.computeSecondDerivative(y, dx, 0, 0, method = 'fourth_order')
        y_prime_6 = floatpy.derivatives.second.computeSecondDerivative(y, dx, 0, 0, method = 'sixth_order')
        
        error_2 = numpy.linalg.norm(y_prime_exact[0, :] - y_prime_2[:])
        error_4 = numpy.linalg.norm(y_prime_exact[0, :] - y_prime_4[:])
        error_6 = numpy.linalg.norm(y_prime_exact[0, :] - y_prime_6[:])
        
        self.assertLess(error_2, 1.0e-2, "Incorrect derivative in first direction for second order finite difference")
        self.assertLess(error_4, 1.0e-5, "Incorrect derivative in first direction for fourth order finite difference")
        self.assertLess(error_6, 1.0e-7, "Incorrect derivative in first direction for sixth order finite difference")
    
    
    def TestSecondDirection(self):
        """
        Test the second derivatives in the second direction.
        """
        
        x = numpy.empty([1, 1, 100])
        x[0, 0, :] = numpy.linspace(0, 2*numpy.pi, 100)
        dx = x[0, 0, 1] - x[0, 0, 0]
        
        y = numpy.sin(x)
        y_prime_exact = -numpy.sin(x)
        
        y_prime_2 = floatpy.derivatives.second.computeSecondDerivative(y, dx, 1, 0, method = 'second_order')
        y_prime_4 = floatpy.derivatives.second.computeSecondDerivative(y, dx, 1, 0, method = 'fourth_order')
        y_prime_6 = floatpy.derivatives.second.computeSecondDerivative(y, dx, 1, 0, method = 'sixth_order')
        
        error_2 = numpy.linalg.norm(y_prime_exact[0, 0, :] - y_prime_2[0, :])
        error_4 = numpy.linalg.norm(y_prime_exact[0, 0, :] - y_prime_4[0, :])
        error_6 = numpy.linalg.norm(y_prime_exact[0, 0, :] - y_prime_6[0, :])
        
        self.assertLess(error_2, 1.0e-2, "Incorrect derivative in second direction for second order finite difference")
        self.assertLess(error_4, 1.0e-5, "Incorrect derivative in second direction for fourth order finite difference")
        self.assertLess(error_6, 1.0e-7, "Incorrect derivative in second direction for sixth order finite difference")
    
    
    def TestThirdDirection(self):
        """
        Test the second derivatives in the thrid direction.
        """
        
        x = numpy.empty([1, 1, 1, 100])
        x[0, 0, 0, :] = numpy.linspace(0, 2*numpy.pi, 100)
        dx = x[0, 0, 0, 1] - x[0, 0, 0, 0]
        
        y = numpy.sin(x)
        y_prime_exact = -numpy.sin(x)
        
        y_prime_2 = floatpy.derivatives.second.computeSecondDerivative(y, dx, 2, 0, method = 'second_order')
        y_prime_4 = floatpy.derivatives.second.computeSecondDerivative(y, dx, 2, 0, method = 'fourth_order')
        y_prime_6 = floatpy.derivatives.second.computeSecondDerivative(y, dx, 2, 0, method = 'sixth_order')
        
        error_2 = numpy.linalg.norm(y_prime_exact[0, 0, 0, :] - y_prime_2[0, 0, :])
        error_4 = numpy.linalg.norm(y_prime_exact[0, 0, 0, :] - y_prime_4[0, 0, :])
        error_6 = numpy.linalg.norm(y_prime_exact[0, 0, 0, :] - y_prime_6[0, 0, :])
        
        self.assertLess(error_2, 1.0e-2, "Incorrect derivative in third direction for second order finite difference")
        self.assertLess(error_4, 1.0e-5, "Incorrect derivative in third direction for fourth order finite difference")
        self.assertLess(error_6, 1.0e-7, "Incorrect derivative in third direction for sixth order finite difference")
    
    
if __name__ == '__main__':
    unittest.main()
