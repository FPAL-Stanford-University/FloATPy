import numpy
import unittest

import floatpy.derivatives.first_order_derivative

class TestDerivativesFirst(unittest.TestCase):
    
    def testFirstDirection(self):
        """
        Test the first derivatives in the first direction.
        """
        
        x = numpy.empty([1, 100])
        x[0, :] = numpy.linspace(0, 2*numpy.pi, 100)
        dx = x[0, 1] - x[0, 0]
        
        y = numpy.sin(x)
        y_prime_exact = numpy.cos(x)
        
        dydx_2 = floatpy.derivatives.first_order_derivative.FirstOrderDerivative('second_order', direction=0, data_order='C')
        dydx_4 = floatpy.derivatives.first_order_derivative.FirstOrderDerivative('fourth_order', direction=0, data_order='C')
        dydx_6 = floatpy.derivatives.first_order_derivative.FirstOrderDerivative('sixth_order', direction=0, data_order='C')
        
        y_prime_2 = dydx_2.differentiate(y, dx, 0, True)
        y_prime_4 = dydx_4.differentiate(y, dx, 0, True)
        y_prime_6 = dydx_6.differentiate(y, dx, 0, True)
        
        error_2 = numpy.linalg.norm(y_prime_exact[0, :] - y_prime_2[:])
        error_4 = numpy.linalg.norm(y_prime_exact[0, :] - y_prime_4[:])
        error_6 = numpy.linalg.norm(y_prime_exact[0, :] - y_prime_6[:])
        
        self.assertLess(error_2, 1.0e-2, "Incorrect derivative in first direction for second order finite difference!")
        self.assertLess(error_4, 1.0e-5, "Incorrect derivative in first direction for fourth order finite difference!")
        self.assertLess(error_6, 1.0e-7, "Incorrect derivative in first direction for sixth order finite difference!")
    
    
    def testSecondDirection(self):
        """
        Test the first derivatives in the second direction.
        """
        
        x = numpy.empty([1, 1, 100])
        x[0, 0, :] = numpy.linspace(0, 2*numpy.pi, 100)
        dx = x[0, 0, 1] - x[0, 0, 0]
        
        y = numpy.sin(x)
        y_prime_exact = numpy.cos(x)
        
        dydx_2 = floatpy.derivatives.first_order_derivative.FirstOrderDerivative('second_order', direction=1, data_order='C')
        dydx_4 = floatpy.derivatives.first_order_derivative.FirstOrderDerivative('fourth_order', direction=1, data_order='C')
        dydx_6 = floatpy.derivatives.first_order_derivative.FirstOrderDerivative('sixth_order', direction=1, data_order='C')
        
        y_prime_2 = dydx_2.differentiate(y, dx, 0, True)
        y_prime_4 = dydx_4.differentiate(y, dx, 0, True)
        y_prime_6 = dydx_6.differentiate(y, dx, 0, True)
        
        error_2 = numpy.linalg.norm(y_prime_exact[0, 0, :] - y_prime_2[0, :])
        error_4 = numpy.linalg.norm(y_prime_exact[0, 0, :] - y_prime_4[0, :])
        error_6 = numpy.linalg.norm(y_prime_exact[0, 0, :] - y_prime_6[0, :])
        
        self.assertLess(error_2, 1.0e-2, "Incorrect derivative in second direction for second order finite difference!")
        self.assertLess(error_4, 1.0e-5, "Incorrect derivative in second direction for fourth order finite difference!")
        self.assertLess(error_6, 1.0e-7, "Incorrect derivative in second direction for sixth order finite difference!")
    
    
    def testThirdDirection(self):
        """
        Test the first derivatives in the thrid direction.
        """
        
        x = numpy.empty([1, 1, 1, 100])
        x[0, 0, 0, :] = numpy.linspace(0, 2*numpy.pi, 100)
        dx = x[0, 0, 0, 1] - x[0, 0, 0, 0]
        
        y = numpy.sin(x)
        y_prime_exact = numpy.cos(x)
        
        dydx_2 = floatpy.derivatives.first_order_derivative.FirstOrderDerivative('second_order', direction=2, data_order='C')
        dydx_4 = floatpy.derivatives.first_order_derivative.FirstOrderDerivative('fourth_order', direction=2, data_order='C')
        dydx_6 = floatpy.derivatives.first_order_derivative.FirstOrderDerivative('sixth_order', direction=2, data_order='C')
        
        y_prime_2 = dydx_2.differentiate(y, dx, 0, True)
        y_prime_4 = dydx_4.differentiate(y, dx, 0, True)
        y_prime_6 = dydx_6.differentiate(y, dx, 0, True)
        
        error_2 = numpy.linalg.norm(y_prime_exact[0, 0, 0, :] - y_prime_2[0, 0, :])
        error_4 = numpy.linalg.norm(y_prime_exact[0, 0, 0, :] - y_prime_4[0, 0, :])
        error_6 = numpy.linalg.norm(y_prime_exact[0, 0, 0, :] - y_prime_6[0, 0, :])
        
        self.assertLess(error_2, 1.0e-2, "Incorrect derivative in third direction for second order finite difference")
        self.assertLess(error_4, 1.0e-5, "Incorrect derivative in third direction for fourth order finite difference")
        self.assertLess(error_6, 1.0e-7, "Incorrect derivative in third direction for sixth order finite difference")
    
    
if __name__ == '__main__':
    unittest.main()
