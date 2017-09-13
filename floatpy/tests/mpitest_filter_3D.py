from mpi4py import MPI
import numpy
import unittest

from floatpy.parallel import _t3dmod
from floatpy.filters import Filter

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


class TestFilter(unittest.TestCase):
    
    def setUp(self):
        self.nx, self.ny, self.nz = 32, 32, 32
        self.omega = 12.

        self.comm = MPI.COMM_WORLD
        self.fcomm = self.comm.py2f()
        self.periodic = numpy.array([True, True, True])
        self.filter_type = ('compact', 'compact', 'compact')

        self.grid_partition = _t3dmod.t3d(self.fcomm, self.nx, self.ny, self.nz, self.periodic )

        self.chunk_3d_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.grid_partition.get_sz3d(self.chunk_3d_size)
        self.grid_partition.get_st3d(self.chunk_3d_lo)
        self.grid_partition.get_en3d(self.chunk_3d_hi)
        self.chunk_3d_lo = self.chunk_3d_lo - 1 # Convert to 0 based indexing
        self.chunk_3d_hi = self.chunk_3d_hi - 1 # Convert to 0 based indexing

        self.fil = Filter( self.grid_partition, self.filter_type, periodic_dimensions=self.periodic )

        self.x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[self.chunk_3d_lo[0]:self.chunk_3d_hi[0]+1]
        self.y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[self.chunk_3d_lo[1]:self.chunk_3d_hi[1]+1]
        self.z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[self.chunk_3d_lo[2]:self.chunk_3d_hi[2]+1]

        self.x, self.y, self.z = numpy.meshgrid(self.x, self.y, self.z, indexing='ij')
        self.x = numpy.asfortranarray(self.x)
        self.y = numpy.asfortranarray(self.y)
        self.z = numpy.asfortranarray(self.z)

        self.f = numpy.sin(self.omega*self.x) * numpy.cos(self.omega*self.y) * numpy.cos(self.omega*self.z)

        self.dx, self.dy, self.dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz
        k_norm_x = self.omega * self.dx
        TF_x = getTransferFunction(k_norm_x)
        k_norm_y = self.omega * self.dy
        TF_y = getTransferFunction(k_norm_y)
        k_norm_z = self.omega * self.dz
        TF_z = getTransferFunction(k_norm_z)

        self.f_tilde_x_exact = TF_x * self.f
        self.f_tilde_y_exact = TF_y * self.f
        self.f_tilde_z_exact = TF_z * self.f

        self.f_tilde_exact = TF_x * TF_y * TF_z * self.f


    def testFilterPeriodicX(self):
        """
        Test the filter in the X direction.
        """

        f_tilde = numpy.empty( self.chunk_3d_size, dtype=numpy.float64, order='F' )
        self.fil.filter_x(self.f, f_tilde)

        myerror = numpy.zeros(1)
        myerror[0] = numpy.absolute(self.f_tilde_x_exact - f_tilde).max()

        error = numpy.zeros(1)
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 5.0e-14, "Incorrect filter in X direction!")


    def testFilterPeriodicY(self):
        """
        Test the filter in the Y direction.
        """

        f_tilde = numpy.empty( self.chunk_3d_size, dtype=numpy.float64, order='F' )
        self.fil.filter_y(self.f, f_tilde)

        myerror = numpy.zeros(1)
        myerror[0] = numpy.absolute(self.f_tilde_y_exact - f_tilde).max()

        error = numpy.zeros(1)
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 5.0e-14, "Incorrect filter in Y direction!")
    
    
    def testFilterPeriodicZ(self):
        """
        Test the filter in the Z direction.
        """

        f_tilde = numpy.empty( self.chunk_3d_size, dtype=numpy.float64, order='F' )
        self.fil.filter_z(self.f, f_tilde)

        myerror = numpy.zeros(1)
        myerror[0] = numpy.absolute(self.f_tilde_z_exact - f_tilde).max()

        error = numpy.zeros(1)
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 5.0e-14, "Incorrect compact filter in Z direction!")
    
    
    def testFilterPeriodic3D(self):
        """
        Test the periodic 3D filter.
        """

        f_tilde = self.fil.filter_all(self.f)

        myerror = numpy.zeros(1)
        myerror[0] = numpy.absolute(self.f_tilde_exact - f_tilde).max()

        error = numpy.zeros(1)
        self.comm.Allreduce(myerror, error, op=MPI.MAX)

        self.assertLess(error[0], 1.0e-13, "Incorrect 3D filter!")
    
    
if __name__ == '__main__':
    unittest.main()
