from mpi4py import MPI
import numpy
import unittest

from floatpy.parallel import t3dmod
import floatpy.utilities.reduction as red


class TestReduction(unittest.TestCase):
    
    def setUp(self):
        self.nx, self.ny, self.nz = 64, 32, 16
        self.omega_x, self.omega_y, self.omega_z = 1., 2., 3.

        self.comm = MPI.COMM_WORLD
        self.fcomm = self.comm.py2f()
        self.periodic = numpy.array([True, True, True])

        self.dx, self.dy, self.dz = 2.*numpy.pi / self.nx, 2.*numpy.pi / self.ny, 2.*numpy.pi / self.nz

        self.grid_partition = t3dmod.t3d(self.fcomm, self.nx, self.ny, self.nz, self.periodic )

        self.chunk_3d_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.grid_partition.get_sz3d(self.chunk_3d_size)
        self.grid_partition.get_st3d(self.chunk_3d_lo)
        self.grid_partition.get_en3d(self.chunk_3d_hi)
        self.chunk_3d_lo = self.chunk_3d_lo - 1 # Convert to 0 based indexing
        self.chunk_3d_hi = self.chunk_3d_hi - 1 # Convert to 0 based indexing

        self.x = numpy.linspace(0., 2.*numpy.pi, num=self.nx+1)[self.chunk_3d_lo[0]:self.chunk_3d_hi[0]+1]
        self.y = numpy.linspace(0., 2.*numpy.pi, num=self.ny+1)[self.chunk_3d_lo[1]:self.chunk_3d_hi[1]+1]
        self.z = numpy.linspace(0., 2.*numpy.pi, num=self.nz+1)[self.chunk_3d_lo[2]:self.chunk_3d_hi[2]+1]

        self.x, self.y, self.z = numpy.meshgrid(self.x, self.y, self.z, indexing='ij')
        self.x = numpy.asfortranarray(self.x)
        self.y = numpy.asfortranarray(self.y)
        self.z = numpy.asfortranarray(self.z)

        self.f = numpy.sin(self.omega_x*self.x) * numpy.cos(self.omega_y*self.y) * numpy.cos(self.omega_z*self.z)

        self.avg_x   = red.Reduction(self.grid_partition, ( True, False, False))
        self.avg_y   = red.Reduction(self.grid_partition, (False,  True, False))
        self.avg_z   = red.Reduction(self.grid_partition, (False, False,  True))
        self.avg_xy  = red.Reduction(self.grid_partition, ( True,  True, False))
        self.avg_yz  = red.Reduction(self.grid_partition, (False,  True,  True))
        self.avg_xz  = red.Reduction(self.grid_partition, ( True, False,  True))
        self.avg_xyz = red.Reduction(self.grid_partition, ( True,  True,  True))

    def testAverageX(self):
        """
        Test the averaging in the X direction.
        """

        avg = self.avg_x.average(self.f)

        myerror = numpy.array([ numpy.absolute(avg).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)
        # print "avg_x error = %g" %error[0]

        self.assertLess(error[0], 5.0e-14, "Incorrect average in X!")

        avg[:,:,:] = self.f[8:9,:,:]
        avg_full = self.avg_x.allgather(avg)
        error = avg_full[0:1, self.chunk_3d_lo[1]:self.chunk_3d_hi[1]+1, self.chunk_3d_lo[2]:self.chunk_3d_hi[2]+1] \
                - avg
        error = numpy.absolute(error).max()
        print "allgather X error: ", error
        self.assertLess(error, 5.0e-14, "Incorrect allgather in X!")


    def testAverageY(self):
        """
        Test the averaging in the Y direction.
        """

        avg = self.avg_y.average(self.f)

        myerror = numpy.array([ numpy.absolute(avg).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)
        # print "avg_y error = %g" %error[0]

        self.assertLess(error[0], 5.0e-14, "Incorrect average in Y!")

        avg[:,:,:] = self.f[:,0:1,:]
        avg_full = self.avg_y.allgather(avg)
        error = avg_full[self.chunk_3d_lo[0]:self.chunk_3d_hi[0]+1, 0:1, self.chunk_3d_lo[2]:self.chunk_3d_hi[2]+1] \
                - avg
        error = numpy.absolute(error).max()
        print "allgather Y error: ", error
        self.assertLess(error, 5.0e-14, "Incorrect allgather in Y!")


    def testAverageZ(self):
        """
        Test the averaging in the Z direction.
        """

        avg = self.avg_z.average(self.f)

        myerror = numpy.array([ numpy.absolute(avg).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)
        print "avg_z error = %g" %error[0]

        self.assertLess(error[0], 5.0e-14, "Incorrect average in Z!")

        avg[:,:,:] = self.f[:,:,0:1]
        avg_full = self.avg_z.allgather(avg)
        error = avg_full[self.chunk_3d_lo[0]:self.chunk_3d_hi[0]+1, self.chunk_3d_lo[1]:self.chunk_3d_hi[1]+1, 0:1] \
                - avg
        error = numpy.absolute(error).max()
        print "allgather Z error: ", error
        self.assertLess(error, 5.0e-14, "Incorrect allgather in Z!")


    def testAverageXY(self):
        """
        Test the averaging in the X-Y directions.
        """

        avg = self.avg_xy.average(self.f)

        myerror = numpy.array([ numpy.absolute(avg).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)
        # print "avg_xy error = %g" %error[0]

        self.assertLess(error[0], 5.0e-14, "Incorrect average in X-Y!")


    def testAverageYZ(self):
        """
        Test the averaging in the Y-Z directions.
        """

        avg = self.avg_yz.average(self.f)

        myerror = numpy.array([ numpy.absolute(avg).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)
        # print "avg_yz error = %g" %error[0]

        self.assertLess(error[0], 5.0e-14, "Incorrect average in Y-Z!")


    def testAverageXZ(self):
        """
        Test the averaging in the X-Z directions.
        """

        avg = self.avg_xz.average(self.f)

        myerror = numpy.array([ numpy.absolute(avg).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)
        # print "avg_xz error = %g" %error[0]

        self.assertLess(error[0], 5.0e-14, "Incorrect average in X-Z!")


    def testAverageXYZ(self):
        """
        Test the averaging in the X-Y-Z directions.
        """

        avg = self.avg_xyz.average(self.f)

        myerror = numpy.array([ numpy.absolute(avg).max() ])
        error = numpy.array([ myerror[0] ])
        self.comm.Allreduce(myerror, error, op=MPI.MAX)
        # print "avg_xyz error = %g" %error[0]

        self.assertLess(error[0], 5.0e-14, "Incorrect average in X-Y-Z!")


if __name__ == '__main__':
    unittest.main()
