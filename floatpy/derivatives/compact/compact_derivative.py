import numpy
import pycd06
import pycd10

from floatpy.parallel import _t3dmod

class CompactDerivative(object):
    """
    Class to perform parallel compact finite difference routines
    """

    def __init__(self, grid_partition, grid_spacing, order, periodic=(False, False, False)):
        """
        Constructor of the class

        grid_partition : t3d object or the grid_partition property of the parallel data reader class
        grid_spacing : float tuple of size 3 with the grid spacing in each direction
        order : integer tuple of size 3 with each value in {6,10} representing the order of accuracy 
                of the derivatives in each direction
        periodic : boolean tuple of size 3
        """

        if not isinstance(grid_partition, _t3dmod.t3d):
            raise RuntimeError("The given grid partition object is not an instance of the t3d class!")

        for i in range(3):
            if order[i] not in [6,10]:
                raise RuntimeError("order[%d] has to be one of {6,10}" %i)

        self._order = order
        self._grid_partition = grid_partition

        self._periodic = tuple(periodic)

        self._chunk_3d_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_3d_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_3d_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._grid_partition.get_sz3d(self._chunk_3d_size)
        self._grid_partition.get_st3d(self._chunk_3d_lo)
        self._grid_partition.get_en3d(self._chunk_3d_hi)
        self._chunk_3d_lo = self._chunk_3d_lo - 1 # Convert to 0 based indexing
        self._chunk_3d_hi = self._chunk_3d_hi - 1 # Convert to 0 based indexing

        self._chunk_x_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_x_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_x_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._grid_partition.get_szx(self._chunk_x_size)
        self._grid_partition.get_stx(self._chunk_x_lo)
        self._grid_partition.get_enx(self._chunk_x_hi)
        self._chunk_x_lo = self._chunk_x_lo - 1 # Convert to 0 based indexing
        self._chunk_x_hi = self._chunk_x_hi - 1 # Convert to 0 based indexing

        self._chunk_y_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_y_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_y_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._grid_partition.get_szy(self._chunk_y_size)
        self._grid_partition.get_sty(self._chunk_y_lo)
        self._grid_partition.get_eny(self._chunk_y_hi)
        self._chunk_y_lo = self._chunk_y_lo - 1 # Convert to 0 based indexing
        self._chunk_y_hi = self._chunk_y_hi - 1 # Convert to 0 based indexing

        self._chunk_z_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_z_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_z_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._grid_partition.get_szz(self._chunk_z_size)
        self._grid_partition.get_stz(self._chunk_z_lo)
        self._grid_partition.get_enz(self._chunk_z_hi)
        self._chunk_z_lo = self._chunk_z_lo - 1 # Convert to 0 based indexing
        self._chunk_z_hi = self._chunk_z_hi - 1 # Convert to 0 based indexing

        self._nx = self._chunk_x_size[0]
        self._ny = self._chunk_y_size[1]
        self._nz = self._chunk_z_size[2]

        self._dx = grid_spacing[0]
        self._dy = grid_spacing[1]
        self._dz = grid_spacing[2]

        if self._order[0] == 6:
            self._xder = pycd06.cd06stuff.cd06( self._nx, self._dx, self._periodic[0], 0, 0 )
        elif self._order[0] == 10:
            self._xder = pycd10.cd10stuff.cd10( self._nx, self._dx, self._periodic[0], 0, 0 )

        if self._order[1] == 6:
            self._yder = pycd06.cd06stuff.cd06( self._ny, self._dy, self._periodic[1], 0, 0 )
        elif self._order[1] == 10:
            self._yder = pycd10.cd10stuff.cd10( self._ny, self._dy, self._periodic[1], 0, 0 )

        if self._order[2] == 6:
            self._zder = pycd06.cd06stuff.cd06( self._nz, self._dz, self._periodic[2], 0, 0 )
        elif self._order[2] == 10:
            self._zder = pycd10.cd10stuff.cd10( self._nz, self._dz, self._periodic[2], 0, 0 )


    def ddx(self, f, dfdx, bc=(0,0)):
        """
        Function to compute the X derivative of f

        f : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
            This array must be in the 3D decomposition
        dfdx : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
               This array is in the 3D decomposition
        bc : integer tuple of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        if (f.shape[0] != self._chunk_3d_size[0]) or (f.shape[1] != self._chunk_3d_size[1]) or (f.shape[2] != self._chunk_3d_size[2]):
            raise RuntimeError("Make sure f is of the same size as in grid_partition!")

        f_x    = numpy.empty( self._chunk_x_size, dtype=numpy.float64, order='F' )
        dfdx_x = numpy.empty( self._chunk_x_size, dtype=numpy.float64, order='F' )

        self._grid_partition.transpose_3d_to_x( f, f_x )
        if self._order[0] == 6:
            self._xder.dd1(f_x, dfdx_x, self._chunk_x_size[1], self._chunk_x_size[2]) # symmetry BC only supported in 10th order for now
        elif self._order[0] == 10:
            self._xder.dd1(f_x, dfdx_x, self._chunk_x_size[1], self._chunk_x_size[2], bc1_=bc[0], bcn_=bc[1])
        self._grid_partition.transpose_x_to_3d( dfdx_x, dfdx )


    def ddy(self, f, dfdy, bc=(0,0)):
        """
        Function to compute the Y derivative of f

        f : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
            This array must be in the 3D decomposition
        dfdy : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
               This array is in the 3D decomposition
        bc : integer tuple of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        if (f.shape[0] != self._chunk_3d_size[0]) or (f.shape[1] != self._chunk_3d_size[1]) or (f.shape[2] != self._chunk_3d_size[2]):
            raise RuntimeError("Make sure f is of the same size as in grid_partition!")

        f_y    = numpy.empty( self._chunk_y_size, dtype=numpy.float64, order='F' )
        dfdy_y = numpy.empty( self._chunk_y_size, dtype=numpy.float64, order='F' )

        self._grid_partition.transpose_3d_to_y( f, f_y )
        if self._order[1] == 6:
            self._yder.dd2(f_y, dfdy_y, self._chunk_y_size[0], self._chunk_y_size[2]) # symmetry BC only supported in 10th order for now
        elif self._order[1] == 10:
            self._yder.dd2(f_y, dfdy_y, self._chunk_y_size[0], self._chunk_y_size[2], bc1_=bc[0], bcn_=bc[1])
        self._grid_partition.transpose_y_to_3d( dfdy_y, dfdy )


    def ddz(self, f, dfdz, bc=(0,0)):
        """
        Function to compute the Z derivative of f

        f : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
            This array must be in the 3D decomposition
        dfdy : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
               This array is in the 3D decomposition
        bc : integer tuple of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        if (f.shape[0] != self._chunk_3d_size[0]) or (f.shape[1] != self._chunk_3d_size[1]) or (f.shape[2] != self._chunk_3d_size[2]):
            raise RuntimeError("Make sure f is of the same size as in grid_partition!")

        f_z    = numpy.empty( self._chunk_z_size, dtype=numpy.float64, order='F' )
        dfdz_z = numpy.empty( self._chunk_z_size, dtype=numpy.float64, order='F' )

        self._grid_partition.transpose_3d_to_z( f, f_z )
        if self._order[2] == 6:
            self._zder.dd3(f_z, dfdz_z, self._chunk_z_size[0], self._chunk_z_size[1]) # symmetry BC only supported in 10th order for now
        elif self._order[2] == 10:
            self._zder.dd3(f_z, dfdz_z, self._chunk_z_size[0], self._chunk_z_size[1], bc1_=bc[0], bcn_=bc[1])
        self._grid_partition.transpose_z_to_3d( dfdz_z, dfdz )


    def d2dx2(self, f, d2fdx2, bc=(0,0)):
        """
        Function to compute the X second derivative of f

        f : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
            This array must be in the 3D decomposition
        d2fdx2 : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
               This array is in the 3D decomposition
        bc : integer tuple of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        if (f.shape[0] != self._chunk_3d_size[0]) or (f.shape[1] != self._chunk_3d_size[1]) or (f.shape[2] != self._chunk_3d_size[2]):
            raise RuntimeError("Make sure f is of the same size as in grid_partition!")

        f_x      = numpy.empty( self._chunk_x_size, dtype=numpy.float64, order='F' )
        d2fdx2_x = numpy.empty( self._chunk_x_size, dtype=numpy.float64, order='F' )

        self._grid_partition.transpose_3d_to_x( f, f_x )
        if self._order[0] == 6:
            raise NotImplementedError("6th order 2nd derivatives are not implemented yet. Sorry!")
        elif self._order[0] == 10:
            self._xder.d2d1(f_x, d2fdx2_x, self._chunk_x_size[1], self._chunk_x_size[2], bc1_=bc[0], bcn_=bc[1])
        self._grid_partition.transpose_x_to_3d( d2fdx2_x, d2fdx2 )


    def d2dy2(self, f, d2fdy2, bc=(0,0)):
        """
        Function to compute the Y second derivative of f

        f : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
            This array must be in the 3D decomposition
        d2fdy2 : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
               This array is in the 3D decomposition
        bc : integer tuple of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        if (f.shape[0] != self._chunk_3d_size[0]) or (f.shape[1] != self._chunk_3d_size[1]) or (f.shape[2] != self._chunk_3d_size[2]):
            raise RuntimeError("Make sure f is of the same size as in grid_partition!")

        f_y      = numpy.empty( self._chunk_y_size, dtype=numpy.float64, order='F' )
        d2fdy2_y = numpy.empty( self._chunk_y_size, dtype=numpy.float64, order='F' )

        self._grid_partition.transpose_3d_to_y( f, f_y )
        if self._order[1] == 6:
            raise NotImplementedError("6th order 2nd derivatives are not implemented yet. Sorry!")
        elif self._order[1] == 10:
            self._yder.d2d2(f_y, d2fdy2_y, self._chunk_y_size[0], self._chunk_y_size[2], bc1_=bc[0], bcn_=bc[1])
        self._grid_partition.transpose_y_to_3d( d2fdy2_y, d2fdy2 )


    def d2dz2(self, f, d2fdz2, bc=(0,0)):
        """
        Function to compute the Z second derivative of f

        f : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
            This array must be in the 3D decomposition
        d2fdz2 : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
               This array is in the 3D decomposition
        bc : integer tuple of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        if (f.shape[0] != self._chunk_3d_size[0]) or (f.shape[1] != self._chunk_3d_size[1]) or (f.shape[2] != self._chunk_3d_size[2]):
            raise RuntimeError("Make sure f is of the same size as in grid_partition!")

        f_z      = numpy.empty( self._chunk_z_size, dtype=numpy.float64, order='F' )
        d2fdz2_z = numpy.empty( self._chunk_z_size, dtype=numpy.float64, order='F' )

        self._grid_partition.transpose_3d_to_z( f, f_z )
        if self._order[2] == 6:
            raise NotImplementedError("6th order 2nd derivatives are not implemented yet. Sorry!")
        elif self._order[2] == 10:
            self._zder.d2d3(f_z, d2fdz2_z, self._chunk_z_size[0], self._chunk_z_size[1], bc1_=bc[0], bcn_=bc[1])
        self._grid_partition.transpose_z_to_3d( d2fdz2_z, d2fdz2 )


    def gradient(self, f, x_bc=(0,0), y_bc=(0,0), z_bc=(0,0)):
        """
        Function to compute the gradient of f

        f : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
            This array must be in the 3D decomposition
        dfd* : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
               This array is in the 3D decomposition
        *_bc : integer tuple of size 2 with the boundary condition at the left and right.
               0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        dfdx = numpy.empty( self._chunk_3d_size, dtype=numpy.float64, order='F' )
        dfdy = numpy.empty( self._chunk_3d_size, dtype=numpy.float64, order='F' )
        dfdz = numpy.empty( self._chunk_3d_size, dtype=numpy.float64, order='F' )

        self.ddx(f, dfdx, bc=x_bc)
        self.ddy(f, dfdy, bc=y_bc)
        self.ddz(f, dfdz, bc=z_bc)

        return dfdx, dfdy, dfdz


    def divergence(self, u, v, w, x_bc=(0,0), y_bc=(0,0), z_bc=(0,0)):
        """
        Function to compute the gradient of a vector (u,v,w)

        u, v, w : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
                  These arrays must be in the 3D decomposition
        divergence : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
                     This array is in the 3D decomposition
        *_bc : integer tuple of size 2 with the boundary condition at the left and right.
               0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        dudx = numpy.empty( self._chunk_3d_size, dtype=numpy.float64, order='F' )
        dvdy = numpy.empty( self._chunk_3d_size, dtype=numpy.float64, order='F' )
        dwdz = numpy.empty( self._chunk_3d_size, dtype=numpy.float64, order='F' )

        self.ddx(u, dudx, bc=x_bc)
        self.ddy(v, dvdy, bc=y_bc)
        self.ddz(w, dwdz, bc=z_bc)

        return (dudx + dvdy + dwdz)

