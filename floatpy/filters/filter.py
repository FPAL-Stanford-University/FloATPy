from mpi4py import MPI
import numpy
import pycf90
import pygaussian

from floatpy.parallel import _t3dmod

class Filter(object):
    """
    Class to perform parallel filter operations
    """

    def __init__(self, grid_partition, filter_type, periodic=(False, False, False)):
        """
        Constructor of the class

        grid_partition : t3d object or the grid_partition property of the parallel data reader class
        filter_type : string tuple of size 3 with the type of filter to use in each direction
                      options are "compact" and "gaussian"
        periodic : boolean tuple of size 3
        """

        if not isinstance(grid_partition, _t3dmod.t3d):
            raise RuntimeError("The given grid partition object is not an instance of the t3d class!")
        
        for i in range(3):
            if filter_type[i] not in ['compact','gaussian']:
                raise RuntimeError("filter_type[%d] has to be one of {'compact','gaussian'}" %i)

        self._filter_type = tuple(filter_type)
        self._grid_partition = grid_partition

        self._periodic = tuple(periodic)

        self._chunk_3d_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_3d_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._chunk_3d_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self._grid_partition.get_sz3d(self._chunk_3d_size)
        self._grid_partition.get_st3d(self._chunk_3d_lo)
        self._grid_partition.get_en3d(self._chunk_3d_hi)
        self._chunk_3d_lo = self._chunk_3d_lo - 1 # Convert to 0 based indexi
        self._chunk_3d_hi = self._chunk_3d_hi - 1 # Convert to 0 based indexi
        
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
        
        if self._filter_type[0] == 'compact':
            self._xfil = pycf90.cf90stuff.cf90( self._nx, self._periodic[0] )
        elif self._filter_type[0] == 'gaussian':
            self._xfil = pygaussian.gaussianstuff.gaussian( self._nx, self._periodic[0] )
        
        if self._filter_type[1] == 'compact':
            self._yfil = pycf90.cf90stuff.cf90( self._ny, self._periodic[1] )
        elif self._filter_type[1] == 'gaussian':
            self._yfil = pygaussian.gaussianstuff.gaussian( self._ny, self._periodic[1] )
        
        if self._filter_type[2] == 'compact':
            self._zfil = pycf90.cf90stuff.cf90( self._nz, self._periodic[2] )
        elif self._filter_type[2] == 'gaussian':
            self._zfil = pygaussian.gaussianstuff.gaussian( self._nz, self._periodic[2] )


    def filter_x(self, f, f_tilde, bc=(0,0)):
        """
        Function to filter the field f in the X direction

        f : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
            This array must be in the 3D decomposition
        f_tilde : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
                  This array is in the 3D decomposition
        bc : integer tuple of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        if (f.shape[0] != self._chunk_3d_size[0]) or (f.shape[1] != self._chunk_3d_size[1]) or (f.shape[2] != self._chunk_3d_size[2]):
            raise RuntimeError("Make sure f is of the same size as in grid_partition!")
        
        f_x       = numpy.empty( self._chunk_x_size, dtype=numpy.float64, order='F' )
        f_tilde_x = numpy.empty( self._chunk_x_size, dtype=numpy.float64, order='F' )
        
        self._grid_partition.transpose_3d_to_x( f, f_x )
        self._xfil.filter1(f_x, f_tilde_x, self._chunk_x_size[1], self._chunk_x_size[2], bc1_=bc[0], bcn_=bc[1])
        self._grid_partition.transpose_x_to_3d( f_tilde_x, f_tilde )


    def filter_y(self, f, f_tilde, bc=(0,0)):
        """
        Function to filter the field f in the Y direction

        f : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
            This array must be in the 3D decomposition
        f_tilde : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
                  This array is in the 3D decomposition
        bc : integer tuple of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        if (f.shape[0] != self._chunk_3d_size[0]) or (f.shape[1] != self._chunk_3d_size[1]) or (f.shape[2] != self._chunk_3d_size[2]):
            raise RuntimeError("Make sure f is of the same size as in grid_partition!")
        
        f_y       = numpy.empty( self._chunk_y_size, dtype=numpy.float64, order='F' )
        f_tilde_y = numpy.empty( self._chunk_y_size, dtype=numpy.float64, order='F' )
        
        self._grid_partition.transpose_3d_to_y( f, f_y )
        self._yfil.filter2(f_y, f_tilde_y, self._chunk_y_size[0], self._chunk_y_size[2], bc1_=bc[0], bcn_=bc[1])
        self._grid_partition.transpose_y_to_3d( f_tilde_y, f_tilde )


    def filter_z(self, f, f_tilde, bc=(0,0)):
        """
        Function to filter the field f in the Z direction

        f : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
            This array must be in the 3D decomposition
        f_tilde : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
                  This array is in the 3D decomposition
        bc : integer tuple of size 2 with the boundary condition at the left and right.
             0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        if (f.shape[0] != self._chunk_3d_size[0]) or (f.shape[1] != self._chunk_3d_size[1]) or (f.shape[2] != self._chunk_3d_size[2]):
            raise RuntimeError("Make sure f is of the same size as in grid_partition!")
        
        f_z       = numpy.empty( self._chunk_z_size, dtype=numpy.float64, order='F' )
        f_tilde_z = numpy.empty( self._chunk_z_size, dtype=numpy.float64, order='F' )
        
        self._grid_partition.transpose_3d_to_z( f, f_z )
        self._zfil.filter3(f_z, f_tilde_z, self._chunk_z_size[0], self._chunk_z_size[1], bc1_=bc[0], bcn_=bc[1])
        self._grid_partition.transpose_z_to_3d( f_tilde_z, f_tilde )

    def filter(self, f, x_bc=(0,0), y_bc=(0,0), z_bc=(0,0), ntimes=1):
        """
        Function to filter the field f in all direction

        f : input 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
            This array must be in the 3D decomposition
        f_tilde : output 3D numpy array in Fortran contiguous layout with the 1st index being X and last being Z
                  This array is in the 3D decomposition
        *_bc : integer tuple of size 2 with the boundary condition at the left and right.
               0 is general, 1 is symmetric, -1 is anti-symmetric. Only required if non-periodic
        """

        tmp     = numpy.empty( self._chunk_3d_size, dtype=numpy.float64, order='F' )
        f_tilde = numpy.empty( self._chunk_3d_size, dtype=numpy.float64, order='F' )

        self.filter_x(f, f_tilde, bc=x_bc)
        for i in range(ntimes-1):
            tmp[:,:,:] = f_tilde[:,:,:]
            self.filter_x(tmp, f_tilde, bc=x_bc)

        self.filter_y(f_tilde, tmp, bc=y_bc)
        for i in range(ntimes-1):
            f_tilde[:,:,:] = tmp[:,:,:]
            self.filter_y(f_tilde, tmp, bc=y_bc)

        self.filter_z(tmp, f_tilde, bc=z_bc)
        for i in range(ntimes-1):
            tmp[:,:,:] = f_tilde[:,:,:]
            self.filter_z(tmp, f_tilde, bc=z_bc)

        return f_tilde


