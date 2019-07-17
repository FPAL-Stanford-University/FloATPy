import numpy as np
import pyfftw
from scipy.sparse import csr_matrix, lil_matrix, SparseEfficiencyWarning
from scipy.sparse.linalg import spsolve, lgmres
import scipy
import scipy.io as sio
import warnings
from numpy.linalg import det


class PoissonSol:
    def __init__ (self,numSet):
        """
        Initialize Poisson solver.
        args:
            numSet - object of NumSetting
        """
        warnings.simplefilter("ignore", SparseEfficiencyWarning)
        self.numSet = numSet
        self.cout = numSet.cout
        self.NX = numSet.NX
        self.NY = numSet.NY
        self.NZ = numSet.NZ
        self.LX = numSet.LX
        self.LY = numSet.LY
        self.LZ = numSet.LZ
        self.size_3d = numSet.chunk_3d_size
        self.size_x = numSet.chunk_x_size
        self.size_y = numSet.chunk_y_size
        self.size_z = numSet.chunk_z_size
        self.chunk_x_ind = range(numSet.chunk_z_lo[0],numSet.chunk_z_hi[0]+1)
        self.chunk_y_ind = range(numSet.chunk_z_lo[1],numSet.chunk_z_hi[1]+1)
        self.chunk_z_ind = range(numSet.chunk_z_lo[2],numSet.chunk_z_hi[2]+1)
        self.dx = numSet.dx
        self.dy = numSet.dy
        self.dz = numSet.dz
        self._partition = numSet.grid_partition
        self._initFFTW()
        if numSet.order[2] == 6:
            self._initMAT4CompactDiff_6()
        elif numSet.order[2] == 10:
            self._initMAT4CompactDiff_10()
        else:
            raise RuntimeError("Order of accuracy must be either 6 or 10.")


    def _initFFTW(self):
        """
        Initialize all settings for pyFFTW.
        """
        f2dx = pyfftw.empty_aligned(self.size_x[[0,1]], dtype='complex128')
        F2dx = pyfftw.empty_aligned(self.size_x[[0,1]], dtype='complex128')
        self.fft2dx = pyfftw.FFTW(input_array=f2dx, output_array=F2dx, axes=(0,), 
                direction='FFTW_FORWARD')
        f2dy = pyfftw.empty_aligned(self.size_y[[0,1]], dtype='complex128')
        F2dy = pyfftw.empty_aligned(self.size_y[[0,1]], dtype='complex128')
        self.fft2dy = pyfftw.FFTW(input_array=f2dy, output_array=F2dy, axes=(1,), 
                direction='FFTW_FORWARD')

        fx = pyfftw.empty_aligned(self.size_x, dtype='complex128')
        Fx = pyfftw.empty_aligned(self.size_x, dtype='complex128')
        self.fftx = pyfftw.FFTW(input_array=fx, output_array=Fx, axes=(0,), 
                direction='FFTW_FORWARD')
        fy = pyfftw.empty_aligned(self.size_y, dtype='complex128')
        Fy = pyfftw.empty_aligned(self.size_y, dtype='complex128')
        self.ffty = pyfftw.FFTW(input_array=fy, output_array=Fy, axes=(1,), 
                direction='FFTW_FORWARD')
        fz = pyfftw.empty_aligned(self.size_z, dtype='complex128')
        Fz = pyfftw.empty_aligned(self.size_z, dtype='complex128')
        self.fftz = pyfftw.FFTW(input_array=fz, output_array=Fz, axes=(2,), 
                direction='FFTW_FORWARD')

        gx = pyfftw.empty_aligned(self.size_x, dtype='complex128')
        Gx = pyfftw.empty_aligned(self.size_x, dtype='complex128')
        self.ifftx = pyfftw.FFTW(input_array=Gx, output_array=gx, axes=(0,), 
                direction='FFTW_BACKWARD')
        gy = pyfftw.empty_aligned(self.size_y, dtype='complex128')
        Gy = pyfftw.empty_aligned(self.size_y, dtype='complex128')
        self.iffty = pyfftw.FFTW(input_array=Gy, output_array=gy, axes=(1,), 
                direction='FFTW_BACKWARD')
        gz = pyfftw.empty_aligned(self.size_z, dtype='complex128')
        Gz = pyfftw.empty_aligned(self.size_z, dtype='complex128')
        self.ifftz = pyfftw.FFTW(input_array=Gz, output_array=gz, axes=(2,), 
                direction='FFTW_BACKWARD')

        kx_x = np.arange(-self.NX//2, self.NX//2, dtype=np.float64) * (2*np.pi/self.LX)
        ky_y = np.arange(-self.NY//2, self.NY//2, dtype=np.float64) * (2*np.pi/self.LY)
        kz_z = np.arange(-self.NZ//2, self.NZ//2, dtype=np.float64) * (2*np.pi/self.LZ)
        kx_x = np.fft.fftshift(kx_x)
        ky_y = np.fft.fftshift(ky_y)
        kz_z = np.fft.fftshift(kz_z)
        self.kx = kx_x[self.chunk_x_ind]
        self.ky = ky_y[self.chunk_y_ind]
        self.kz = kz_z[self.chunk_z_ind]

    def _initMAT4CompactDiff_6(self):
        """
        Initialize matrices for 6th-order compact finite difference method:
        Af''= (1/h^2)Bf
        in z direction.
        return:
            A - matrix A
            B - matrix B
        """
        NZ = self.size_z[2]
        alpi = 2.0 / 11.0
        qi = 12.0 / 11.0
        ri = 3.0 / 44.0
        alp1 = 1.0 / 10.0
        q1 = 6.0 / 5.0
        alp0 = 11.0
        p = 13.0
        q = -27.0
        r = 15.0
        s = -1.0
        #ind = range(2,NZ-3)
        A = lil_matrix((NZ,NZ),dtype=np.float64)
        B = lil_matrix((NZ,NZ),dtype=np.float64)

        A.setdiag(1.0, k=0)
        A.setdiag(alpi, k=1)
        A.setdiag(alpi, k=-1)
        A[0,1], A[NZ-1,NZ-2] = alp0, alp0
        A[1,[0,2]], A[NZ-2,[NZ-1,NZ-3]] = alp1, alp1

        B.setdiag(-2.0*(qi+ri), k=0)
        B.setdiag(qi, k=1)
        B.setdiag(qi, k=-1)
        B.setdiag(ri, k=2)
        B.setdiag(ri, k=-2)
        B[0,:] = B[1,:] = B[NZ-1,:] = B[NZ-2,:] = 0.0 # clear the boundary and near boundary rows
        B[0,[0,1,2,3]] = [p,q,r,s]
        B[NZ-1,[NZ-4,NZ-3,NZ-2,NZ-1]] = [s,r,q,p]
        B[1,[0,1,2]] = [q1, -2*q1, q1]
        B[NZ-2,[NZ-3,NZ-2,NZ-1]] = [q1, -2*q1, q1]

        self.A = csr_matrix(A)
        self.B = csr_matrix(B)
        self.A.eliminate_zeros()
        self.B.eliminate_zeros()

        derCoef = [-1.5, 2.0, -0.5, 0.0]
        self.bcMAT = {  # First and last row modifications for matrix
            'DL':[1.0, 0.0, 0.0, 0.0],
            'DH':[0.0, 0.0, 0.0, 1.0],
            'NL':[x/self.dz for x in derCoef],
            'NH':[-x/self.dz for x in reversed(derCoef)]
        }


    def _initMAT4CompactDiff_10(self):
        """
        Initialize matrices for 6th-order compact finite difference method:
        Af''= (1/h^2)Bf
        in z direction.
        return:
            A - matrix A
            B - matrix B
        """
        NZ = self.size_z[2]
        alpi = 334.0 / 899.0
        beti = 43.0 / 1798.0
        qi = 1065.0 / 1798.0
        ri = 1038.0 / 899.0 / 4.0
        si = 79.0 / 1798.0 / 9.0
        alp1 = 1.0 / 10.0
        q1 = 6.0 / 5.0
        alp2 = 344.0 / 1179.0
        bet2 = 23.0 / 2358.0
        q2 = 320.0 / 393.0
        r2 = 310.0 / 393.0 / 4.0
        alp0 = 11.0
        p = 13.0
        q = -27.0
        r = 15.0
        s = -1.0
        #alp0, p, q, r, s = 0., 0., 0., 0., 0.
        #ind = range(2,NZ-3)
        A = lil_matrix((NZ,NZ),dtype=np.float64)
        B = lil_matrix((NZ,NZ),dtype=np.float64)

        A.setdiag(alpi, k=1)
        A.setdiag(alpi, k=-1)
        A.setdiag(beti, k=2)
        A.setdiag(beti, k=-2)
        A[0,:] = A[1,:] = A[2,:] = A[NZ-1,:] = A[NZ-2,:] = A[NZ-3,:] = 0.0
        A.setdiag(1.0, k=0)
        A[0,1], A[NZ-1,NZ-2] = alp0, alp0
        A[1,[0,2]], A[NZ-2,[NZ-1,NZ-3]] = alp1, alp1
        A[2,[1,3]], A[NZ-3,[NZ-2,NZ-4]], A[2,[0,4]], A[NZ-3,[NZ-1,NZ-5]] = alp2, alp2, bet2, bet2

        B.setdiag(-2.0*(qi+ri+si), k=0)
        B.setdiag(qi, k=1)
        B.setdiag(qi, k=-1)
        B.setdiag(ri, k=2)
        B.setdiag(ri, k=-2)
        B.setdiag(si, k=3)
        B.setdiag(si, k=-3)
        B[0,:] = B[1,:] = B[NZ-1,:] = B[NZ-2,:] = B[2,:] = B[NZ-3,:] = 0.0 # clear the boundary and near boundary rows
        B[0,[0,1,2,3]] = [p,q,r,s]
        B[NZ-1,[NZ-4,NZ-3,NZ-2,NZ-1]] = [s,r,q,p]
        B[1,[0,1,2]] = [q1, -2*q1, q1]
        B[NZ-2,[NZ-3,NZ-2,NZ-1]] = [q1, -2*q1, q1]
        B[2,[0,1,2,3,4]] = [r2, q2, -2*(r2+q2), q2, r2]
        B[NZ-3,[NZ-1,NZ-2,NZ-3,NZ-4,NZ-5]] = [r2, q2, -2*(r2+q2), q2, r2]

        self.A = csr_matrix(A)
        self.B = csr_matrix(B)
        self.A.eliminate_zeros()
        self.B.eliminate_zeros()

        derCoef = [-11.0/6.0, 3.0, -3.0/2.0, 1.0/3.0]
        self.bcMAT = {  # First and last row modifications for matrix
            'DL':[1.0, 0.0, 0.0, 0.0],
            'DH':[0.0, 0.0, 0.0, 1.0],
            'NL':[x/self.dz for x in derCoef],
            'NH':[-x/self.dz for x in reversed(derCoef)]
        }

    def _fftX(self,f3d):
        """
        Fourier transform in x direction.
        args:
            f3d - field with chunk_3d_size to be transformed
        return:
            F3d - complex128 fields after transform with chunk_3d_size
        """
        fxR = np.empty(self.size_x, dtype=np.dtype(np.float64, align=True), order='F')
        fxI = np.empty(self.size_x, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_3d_to_x(np.real(f3d),fxR)
        self._partition.transpose_3d_to_x(np.imag(f3d),fxI)

        FxC = pyfftw.empty_aligned(self.size_x, dtype='complex128')
        fxC = pyfftw.empty_aligned(self.size_x, dtype='complex128')
        fxC[:] = fxR + 1j * fxI
        self.fftx(input_array=fxC, output_array=FxC)

        F3dR = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        F3dI = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_x_to_3d(np.real(FxC), F3dR)
        self._partition.transpose_x_to_3d(np.imag(FxC), F3dI)

        F3d = np.empty(self.size_3d, dtype=np.dtype(np.complex128, align=True), order='F')
        F3d[:] = F3dR + 1j * F3dI
        return F3d



    def _fftY(self,f3d):
        """
        Fourier transform in y direction.
        args:
            f3d - field with chunk_3d_size to be transformed
        return:
            F3d - complex128 fields after transform with chunk_3d_size
        """
        fyR = np.empty(self.size_y, dtype=np.dtype(np.float64, align=True), order='F')
        fyI = np.empty(self.size_y, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_3d_to_y(np.real(f3d),fyR)
        self._partition.transpose_3d_to_y(np.imag(f3d),fyI)

        FyC = pyfftw.empty_aligned(self.size_y, dtype='complex128')
        fyC = pyfftw.empty_aligned(self.size_y, dtype='complex128')
        fyC[:] = fyR + 1j * fyI
        self.ffty(input_array=fyC, output_array=FyC)

        F3dR = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        F3dI = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_y_to_3d(np.real(FyC), F3dR)
        self._partition.transpose_y_to_3d(np.imag(FyC), F3dI)

        F3d = np.empty(self.size_3d, dtype=np.dtype(np.complex128, align=True), order='F')
        F3d[:] = F3dR + 1j * F3dI
        return F3d


    def _fftZ(self,f3d):
        """
        Fourier transform in z direction.
        args:
            f3d - field with chunk_3d_size to be transformed
        return:
            F3d - complex128 fields after transform with chunk_3d_size
        """
        fzR = np.empty(self.size_z, dtype=np.dtype(np.float64, align=True), order='F')
        fzI = np.empty(self.size_z, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_3d_to_z(np.real(f3d),fzR)
        self._partition.transpose_3d_to_z(np.imag(f3d),fzI)

        FzC = pyfftw.empty_aligned(self.size_z, dtype='complex128')
        fzC = pyfftw.empty_aligned(self.size_z, dtype='complex128')
        fzC[:] = fzR + 1j * fzI
        self.fftz(input_array=fzC, output_array=FzC)

        F3dR = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        F3dI = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_z_to_3d(np.real(FzC), F3dR)
        self._partition.transpose_z_to_3d(np.imag(FzC), F3dI)

        F3d = np.empty(self.size_3d, dtype=np.dtype(np.complex128, align=True), order='F')
        F3d[:] = F3dR + 1j * F3dI
        return F3d


    def _fftXY(self,f3d):
        """
        2D Fourier transform in xy plane.
        args:
            f3d - fields that with chunk_3d_size to be transformed
        return:
            F3d - complex128 fields after transform with chunk_3d_size
        """
        F3dX = np.empty(self.size_3d, dtype=np.dtype(np.complex128, align=True), order='F')
        F3d = np.empty(self.size_3d, dtype=np.dtype(np.complex128, align=True), order='F')
        F3dX[:] = self._fftX(f3d)
        F3d[:] = self._fftY(F3dX)
        return F3d



    def _ifftX(self,F3d):
        """
        Inverse Fourier transform in x direction.
        args:
            F3d - FT field with chunk_3d_size
        return:
            f3d - IFT field with chunk_3d_size
        """
        FxR = np.empty(self.size_x, dtype=np.dtype(np.float64, align=True), order='F')
        FxI = np.empty(self.size_x, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_3d_to_x(np.real(F3d), FxR)
        self._partition.transpose_3d_to_x(np.imag(F3d), FxI)

        FxC = pyfftw.empty_aligned(self.size_x, dtype='complex128')
        fxC = pyfftw.empty_aligned(self.size_x, dtype='complex128')
        FxC[:] = FxR + 1j * FxI
        self.ifftx(input_array=FxC, output_array=fxC)

        f3dR = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        f3dI = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_x_to_3d(np.real(fxC), f3dR)
        self._partition.transpose_x_to_3d(np.imag(fxC), f3dI)

        f3d = np.empty(self.size_3d, dtype=np.dtype(np.complex128, align=True), order='F')
        f3d[:] = f3dR + 1j * f3dI
        return f3d



    def _ifftY(self,F3d):
        """
        Inverse Fourier transform in y direction.
        args:
            F3d - FT field with chunk_3d_size
        return:
            f3d - IFT field with chunk_3d_size
        """
        FyR = np.empty(self.size_y, dtype=np.dtype(np.float64, align=True), order='F')
        FyI = np.empty(self.size_y, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_3d_to_y(np.real(F3d), FyR)
        self._partition.transpose_3d_to_y(np.imag(F3d), FyI)
        FyC = pyfftw.empty_aligned(self.size_y, dtype='complex128')
        fyC = pyfftw.empty_aligned(self.size_y, dtype='complex128')
        FyC[:] = FyR + 1j * FyI
        self.iffty(input_array=FyC, output_array=fyC)

        f3dR = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        f3dI = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_y_to_3d(np.real(fyC), f3dR)
        self._partition.transpose_y_to_3d(np.imag(fyC), f3dI)

        f3d = np.empty(self.size_3d, dtype=np.dtype(np.complex128, align=True), order='F')
        f3d[:] = f3dR + 1j * f3dI
        return f3d


    def _ifftZ(self,F3d):
        """
        Inverse Fourier transform in z direction.
        args:
            F3d - FT field with chunk_3d_size
        return:
            f3d - IFT field with chunk_3d_size
        """
        FzR = np.empty(self.size_z, dtype=np.dtype(np.float64, align=True), order='F')
        FzI = np.empty(self.size_z, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_3d_to_z(np.real(F3d), FzR)
        self._partition.transpose_3d_to_z(np.imag(F3d), FzI)
        FzC = pyfftw.empty_aligned(self.size_z, dtype='complex128')
        fzC = pyfftw.empty_aligned(self.size_z, dtype='complex128')
        FzC[:] = FzR + 1j * FzI
        self.ifftz(input_array=FzC, output_array=fzC)

        f3dR = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        f3dI = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_z_to_3d(np.real(fzC), f3dR)
        self._partition.transpose_z_to_3d(np.imag(fzC), f3dI)

        f3d = np.empty(self.size_3d, dtype=np.dtype(np.complex128, align=True), order='F')
        f3d[:] = f3dR + 1j * f3dI
        return f3d

    def _ifftXY(self,F3dXY):
        """
        2D inverse Fourier transform in xy plane.
        args:
            F3dXY - 2D FT field with chunk_3d_size
        return:
            f3dR - 2D IFT field with chunk_3d_size and real data
        """
        F3dX = self._ifftY(F3dXY)
        f3dC = self._ifftX(F3dX)
        f3dR = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        f3dR[:] = np.real(f3dC)
        return f3dR



    def _fft2dX(self,fxR,fxI):
        """
        Fourier transform in x direction for 2D array.
        args:
            # f_chunk - field with chunk_3d_size[[0,1]] to be transformed
            fxR - real part of the field with size_x[[0,1]] to be transformed
            fxI - imag part of the field with size_x[[0,1]] to be transformed
        return:
            # F_chunk - complex128 fields after transform with chunk_3d_size[[0,1]]
            FxC
        """
        FxC = pyfftw.empty_aligned(self.size_x[[0,1]], dtype='complex128')
        fxC = pyfftw.empty_aligned(self.size_x[[0,1]], dtype='complex128')
        fxC[:] = fxR + 1j * fxI
        self.fft2dx(input_array=fxC, output_array=FxC)

        return FxC



    def _fft2dY(self,fyR, fyI):
        """
        Fourier transform in y direction for 2D array.
        args:
            fyR - real part of the field with size_y[[0,1]] to be transformed
            fyI - imag part of the field with size_y[[0,1]] to be transformed
        return:
            FyC - complex128 fields after transform with chunk_3d_size[[0,1]]
        """
        FyC = pyfftw.empty_aligned(self.size_y[[0,1]], dtype='complex128')
        fyC = pyfftw.empty_aligned(self.size_y[[0,1]], dtype='complex128')
        fyC[:] = fyR + 1j * fyI
        self.fft2dy(input_array=fyC, output_array=FyC)

        return FyC


    def _fftXY2Layer(self,f_chunk_l1, f_chunk_l2):
        """
        2D FFT transform for 2D np array on xy plane.
        args:
            f_chunk_l1 - 2D array (layer 1 of the chunk)
            f_chunk_l2 - 2D array (layer 2 of the chunk)
        return:
            F_chunk_l1
            F_chunk_l2
        """

        """
        fC_chunk = np.empty(self.size_3d, dtype=np.dtype(np.complex128, align=True), order='F')
        fR_chunk = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        fI_chunk = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')

        # FFT in x direction
        fxR3d = np.empty(self.size_x, dtype=np.dtype(np.float64, align=True), order='F')
        fxI3d = np.zeros(self.size_x, dtype=np.dtype(np.float64, align=True), order='F')
        fR_chunk[:,:,0] = f_chunk_l1
        fR_chunk[:,:,-1] = f_chunk_l2
        self._partition.transpose_3d_to_x(fR_chunk,fxR3d)

        fxR1 = np.empty(self.size_x[[0,1]], dtype=np.dtype(np.float64, align=True), order='F')
        fxI1 = np.zeros(self.size_x[[0,1]], dtype=np.dtype(np.float64, align=True), order='F')
        fxR2 = np.empty(self.size_x[[0,1]], dtype=np.dtype(np.float64, align=True), order='F')
        fxI2 = np.zeros(self.size_x[[0,1]], dtype=np.dtype(np.float64, align=True), order='F')
        fxR1[:] = fxR3d[:,:,0]
        fxR2[:] = fxR3d[:,:,-1]
        FxC1 = self._fft2dX(fxR1,fxI1)
        FxC2 = self._fft2dX(fxR2,fxI2)
        fxR3d[:,:,0] = np.real(FxC1)
        fxI3d[:,:,0] = np.imag(FxC1)
        fxR3d[:,:,-1] = np.real(FxC2)
        fxI3d[:,:,-1] = np.imag(FxC2)
        self._partition.transpose_x_to_3d(fxR3d,fR_chunk)
        self._partition.transpose_x_to_3d(fxI3d,fI_chunk)


        # FFT in y direction
        fyR3d = np.empty(self.size_y, dtype=np.dtype(np.float64, align=True), order='F')
        fyI3d = np.empty(self.size_y, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_3d_to_y(fR_chunk, fyR3d)
        self._partition.transpose_3d_to_y(fI_chunk, fyI3d)

        fyR1 = np.empty(self.size_y[[0,1]], dtype=np.dtype(np.float64, align=True), order='F')
        fyI1 = np.empty(self.size_y[[0,1]], dtype=np.dtype(np.float64, align=True), order='F')
        fyR2 = np.empty(self.size_y[[0,1]], dtype=np.dtype(np.float64, align=True), order='F')
        fyI2 = np.empty(self.size_y[[0,1]], dtype=np.dtype(np.float64, align=True), order='F')
        fyR1[:] = fyR3d[:,:,0]
        fyI1[:] = fyI3d[:,:,0]
        fyR2[:] = fyR3d[:,:,-1]
        fyI2[:] = fyI3d[:,:,-1]
        FyC1 = self._fft2dY(fyR1,fyI1)
        FyC2 = self._fft2dY(fyR2,fyI2)
        fyR3d[:,:,0] = np.real(FyC1)
        fyI3d[:,:,0] = np.imag(FyC1)
        fyR3d[:,:,-1] = np.real(FyC2)
        fyI3d[:,:,-1] = np.imag(FyC2)
        self._partition.transpose_y_to_3d(fyR3d,fR_chunk)
        self._partition.transpose_y_to_3d(fyI3d,fI_chunk)
        """

        fC_chunk = np.empty(self.size_3d, dtype=np.dtype(np.complex128, align=True), order='F')
        fR_chunk = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        fI_chunk = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')

        fC_chunk[:]

        # return values
        F_chunk_l1 = np.empty(self.size_z[[0,1]], dtype=np.dtype(np.complex128, align=True), order='F')
        F_chunk_l2 = np.empty(self.size_z[[0,1]], dtype=np.dtype(np.complex128, align=True), order='F')
        F_chunk_l1[:] = fR_chunk[:,:,0] + 1j * fI_chunk[:,:,0]
        F_chunk_l2[:] = fR_chunk[:,:,-1] + 1j * fI_chunk[:,:,-1]


        return F_chunk_l1, F_chunk_l2

    def _fftXYLayer(self,f_chunk_l):
        """
        2D FFT transform for 2D np array on xy plane.
        args:
            f_chunk_l - 2D array (layer of the chunk)

        return:
            F_chunk_l - 2D array (FFTXYed layer of the size_z)
        """
        f3d = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        for k in range(self.size_3d[2]):
            f3d[:,:,k] = f_chunk_l
        F3d = self._fftXY(f3d)

        buffer_3d = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        buffer_z  = np.empty(self.size_z, dtype=np.dtype(np.float64, align=True), order='F')

        buffer_3d[:] = np.real(F3d)
        self._partition.transpose_3d_to_z(buffer_3d, buffer_z)
        F_chunk_l = np.empty(self.size_z[[0,1]], dtype=np.dtype(np.complex128, align=True), order='F')
        F_chunk_l[:,:] = buffer_z[:,:,0]

        buffer_3d[:] = np.imag(F3d)
        self._partition.transpose_3d_to_z(buffer_3d, buffer_z)
        F_chunk_l[:,:] = F_chunk_l + 1j * buffer_z[:,:,0]

        return F_chunk_l


    def _solveZ(self, rhs_3d, bcLow, bcHigh, bcLowType='D', bcHighType='D'):
        """
        Solve 2nd order ODE for each vertical column in FFTXY space.
        args:
            rhs_3d - right-hand-side with chunk_3d_size in FFTXY space
            bcLow - boundary value at z=ZMIN in FT space
            bcHigh - boundary value at z=ZMAX in FT space
            bcLowType - boundary condition type at Z=ZMIN: 'D' for Dirichlet, 'N' for Neumann
            bcHighType - boundary condition type at Z=ZMAX: 'D' for Dirichlet, 'N' for Neumann
        return:
            solnC_3d - solution in FFTXY space with complex128 dtype and chunk_3d_size
        """
        rhsC_z = np.empty(self.size_z, dtype=np.dtype(np.complex128, align=True), order='F')
        rhsR_z = np.empty(self.size_z, dtype=np.dtype(np.float64, align=True), order='F')
        rhsI_z = np.empty(self.size_z, dtype=np.dtype(np.float64, align=True), order='F')
        #rhs_3d_buffer= np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')

        #rhs_3d_buffer[:] =
        self._partition.transpose_3d_to_z(np.real(rhs_3d), rhsR_z)
        #rhs_3d_buffer[:] =
        self._partition.transpose_3d_to_z(np.imag(rhs_3d), rhsI_z)
        rhsC_z[:] = rhsR_z + 1j * rhsI_z

        NZ = self.size_z[2]
        solnC_z = np.empty(self.size_z, dtype=np.dtype(np.complex128, align=True), order='F')
        for j in range(self.size_z[1]):
            for i in range(self.size_z[0]):
                # Build sparse linear system for interior mesh points
                kxy2 = self.kx[i]**2 + self.ky[j]**2
                M = self.B / (self.dz**2) - self.A * kxy2
                M = csr_matrix(M, dtype=np.complex128)
                #rhs_vec = np.zeros(NZ, dtype=np.dtype(np.complex128, align=True), order='F')
                rhs_vec = self.A.dot(rhsC_z[i,j,:])
                # Impose boundary conditions on bottom and top
                M[0,[0,1,2,3]] = self.bcMAT[bcLowType+'L']
                M[NZ-1,[NZ-4,NZ-3,NZ-2,NZ-1]] = self.bcMAT[bcHighType+'H']
                rhs_vec[0] = bcLow[i,j]
                rhs_vec[NZ-1] = bcHigh[i,j]
                M.eliminate_zeros()
                # Solve the system in z direction
                solnC_z[i,j,:] = spsolve(M,rhs_vec)

        solnR_3d = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        solnI_3d = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_z_to_3d(np.real(solnC_z), solnR_3d)
        self._partition.transpose_z_to_3d(np.imag(solnC_z), solnI_3d)
        solnC_3d = np.empty(self.size_3d, dtype=np.dtype(np.complex128, align=True), order='F')
        solnC_3d[:] = solnR_3d + 1j * solnI_3d
        return solnC_3d



    def solve(self, rhs, bcLow, bcHigh, bcLowType='D', bcHighType='D'):
        """
        Completely solve Poisson equation for periodic domain in x and y directions.
        args:
            rhs - right-hand-side of the Poisson equation with chunk_3d_size
            bcLow - boundary value at z=ZMIN in physical space
            bcHigh - boundary value at z=ZMAX in physical space
            bcLowType - boundary condition type at Z=ZMIN: 'D' for Dirichlet, 'N' for Neumann
            bcHighType - boundary condition type at Z=ZMAX: 'D' for Dirichlet, 'N' for Neumann
        """
        BC_Low =  self._fftXYLayer(bcLow)
        BC_High=  self._fftXYLayer(bcHigh)

        RHS = self._fftXY(rhs)
        solnC = self._solveZ(RHS, bcLow=BC_Low, bcHigh=BC_High, bcLowType=bcLowType, bcHighType=bcHighType)
        return self._ifftXY(solnC)


    def transposeTest(self, field_3d):
        field_z = np.empty(self.size_z, dtype=np.dtype(np.float64, align=True), order='F')
        z = np.linspace(self.numSet.ZMIN, self.numSet.ZMAX, num=self.numSet.NZ)
        for k in range(self.numSet.NZ):
            field_z[:,:,k] = z[k]**2
        field_3d_new = np.empty(self.size_3d, dtype=np.dtype(np.float64, align=True), order='F')
        self._partition.transpose_z_to_3d(np.absolute(field_z),field_3d_new)
        err = np.max(np.abs(field_3d - field_3d_new))
        self.cout("err = {:5E}".format(err))
