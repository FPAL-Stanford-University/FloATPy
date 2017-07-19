from __future__ import print_function, absolute_import, division
import _pyt3d
import f90wrap.runtime
import logging

class T3Dmod(f90wrap.runtime.FortranModule):
    """
    Module t3dmod
    
    
    Defined at t3dMod.F90 lines 1-1577
    
    """
    class t3d(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t3d)
        
        
        Defined at t3dMod.F90 lines 18-89
        
        """
        def __init__(self, comm3d, nx, ny, nz, px, py, pz, periodic_, reorder, fail, \
            nghosts=None, createcrosscommunicators=None, handle=None):
            """
            self = T3D(comm3d, nx, ny, nz, px, py, pz, periodic_, reorder, fail[, nghosts, \
                createcrosscommunicators])
            
            
            Defined at t3dMod.F90 lines 97-505
            
            Parameters
            ----------
            comm3d : int
            nx : int
            ny : int
            nz : int
            px : int
            py : int
            pz : int
            periodic_ : bool array
            reorder : bool
            fail : bool
            nghosts : int array
            createcrosscommunicators : bool
            
            Returns
            -------
            this : T3D
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            self._handle = _pyt3d.f90wrap_init(comm3d=comm3d, nx=nx, ny=ny, nz=nz, px=px, \
                py=py, pz=pz, periodic_=periodic_, reorder=reorder, fail=fail, \
                nghosts=nghosts, createcrosscommunicators=createcrosscommunicators)
        
        def __del__(self):
            """
            Destructor for class T3D
            
            
            Defined at t3dMod.F90 lines 507-523
            
            Parameters
            ----------
            this : T3D
            
            """
            if self._alloc:
                _pyt3d.f90wrap_destroy(this=self._handle)
        
        def transpose_3d_to_x(self, input, output):
            """
            transpose_3d_to_x(self, input, output)
            
            
            Defined at t3dMod.F90 lines 525-577
            
            Parameters
            ----------
            this : T3D
            input : float array
            output : float array
            
            """
            _pyt3d.f90wrap_transpose_3d_to_x(this=self._handle, input=input, output=output)
        
        def transpose_x_to_3d(self, input, output):
            """
            transpose_x_to_3d(self, input, output)
            
            
            Defined at t3dMod.F90 lines 579-620
            
            Parameters
            ----------
            this : T3D
            input : float array
            output : float array
            
            """
            _pyt3d.f90wrap_transpose_x_to_3d(this=self._handle, input=input, output=output)
        
        def transpose_3d_to_y(self, input, output):
            """
            transpose_3d_to_y(self, input, output)
            
            
            Defined at t3dMod.F90 lines 622-674
            
            Parameters
            ----------
            this : T3D
            input : float array
            output : float array
            
            """
            _pyt3d.f90wrap_transpose_3d_to_y(this=self._handle, input=input, output=output)
        
        def transpose_y_to_3d(self, input, output):
            """
            transpose_y_to_3d(self, input, output)
            
            
            Defined at t3dMod.F90 lines 676-718
            
            Parameters
            ----------
            this : T3D
            input : float array
            output : float array
            
            """
            _pyt3d.f90wrap_transpose_y_to_3d(this=self._handle, input=input, output=output)
        
        def transpose_3d_to_z(self, input, output):
            """
            transpose_3d_to_z(self, input, output)
            
            
            Defined at t3dMod.F90 lines 720-772
            
            Parameters
            ----------
            this : T3D
            input : float array
            output : float array
            
            """
            _pyt3d.f90wrap_transpose_3d_to_z(this=self._handle, input=input, output=output)
        
        def transpose_z_to_3d(self, input, output):
            """
            transpose_z_to_3d(self, input, output)
            
            
            Defined at t3dMod.F90 lines 774-816
            
            Parameters
            ----------
            this : T3D
            input : float array
            output : float array
            
            """
            _pyt3d.f90wrap_transpose_z_to_3d(this=self._handle, input=input, output=output)
        
        def fill_halo_x(self, array):
            """
            fill_halo_x(self, array)
            
            
            Defined at t3dMod.F90 lines 1164-1240
            
            Parameters
            ----------
            this : T3D
            array : float array
            
            """
            _pyt3d.f90wrap_fill_halo_x(this=self._handle, array=array)
        
        def fill_halo_y(self, array):
            """
            fill_halo_y(self, array)
            
            
            Defined at t3dMod.F90 lines 1242-1262
            
            Parameters
            ----------
            this : T3D
            array : float array
            
            """
            _pyt3d.f90wrap_fill_halo_y(this=self._handle, array=array)
        
        def fill_halo_z(self, array):
            """
            fill_halo_z(self, array)
            
            
            Defined at t3dMod.F90 lines 1264-1284
            
            Parameters
            ----------
            this : T3D
            array : float array
            
            """
            _pyt3d.f90wrap_fill_halo_z(this=self._handle, array=array)
        
        def get_sz3d(self, sz3d):
            """
            get_sz3d(self, sz3d)
            
            
            Defined at t3dMod.F90 lines 1512-1517
            
            Parameters
            ----------
            this : T3D
            sz3d : int array
            
            """
            _pyt3d.f90wrap_get_sz3d(this=self._handle, sz3d=sz3d)
        
        def get_st3d(self, st3d):
            """
            get_st3d(self, st3d)
            
            
            Defined at t3dMod.F90 lines 1519-1524
            
            Parameters
            ----------
            this : T3D
            st3d : int array
            
            """
            _pyt3d.f90wrap_get_st3d(this=self._handle, st3d=st3d)
        
        def get_en3d(self, en3d):
            """
            get_en3d(self, en3d)
            
            
            Defined at t3dMod.F90 lines 1526-1533
            
            Parameters
            ----------
            this : T3D
            en3d : int array
            
            """
            _pyt3d.f90wrap_get_en3d(this=self._handle, en3d=en3d)
        
        def get_sz3dg(self, sz3dg):
            """
            get_sz3dg(self, sz3dg)
            
            
            Defined at t3dMod.F90 lines 1535-1540
            
            Parameters
            ----------
            this : T3D
            sz3dg : int array
            
            """
            _pyt3d.f90wrap_get_sz3dg(this=self._handle, sz3dg=sz3dg)
        
        def get_st3dg(self, st3dg):
            """
            get_st3dg(self, st3dg)
            
            
            Defined at t3dMod.F90 lines 1542-1547
            
            Parameters
            ----------
            this : T3D
            st3dg : int array
            
            """
            _pyt3d.f90wrap_get_st3dg(this=self._handle, st3dg=st3dg)
        
        def get_en3dg(self, en3dg):
            """
            get_en3dg(self, en3dg)
            
            
            Defined at t3dMod.F90 lines 1549-1556
            
            Parameters
            ----------
            this : T3D
            en3dg : int array
            
            """
            _pyt3d.f90wrap_get_en3dg(this=self._handle, en3dg=en3dg)
        
        def get_szx(self, szx):
            """
            get_szx(self, szx)
            
            
            Defined at t3dMod.F90 lines 1558-1563
            
            Parameters
            ----------
            this : T3D
            szx : int array
            
            """
            _pyt3d.f90wrap_get_szx(this=self._handle, szx=szx)
        
        def get_szy(self, szy):
            """
            get_szy(self, szy)
            
            
            Defined at t3dMod.F90 lines 1565-1570
            
            Parameters
            ----------
            this : T3D
            szy : int array
            
            """
            _pyt3d.f90wrap_get_szy(this=self._handle, szy=szy)
        
        def get_szz(self, szz):
            """
            get_szz(self, szz)
            
            
            Defined at t3dMod.F90 lines 1572-1577
            
            Parameters
            ----------
            this : T3D
            szz : int array
            
            """
            _pyt3d.f90wrap_get_szz(this=self._handle, szz=szz)
        
        _dt_array_initialisers = []
        
    
    f90wrap.runtime.register_class(t3d, "t3d")
    
    @staticmethod
    def optimize_decomposition(comm3d, nx, ny, nz, periodic, nghosts=None):
        """
        this = optimize_decomposition(comm3d, nx, ny, nz, periodic[, nghosts])
        
        
        Defined at t3dMod.F90 lines 1327-1410
        
        Parameters
        ----------
        comm3d : int
        nx : int
        ny : int
        nz : int
        periodic : bool array
        nghosts : int array
        
        Returns
        -------
        this : T3D
        
        """
        this = _pyt3d.f90wrap_optimize_decomposition(comm3d=comm3d, nx=nx, ny=ny, nz=nz, \
            periodic=periodic, nghosts=nghosts)
        this = f90wrap.runtime.lookup_class("t3d").from_handle(this)
        return this
    
    _dt_array_initialisers = []
    

t3dmod = T3Dmod()

