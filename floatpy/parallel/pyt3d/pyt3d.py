from __future__ import print_function, absolute_import, division
import _pyt3d
import f90wrap.runtime
import logging

class T3Dmod(f90wrap.runtime.FortranModule):
    """
    Module t3dmod
    
    
    Defined at t3dMod.F90 lines 1-1644
    
    """
    class t3d(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t3d)
        
        
        Defined at t3dMod.F90 lines 19-90
        
        """
        def init(self, comm3d, nx, ny, nz, px, py, pz, periodic_, reorder, fail, \
            nghosts=None, createcrosscommunicators=None):
            """
            init(self, comm3d, nx, ny, nz, px, py, pz, periodic_, reorder, fail[, nghosts, \
                createcrosscommunicators])
            
            
            Defined at t3dMod.F90 lines 98-506
            
            Parameters
            ----------
            this : T3D
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
            
            """
            _pyt3d.f90wrap_init(this=self._handle, comm3d=comm3d, nx=nx, ny=ny, nz=nz, \
                px=px, py=py, pz=pz, periodic_=periodic_, reorder=reorder, fail=fail, \
                nghosts=nghosts, createcrosscommunicators=createcrosscommunicators)
        
        def __del__(self):
            """
            Destructor for class T3D
            
            
            Defined at t3dMod.F90 lines 508-524
            
            Parameters
            ----------
            this : T3D
            
            """
            if self._alloc:
                _pyt3d.f90wrap_destroy(this=self._handle)
        
        def transpose_3d_to_x(self, input, output):
            """
            transpose_3d_to_x(self, input, output)
            
            
            Defined at t3dMod.F90 lines 526-578
            
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
            
            
            Defined at t3dMod.F90 lines 580-621
            
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
            
            
            Defined at t3dMod.F90 lines 623-675
            
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
            
            
            Defined at t3dMod.F90 lines 677-719
            
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
            
            
            Defined at t3dMod.F90 lines 721-773
            
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
            
            
            Defined at t3dMod.F90 lines 775-817
            
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
            
            
            Defined at t3dMod.F90 lines 1165-1185
            
            Parameters
            ----------
            this : T3D
            array : float array
            
            """
            _pyt3d.f90wrap_fill_halo_x(this=self._handle, array=array)
        
        def fill_halo_y(self, array):
            """
            fill_halo_y(self, array)
            
            
            Defined at t3dMod.F90 lines 1187-1207
            
            Parameters
            ----------
            this : T3D
            array : float array
            
            """
            _pyt3d.f90wrap_fill_halo_y(this=self._handle, array=array)
        
        def fill_halo_z(self, array):
            """
            fill_halo_z(self, array)
            
            
            Defined at t3dMod.F90 lines 1209-1230
            
            Parameters
            ----------
            this : T3D
            array : float array
            
            """
            _pyt3d.f90wrap_fill_halo_z(this=self._handle, array=array)
        
        def __init__(self, comm3d, nx, ny, nz, periodic, nghosts=None, handle=None):
            """
            self = T3D(comm3d, nx, ny, nz, periodic[, nghosts])
            
            
            Defined at t3dMod.F90 lines 1273-1356
            
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
            f90wrap.runtime.FortranDerivedType.__init__(self)
            self._handle = _pyt3d.f90wrap_optimize_decomposition(comm3d=comm3d, nx=nx, \
                ny=ny, nz=nz, periodic=periodic, nghosts=nghosts)
        
        def get_sz3d(self, sz3d):
            """
            get_sz3d(self, sz3d)
            
            
            Defined at t3dMod.F90 lines 1458-1463
            
            Parameters
            ----------
            this : T3D
            sz3d : int array
            
            """
            _pyt3d.f90wrap_get_sz3d(this=self._handle, sz3d=sz3d)
        
        def get_st3d(self, st3d):
            """
            get_st3d(self, st3d)
            
            
            Defined at t3dMod.F90 lines 1465-1470
            
            Parameters
            ----------
            this : T3D
            st3d : int array
            
            """
            _pyt3d.f90wrap_get_st3d(this=self._handle, st3d=st3d)
        
        def get_en3d(self, en3d):
            """
            get_en3d(self, en3d)
            
            
            Defined at t3dMod.F90 lines 1472-1479
            
            Parameters
            ----------
            this : T3D
            en3d : int array
            
            """
            _pyt3d.f90wrap_get_en3d(this=self._handle, en3d=en3d)
        
        def get_sz3dg(self, sz3dg):
            """
            get_sz3dg(self, sz3dg)
            
            
            Defined at t3dMod.F90 lines 1481-1486
            
            Parameters
            ----------
            this : T3D
            sz3dg : int array
            
            """
            _pyt3d.f90wrap_get_sz3dg(this=self._handle, sz3dg=sz3dg)
        
        def get_st3dg(self, st3dg):
            """
            get_st3dg(self, st3dg)
            
            
            Defined at t3dMod.F90 lines 1488-1493
            
            Parameters
            ----------
            this : T3D
            st3dg : int array
            
            """
            _pyt3d.f90wrap_get_st3dg(this=self._handle, st3dg=st3dg)
        
        def get_en3dg(self, en3dg):
            """
            get_en3dg(self, en3dg)
            
            
            Defined at t3dMod.F90 lines 1495-1502
            
            Parameters
            ----------
            this : T3D
            en3dg : int array
            
            """
            _pyt3d.f90wrap_get_en3dg(this=self._handle, en3dg=en3dg)
        
        def get_szx(self, szx):
            """
            get_szx(self, szx)
            
            
            Defined at t3dMod.F90 lines 1504-1509
            
            Parameters
            ----------
            this : T3D
            szx : int array
            
            """
            _pyt3d.f90wrap_get_szx(this=self._handle, szx=szx)
        
        def get_stx(self, stx):
            """
            get_stx(self, stx)
            
            
            Defined at t3dMod.F90 lines 1511-1516
            
            Parameters
            ----------
            this : T3D
            stx : int array
            
            """
            _pyt3d.f90wrap_get_stx(this=self._handle, stx=stx)
        
        def get_enx(self, enx):
            """
            get_enx(self, enx)
            
            
            Defined at t3dMod.F90 lines 1518-1524
            
            Parameters
            ----------
            this : T3D
            enx : int array
            
            """
            _pyt3d.f90wrap_get_enx(this=self._handle, enx=enx)
        
        def get_szy(self, szy):
            """
            get_szy(self, szy)
            
            
            Defined at t3dMod.F90 lines 1526-1531
            
            Parameters
            ----------
            this : T3D
            szy : int array
            
            """
            _pyt3d.f90wrap_get_szy(this=self._handle, szy=szy)
        
        def get_sty(self, sty):
            """
            get_sty(self, sty)
            
            
            Defined at t3dMod.F90 lines 1533-1538
            
            Parameters
            ----------
            this : T3D
            sty : int array
            
            """
            _pyt3d.f90wrap_get_sty(this=self._handle, sty=sty)
        
        def get_eny(self, eny):
            """
            get_eny(self, eny)
            
            
            Defined at t3dMod.F90 lines 1540-1546
            
            Parameters
            ----------
            this : T3D
            eny : int array
            
            """
            _pyt3d.f90wrap_get_eny(this=self._handle, eny=eny)
        
        def get_szz(self, szz):
            """
            get_szz(self, szz)
            
            
            Defined at t3dMod.F90 lines 1548-1553
            
            Parameters
            ----------
            this : T3D
            szz : int array
            
            """
            _pyt3d.f90wrap_get_szz(this=self._handle, szz=szz)
        
        def get_stz(self, stz):
            """
            get_stz(self, stz)
            
            
            Defined at t3dMod.F90 lines 1555-1560
            
            Parameters
            ----------
            this : T3D
            stz : int array
            
            """
            _pyt3d.f90wrap_get_stz(this=self._handle, stz=stz)
        
        def get_enz(self, enz):
            """
            get_enz(self, enz)
            
            
            Defined at t3dMod.F90 lines 1562-1567
            
            Parameters
            ----------
            this : T3D
            enz : int array
            
            """
            _pyt3d.f90wrap_get_enz(this=self._handle, enz=enz)
        
        def comm3d(self):
            """
            comm3d = comm3d(self)
            
            
            Defined at t3dMod.F90 lines 1569-1574
            
            Parameters
            ----------
            this : T3D
            
            Returns
            -------
            comm3d : int
            
            """
            comm3d = _pyt3d.f90wrap_comm3d(this=self._handle)
            return comm3d
        
        def commx(self):
            """
            commx = commx(self)
            
            
            Defined at t3dMod.F90 lines 1576-1581
            
            Parameters
            ----------
            this : T3D
            
            Returns
            -------
            commx : int
            
            """
            commx = _pyt3d.f90wrap_commx(this=self._handle)
            return commx
        
        def commy(self):
            """
            commy = commy(self)
            
            
            Defined at t3dMod.F90 lines 1583-1588
            
            Parameters
            ----------
            this : T3D
            
            Returns
            -------
            commy : int
            
            """
            commy = _pyt3d.f90wrap_commy(this=self._handle)
            return commy
        
        def commz(self):
            """
            commz = commz(self)
            
            
            Defined at t3dMod.F90 lines 1590-1595
            
            Parameters
            ----------
            this : T3D
            
            Returns
            -------
            commz : int
            
            """
            commz = _pyt3d.f90wrap_commz(this=self._handle)
            return commz
        
        def commxy(self):
            """
            commxy = commxy(self)
            
            
            Defined at t3dMod.F90 lines 1597-1602
            
            Parameters
            ----------
            this : T3D
            
            Returns
            -------
            commxy : int
            
            """
            commxy = _pyt3d.f90wrap_commxy(this=self._handle)
            return commxy
        
        def commyz(self):
            """
            commyz = commyz(self)
            
            
            Defined at t3dMod.F90 lines 1604-1609
            
            Parameters
            ----------
            this : T3D
            
            Returns
            -------
            commyz : int
            
            """
            commyz = _pyt3d.f90wrap_commyz(this=self._handle)
            return commyz
        
        def commxz(self):
            """
            commxz = commxz(self)
            
            
            Defined at t3dMod.F90 lines 1611-1616
            
            Parameters
            ----------
            this : T3D
            
            Returns
            -------
            commxz : int
            
            """
            commxz = _pyt3d.f90wrap_commxz(this=self._handle)
            return commxz
        
        def px(self):
            """
            px = px(self)
            
            
            Defined at t3dMod.F90 lines 1618-1623
            
            Parameters
            ----------
            this : T3D
            
            Returns
            -------
            px : int
            
            """
            px = _pyt3d.f90wrap_px(this=self._handle)
            return px
        
        def py(self):
            """
            py = py(self)
            
            
            Defined at t3dMod.F90 lines 1625-1630
            
            Parameters
            ----------
            this : T3D
            
            Returns
            -------
            py : int
            
            """
            py = _pyt3d.f90wrap_py(this=self._handle)
            return py
        
        def pz(self):
            """
            pz = pz(self)
            
            
            Defined at t3dMod.F90 lines 1632-1637
            
            Parameters
            ----------
            this : T3D
            
            Returns
            -------
            pz : int
            
            """
            pz = _pyt3d.f90wrap_pz(this=self._handle)
            return pz
        
        def nprocs(self):
            """
            nprocs = nprocs(self)
            
            
            Defined at t3dMod.F90 lines 1639-1644
            
            Parameters
            ----------
            this : T3D
            
            Returns
            -------
            nprocs : int
            
            """
            nprocs = _pyt3d.f90wrap_nprocs(this=self._handle)
            return nprocs
        
        _dt_array_initialisers = []
        
    
    f90wrap.runtime.register_class(t3d, "t3d")
    
    _dt_array_initialisers = []
    

t3dmod = T3Dmod()

