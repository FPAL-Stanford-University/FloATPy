from __future__ import print_function, absolute_import, division
import _pycd10
import f90wrap.runtime
import logging

class Cd10Stuff(f90wrap.runtime.FortranModule):
    """
    Module cd10stuff
    
    
    Defined at cd10.F90 lines 4-2500
    
    """
    @f90wrap.runtime.register_class("cd10")
    class cd10(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=cd10)
        
        
        Defined at cd10.F90 lines 108-185
        
        """
        def __init__(self, n_, dx_, periodic_, bc1_, bcn_, handle=None):
            """
            self = Cd10(n_, dx_, periodic_, bc1_, bcn_)
            
            
            Defined at cd10.F90 lines 195-315
            
            Parameters
            ----------
            n_ : int
            dx_ : float
            periodic_ : bool
            bc1_ : int
            bcn_ : int
            
            Returns
            -------
            this : Cd10
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            self._handle = _pycd10.f90wrap_init(n_=n_, dx_=dx_, periodic_=periodic_, \
                bc1_=bc1_, bcn_=bcn_)
        
        def __del__(self):
            """
            Destructor for class Cd10
            
            
            Defined at cd10.F90 lines 317-349
            
            Parameters
            ----------
            this : Cd10
            
            """
            if self._alloc:
                _pycd10.f90wrap_destroy(this=self._handle)
        
        def dd1(self, f, df, na, nb, bc1_=None, bcn_=None):
            """
            dd1(self, f, df, na, nb[, bc1_, bcn_])
            
            
            Defined at cd10.F90 lines 2081-2149
            
            Parameters
            ----------
            this : Cd10
            f : float array
            df : float array
            na : int
            nb : int
            bc1_ : int
            bcn_ : int
            
            """
            _pycd10.f90wrap_dd1(this=self._handle, f=f, df=df, na=na, nb=nb, bc1_=bc1_, \
                bcn_=bcn_)
        
        def dd2(self, f, df, na, nb, bc1_=None, bcn_=None):
            """
            dd2(self, f, df, na, nb[, bc1_, bcn_])
            
            
            Defined at cd10.F90 lines 2151-2219
            
            Parameters
            ----------
            this : Cd10
            f : float array
            df : float array
            na : int
            nb : int
            bc1_ : int
            bcn_ : int
            
            """
            _pycd10.f90wrap_dd2(this=self._handle, f=f, df=df, na=na, nb=nb, bc1_=bc1_, \
                bcn_=bcn_)
        
        def dd3(self, f, df, na, nb, bc1_=None, bcn_=None):
            """
            dd3(self, f, df, na, nb[, bc1_, bcn_])
            
            
            Defined at cd10.F90 lines 2221-2289
            
            Parameters
            ----------
            this : Cd10
            f : float array
            df : float array
            na : int
            nb : int
            bc1_ : int
            bcn_ : int
            
            """
            _pycd10.f90wrap_dd3(this=self._handle, f=f, df=df, na=na, nb=nb, bc1_=bc1_, \
                bcn_=bcn_)
        
        def d2d1(self, f, df, na, nb, bc1_=None, bcn_=None):
            """
            d2d1(self, f, df, na, nb[, bc1_, bcn_])
            
            
            Defined at cd10.F90 lines 2291-2359
            
            Parameters
            ----------
            this : Cd10
            f : float array
            df : float array
            na : int
            nb : int
            bc1_ : int
            bcn_ : int
            
            """
            _pycd10.f90wrap_d2d1(this=self._handle, f=f, df=df, na=na, nb=nb, bc1_=bc1_, \
                bcn_=bcn_)
        
        def d2d2(self, f, df, na, nb, bc1_=None, bcn_=None):
            """
            d2d2(self, f, df, na, nb[, bc1_, bcn_])
            
            
            Defined at cd10.F90 lines 2361-2429
            
            Parameters
            ----------
            this : Cd10
            f : float array
            df : float array
            na : int
            nb : int
            bc1_ : int
            bcn_ : int
            
            """
            _pycd10.f90wrap_d2d2(this=self._handle, f=f, df=df, na=na, nb=nb, bc1_=bc1_, \
                bcn_=bcn_)
        
        def d2d3(self, f, df, na, nb, bc1_=None, bcn_=None):
            """
            d2d3(self, f, df, na, nb[, bc1_, bcn_])
            
            
            Defined at cd10.F90 lines 2431-2500
            
            Parameters
            ----------
            this : Cd10
            f : float array
            df : float array
            na : int
            nb : int
            bc1_ : int
            bcn_ : int
            
            """
            _pycd10.f90wrap_d2d3(this=self._handle, f=f, df=df, na=na, nb=nb, bc1_=bc1_, \
                bcn_=bcn_)
        
        _dt_array_initialisers = []
        
    
    _dt_array_initialisers = []
    

cd10stuff = Cd10Stuff()

