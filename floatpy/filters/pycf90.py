from __future__ import print_function, absolute_import, division
import _pycf90
import f90wrap.runtime
import logging

class Cf90Stuff(f90wrap.runtime.FortranModule):
    """
    Module cf90stuff
    
    
    Defined at cf90.F90 lines 4-1390
    
    """
    @f90wrap.runtime.register_class("cf90")
    class cf90(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=cf90)
        
        
        Defined at cf90.F90 lines 53-98
        
        """
        def __init__(self, n_, periodic_, handle=None):
            """
            self = Cf90(n_, periodic_)
            
            
            Defined at cf90.F90 lines 102-170
            
            Parameters
            ----------
            n_ : int
            periodic_ : bool
            
            Returns
            -------
            this : Cf90
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            self._handle = _pycf90.f90wrap_init(n_=n_, periodic_=periodic_)
        
        def __del__(self):
            """
            Destructor for class Cf90
            
            
            Defined at cf90.F90 lines 172-191
            
            Parameters
            ----------
            this : Cf90
            
            """
            if self._alloc:
                _pycf90.f90wrap_destroy(this=self._handle)
        
        def filter1(self, f, df, na, nb, bc1_=None, bcn_=None):
            """
            filter1(self, f, df, na, nb[, bc1_, bcn_])
            
            
            Defined at cf90.F90 lines 1176-1246
            
            Parameters
            ----------
            this : Cf90
            f : float array
            df : float array
            na : int
            nb : int
            bc1_ : int
            bcn_ : int
            
            """
            _pycf90.f90wrap_filter1(this=self._handle, f=f, df=df, na=na, nb=nb, bc1_=bc1_, \
                bcn_=bcn_)
        
        def filter2(self, f, df, na, nb, bc1_=None, bcn_=None):
            """
            filter2(self, f, df, na, nb[, bc1_, bcn_])
            
            
            Defined at cf90.F90 lines 1248-1318
            
            Parameters
            ----------
            this : Cf90
            f : float array
            df : float array
            na : int
            nb : int
            bc1_ : int
            bcn_ : int
            
            """
            _pycf90.f90wrap_filter2(this=self._handle, f=f, df=df, na=na, nb=nb, bc1_=bc1_, \
                bcn_=bcn_)
        
        def filter3(self, f, df, na, nb, bc1_=None, bcn_=None):
            """
            filter3(self, f, df, na, nb[, bc1_, bcn_])
            
            
            Defined at cf90.F90 lines 1320-1390
            
            Parameters
            ----------
            this : Cf90
            f : float array
            df : float array
            na : int
            nb : int
            bc1_ : int
            bcn_ : int
            
            """
            _pycf90.f90wrap_filter3(this=self._handle, f=f, df=df, na=na, nb=nb, bc1_=bc1_, \
                bcn_=bcn_)
        
        _dt_array_initialisers = []
        
    
    _dt_array_initialisers = []
    

cf90stuff = Cf90Stuff()

