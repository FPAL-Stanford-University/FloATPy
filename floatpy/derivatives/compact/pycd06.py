from __future__ import print_function, absolute_import, division
import _pycd06
import f90wrap.runtime
import logging

class Cd06Stuff(f90wrap.runtime.FortranModule):
    """
    Module cd06stuff
    
    
    Defined at cd06.F90 lines 4-852
    
    """
    @f90wrap.runtime.register_class("cd06")
    class cd06(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=cd06)
        
        
        Defined at cd06.F90 lines 70-119
        
        """
        def __init__(self, n_, dx_, periodic_, bc1_, bcn_, handle=None):
            """
            self = Cd06(n_, dx_, periodic_, bc1_, bcn_)
            
            
            Defined at cd06.F90 lines 129-200
            
            Parameters
            ----------
            n_ : int
            dx_ : float
            periodic_ : bool
            bc1_ : int
            bcn_ : int
            
            Returns
            -------
            this : Cd06
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            self._handle = _pycd06.f90wrap_init(n_=n_, dx_=dx_, periodic_=periodic_, \
                bc1_=bc1_, bcn_=bcn_)
        
        def __del__(self):
            """
            Destructor for class Cd06
            
            
            Defined at cd06.F90 lines 202-218
            
            Parameters
            ----------
            this : Cd06
            
            """
            if self._alloc:
                _pycd06.f90wrap_destroy(this=self._handle)
        
        def dd1(self, f, df, na, nb):
            """
            dd1(self, f, df, na, nb)
            
            
            Defined at cd06.F90 lines 774-794
            
            Parameters
            ----------
            this : Cd06
            f : float array
            df : float array
            na : int
            nb : int
            
            """
            _pycd06.f90wrap_dd1(this=self._handle, f=f, df=df, na=na, nb=nb)
        
        def dd2(self, f, df, na, nb):
            """
            dd2(self, f, df, na, nb)
            
            
            Defined at cd06.F90 lines 796-816
            
            Parameters
            ----------
            this : Cd06
            f : float array
            df : float array
            na : int
            nb : int
            
            """
            _pycd06.f90wrap_dd2(this=self._handle, f=f, df=df, na=na, nb=nb)
        
        def dd3(self, f, df, na, nb):
            """
            dd3(self, f, df, na, nb)
            
            
            Defined at cd06.F90 lines 818-838
            
            Parameters
            ----------
            this : Cd06
            f : float array
            df : float array
            na : int
            nb : int
            
            """
            _pycd06.f90wrap_dd3(this=self._handle, f=f, df=df, na=na, nb=nb)
        
        _dt_array_initialisers = []
        
    
    _dt_array_initialisers = []
    

cd06stuff = Cd06Stuff()

