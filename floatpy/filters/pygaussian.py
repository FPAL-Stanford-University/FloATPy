from __future__ import print_function, absolute_import, division
import _pygaussian
import f90wrap.runtime
import logging

class Gaussianstuff(f90wrap.runtime.FortranModule):
    """
    Module gaussianstuff
    
    
    Defined at gaussian.F90 lines 4-789
    
    """
    @f90wrap.runtime.register_class("gaussian")
    class gaussian(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=gaussian)
        
        
        Defined at gaussian.F90 lines 48-68
        
        """
        def __init__(self, n_, periodic_, handle=None):
            """
            self = Gaussian(n_, periodic_)
            
            
            Defined at gaussian.F90 lines 72-91
            
            Parameters
            ----------
            n_ : int
            periodic_ : bool
            
            Returns
            -------
            this : Gaussian
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            self._handle = _pygaussian.f90wrap_init(n_=n_, periodic_=periodic_)
        
        def __del__(self):
            """
            Destructor for class Gaussian
            
            
            Defined at gaussian.F90 lines 93-100
            
            Parameters
            ----------
            this : Gaussian
            
            """
            if self._alloc:
                _pygaussian.f90wrap_destroy(this=self._handle)
        
        def filter1(self, f, fil, nb, nc, bc1_=None, bcn_=None):
            """
            filter1(self, f, fil, nb, nc[, bc1_, bcn_])
            
            
            Defined at gaussian.F90 lines 102-334
            
            Parameters
            ----------
            this : Gaussian
            f : float array
            fil : float array
            nb : int
            nc : int
            bc1_ : int
            bcn_ : int
            
            """
            _pygaussian.f90wrap_filter1(this=self._handle, f=f, fil=fil, nb=nb, nc=nc, \
                bc1_=bc1_, bcn_=bcn_)
        
        def filter2(self, f, fil, na, nc, bc1_=None, bcn_=None):
            """
            filter2(self, f, fil, na, nc[, bc1_, bcn_])
            
            
            Defined at gaussian.F90 lines 336-564
            
            Parameters
            ----------
            this : Gaussian
            f : float array
            fil : float array
            na : int
            nc : int
            bc1_ : int
            bcn_ : int
            
            """
            _pygaussian.f90wrap_filter2(this=self._handle, f=f, fil=fil, na=na, nc=nc, \
                bc1_=bc1_, bcn_=bcn_)
        
        def filter3(self, f, fil, na, nb, bc1_=None, bcn_=None):
            """
            filter3(self, f, fil, na, nb[, bc1_, bcn_])
            
            
            Defined at gaussian.F90 lines 566-789
            
            Parameters
            ----------
            this : Gaussian
            f : float array
            fil : float array
            na : int
            nb : int
            bc1_ : int
            bcn_ : int
            
            """
            _pygaussian.f90wrap_filter3(this=self._handle, f=f, fil=fil, na=na, nb=nb, \
                bc1_=bc1_, bcn_=bcn_)
        
        _dt_array_initialisers = []
        
    
    _dt_array_initialisers = []
    

gaussianstuff = Gaussianstuff()

