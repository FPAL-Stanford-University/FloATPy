from __future__ import print_function, absolute_import, division
import _compact_sixth_order
import f90wrap.runtime
import logging

class Kind_Parameters(f90wrap.runtime.FortranModule):
    """
    Module kind_parameters
    
    
    Defined at kind_parameters.F90 lines 3-16
    
    """
    @property
    def rkind(self):
        """
        Element rkind ftype=integer pytype=int
        
        
        Defined at kind_parameters.F90 line 10
        
        """
        return _compact_sixth_order.f90wrap_kind_parameters__get__rkind()
    
    @property
    def clen(self):
        """
        Element clen ftype=integer pytype=int
        
        
        Defined at kind_parameters.F90 line 12
        
        """
        return _compact_sixth_order.f90wrap_kind_parameters__get__clen()
    
    @property
    def stdin(self):
        """
        Element stdin ftype=integer pytype=int
        
        
        Defined at kind_parameters.F90 line 14
        
        """
        return _compact_sixth_order.f90wrap_kind_parameters__get__stdin()
    
    @property
    def stdout(self):
        """
        Element stdout ftype=integer pytype=int
        
        
        Defined at kind_parameters.F90 line 15
        
        """
        return _compact_sixth_order.f90wrap_kind_parameters__get__stdout()
    
    @property
    def stderr(self):
        """
        Element stderr ftype=integer pytype=int
        
        
        Defined at kind_parameters.F90 line 16
        
        """
        return _compact_sixth_order.f90wrap_kind_parameters__get__stderr()
    
    def __str__(self):
        ret = ['<kind_parameters>{\n']
        ret.append('    rkind : ')
        ret.append(repr(self.rkind))
        ret.append(',\n    clen : ')
        ret.append(repr(self.clen))
        ret.append(',\n    stdin : ')
        ret.append(repr(self.stdin))
        ret.append(',\n    stdout : ')
        ret.append(repr(self.stdout))
        ret.append(',\n    stderr : ')
        ret.append(repr(self.stderr))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

kind_parameters = Kind_Parameters()

class Constants(f90wrap.runtime.FortranModule):
    """
    Module constants
    
    
    Defined at constants.F90 lines 1-31
    
    """
    @property
    def zero(self):
        """
        Element zero ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 6
        
        """
        return _compact_sixth_order.f90wrap_constants__get__zero()
    
    @property
    def one(self):
        """
        Element one ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 7
        
        """
        return _compact_sixth_order.f90wrap_constants__get__one()
    
    @property
    def two(self):
        """
        Element two ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 8
        
        """
        return _compact_sixth_order.f90wrap_constants__get__two()
    
    @property
    def three(self):
        """
        Element three ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 9
        
        """
        return _compact_sixth_order.f90wrap_constants__get__three()
    
    @property
    def four(self):
        """
        Element four ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 10
        
        """
        return _compact_sixth_order.f90wrap_constants__get__four()
    
    @property
    def five(self):
        """
        Element five ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 11
        
        """
        return _compact_sixth_order.f90wrap_constants__get__five()
    
    @property
    def six(self):
        """
        Element six ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 12
        
        """
        return _compact_sixth_order.f90wrap_constants__get__six()
    
    @property
    def seven(self):
        """
        Element seven ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 13
        
        """
        return _compact_sixth_order.f90wrap_constants__get__seven()
    
    @property
    def eight(self):
        """
        Element eight ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 14
        
        """
        return _compact_sixth_order.f90wrap_constants__get__eight()
    
    @property
    def nine(self):
        """
        Element nine ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 15
        
        """
        return _compact_sixth_order.f90wrap_constants__get__nine()
    
    @property
    def ten(self):
        """
        Element ten ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 16
        
        """
        return _compact_sixth_order.f90wrap_constants__get__ten()
    
    @property
    def half(self):
        """
        Element half ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 18
        
        """
        return _compact_sixth_order.f90wrap_constants__get__half()
    
    @property
    def third(self):
        """
        Element third ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 19
        
        """
        return _compact_sixth_order.f90wrap_constants__get__third()
    
    @property
    def fourth(self):
        """
        Element fourth ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 20
        
        """
        return _compact_sixth_order.f90wrap_constants__get__fourth()
    
    @property
    def fifth(self):
        """
        Element fifth ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 21
        
        """
        return _compact_sixth_order.f90wrap_constants__get__fifth()
    
    @property
    def sixth(self):
        """
        Element sixth ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 22
        
        """
        return _compact_sixth_order.f90wrap_constants__get__sixth()
    
    @property
    def seventh(self):
        """
        Element seventh ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 23
        
        """
        return _compact_sixth_order.f90wrap_constants__get__seventh()
    
    @property
    def eighth(self):
        """
        Element eighth ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 24
        
        """
        return _compact_sixth_order.f90wrap_constants__get__eighth()
    
    @property
    def twothird(self):
        """
        Element twothird ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 26
        
        """
        return _compact_sixth_order.f90wrap_constants__get__twothird()
    
    @property
    def fourthird(self):
        """
        Element fourthird ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 27
        
        """
        return _compact_sixth_order.f90wrap_constants__get__fourthird()
    
    @property
    def pi(self):
        """
        Element pi ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 29
        
        """
        return _compact_sixth_order.f90wrap_constants__get__pi()
    
    @property
    def kappa(self):
        """
        Element kappa ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 31
        
        """
        return _compact_sixth_order.f90wrap_constants__get__kappa()
    
    @property
    def eps(self):
        """
        Element eps ftype=real(rkind) pytype=float
        
        
        Defined at constants.F90 line 32
        
        """
        return _compact_sixth_order.f90wrap_constants__get__eps()
    
    def __str__(self):
        ret = ['<constants>{\n']
        ret.append('    zero : ')
        ret.append(repr(self.zero))
        ret.append(',\n    one : ')
        ret.append(repr(self.one))
        ret.append(',\n    two : ')
        ret.append(repr(self.two))
        ret.append(',\n    three : ')
        ret.append(repr(self.three))
        ret.append(',\n    four : ')
        ret.append(repr(self.four))
        ret.append(',\n    five : ')
        ret.append(repr(self.five))
        ret.append(',\n    six : ')
        ret.append(repr(self.six))
        ret.append(',\n    seven : ')
        ret.append(repr(self.seven))
        ret.append(',\n    eight : ')
        ret.append(repr(self.eight))
        ret.append(',\n    nine : ')
        ret.append(repr(self.nine))
        ret.append(',\n    ten : ')
        ret.append(repr(self.ten))
        ret.append(',\n    half : ')
        ret.append(repr(self.half))
        ret.append(',\n    third : ')
        ret.append(repr(self.third))
        ret.append(',\n    fourth : ')
        ret.append(repr(self.fourth))
        ret.append(',\n    fifth : ')
        ret.append(repr(self.fifth))
        ret.append(',\n    sixth : ')
        ret.append(repr(self.sixth))
        ret.append(',\n    seventh : ')
        ret.append(repr(self.seventh))
        ret.append(',\n    eighth : ')
        ret.append(repr(self.eighth))
        ret.append(',\n    twothird : ')
        ret.append(repr(self.twothird))
        ret.append(',\n    fourthird : ')
        ret.append(repr(self.fourthird))
        ret.append(',\n    pi : ')
        ret.append(repr(self.pi))
        ret.append(',\n    kappa : ')
        ret.append(repr(self.kappa))
        ret.append(',\n    eps : ')
        ret.append(repr(self.eps))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

constants = Constants()

class Cd06Stuff(f90wrap.runtime.FortranModule):
    """
    Module cd06stuff
    
    
    Defined at cd06.F90 lines 4-852
    
    """
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
            self._handle = _compact_sixth_order.f90wrap_init(n_=n_, dx_=dx_, periodic_=periodic_, \
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
                _compact_sixth_order.f90wrap_destroy(this=self._handle)
        
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
            _compact_sixth_order.f90wrap_dd1(this=self._handle, f=f, df=df, na=na, nb=nb)
        
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
            _compact_sixth_order.f90wrap_dd2(this=self._handle, f=f, df=df, na=na, nb=nb)
        
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
            _compact_sixth_order.f90wrap_dd3(this=self._handle, f=f, df=df, na=na, nb=nb)
        
        _dt_array_initialisers = []
        
    
    f90wrap.runtime.register_class(cd06, "cd06")
    
    _dt_array_initialisers = []
    

cd06stuff = Cd06Stuff()

