from mpi4py import MPI
import numpy
from floatpy.parallel import t3dmod
#from floatpy.derivatives.compact import CompactDerivative
#import floatpy.derivatives.compact.compact_derivative 
from floatpy.filters.filter import Filter

###############################################################
# Specifies all the features of domain, mesh, time step, etc. #
###############################################################
class NumSetting:
    def __init__(self, comm, grid_partition,
                NX=16, NY=16, NZ=20,
                XMIN=0.0, XMAX=2*numpy.pi,
                YMIN=0.0, YMAX=2*numpy.pi,
                ZMIN=-0.5,ZMAX=0.5,
                order=10):
        # Mesh configurations
        self.XMIN, self.XMAX =               XMIN, XMAX # scanning direction
        self.YMIN, self.YMAX =               YMIN, YMAX # keep this empty for 2D simulation
        self.ZMIN, self.ZMAX =               ZMIN, ZMAX # buoyancy direction
        self.NX, self.NY, self.NZ =          NX, NY, NZ
        
        self.LX = self.XMAX - self.XMIN
        self.LY = self.YMAX - self.YMIN
        self.LZ = self.ZMAX - self.ZMIN
        self.dx = self.LX / numpy.float64(self.NX)
        self.dy = self.LY / numpy.float64(self.NY)
        self.dz = self.LZ / numpy.float64(self.NZ - 1)
        
        
        # Time step
        self.dt = 2e-5
        
        # PadeOpts configurations (copied from mpitest_derivatives_compact.py)
        self.omega = 1.
        self.comm = comm#MPI.COMM_WORLD
        self.fcomm = self.comm.py2f()
        self.periodic = numpy.array([True, True, True])
        self.order = (order, order, order)
        self.grid_partition = grid_partition#t3dmod.t3d(self.fcomm, self.NX, self.NY, self.NZ, self.periodic)
        
        # Numerical settings for 3d chunk
        self.chunk_3d_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.grid_partition.get_sz3d(self.chunk_3d_size)  # size
        self.grid_partition.get_st3d(self.chunk_3d_lo)    # start index
        self.grid_partition.get_en3d(self.chunk_3d_hi)    # end index
        self.chunk_3d_lo = self.chunk_3d_lo - 1 # Convert to 0 based indexing
        self.chunk_3d_hi = self.chunk_3d_hi - 1 # Convert to 0 based indexing
        
        # Numerical settings for x transpose
        self.chunk_x_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_x_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_x_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.grid_partition.get_szx(self.chunk_x_size)    # size
        self.grid_partition.get_stx(self.chunk_x_lo)      # start index
        self.grid_partition.get_enx(self.chunk_x_hi)      # end index
        self.chunk_x_lo = self.chunk_x_lo - 1 # Convert to 0 based indexing
        self.chunk_x_hi = self.chunk_x_hi - 1 # Convert to 0 based indexing
        
        # Numerical settings for y transpose
        self.chunk_y_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_y_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_y_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.grid_partition.get_szy(self.chunk_y_size)    # size
        self.grid_partition.get_sty(self.chunk_y_lo)      # start index
        self.grid_partition.get_eny(self.chunk_y_hi)      # end index
        self.chunk_y_lo = self.chunk_y_lo - 1 # Convert to 0 based indexing
        self.chunk_y_hi = self.chunk_y_hi - 1 # Convert to 0 based indexing
        
        # Numerical settings for z transpose
        self.chunk_z_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_z_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_z_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.grid_partition.get_szz(self.chunk_z_size)
        self.grid_partition.get_stz(self.chunk_z_lo)
        self.grid_partition.get_enz(self.chunk_z_hi)
        self.chunk_z_lo = self.chunk_z_lo - 1 # Convert to 0 based indexing
        self.chunk_z_hi = self.chunk_z_hi - 1 # Convert to 0 based indexing
        
        #self.der = CompactDerivative(self.grid_partition, (self.dx, self.dy, self.dz), self.order, self.periodic)
        self.filter_type = ('compact', 'compact', 'compact')
        self.fil = Filter(self.grid_partition, self.filter_type, periodic_dimensions=self.periodic)
        
        self.x = numpy.linspace(self.XMIN, self.XMAX-self.dx, num=self.NX)[self.chunk_3d_lo[0]:self.chunk_3d_hi[0]+1]
        self.y = numpy.linspace(self.YMIN, self.YMAX-self.dy, num=self.NY)[self.chunk_3d_lo[1]:self.chunk_3d_hi[1]+1]
        self.z = numpy.linspace(self.ZMIN, self.ZMAX, num=self.NZ)[self.chunk_3d_lo[2]:self.chunk_3d_hi[2]+1]
        
        self.CNX = self.chunk_3d_size[0]
        self.CNY = self.chunk_3d_size[1]
        self.CNZ = self.chunk_3d_size[2]
        
        self.PROCESSOR_ID = MPI.COMM_WORLD.rank
        self.cout = cout
        self.cout("Settings: Domain size is {} X {} X {}.  Chunk size is {} X {} X {}".format(NX,NY,NZ,self.CNX,self.CNY,self.CNZ))




########################################
# Contains all the physical properties #
########################################
class Phys:
    def __init__(self):
        """
            Define all constants of physical properties
            """
        # Dimensionless parameters
        self.FROUDE_NUM = 1.0
        self.PECLET_NUM = 1.0
        self.PRANDTL_NUM = 1.0
        self.MARANGONI_NUM = 1.0
        self.WEBER_NUM = 1.0
        
        # Mixture properties
        self.YOUNGS_MODULUS = 1e9
        self.POISSON_RATIO = 0.6
        self.DENSITY_SLD_REF = 1.0
        self.DENSITY_GAS_REF = 2.23e-4     # rho_gas / rho_solid
        self.DENSITY_LIQ_REF = 1.0
        self.CONDUCTIVITY_SLD = 1.0
        self.CONDUCTIVITY_GAS = 0.016/30.0 # compare to SS
        self.CONDUCTIVITY_LIQ = 1.0
        self.CP_SLD = 1.0
        self.CP_LIQ = 1.0
        self.CP_GAS = 523.0/500.0
        self.TEMP_REF = 1.0                # melting temperature
        self.LATENT_HEAT = 0.1
        self.ENTHALPY_GAS_REF = 1.0                                       # reference enthalpy of gas at metal melting temp
        self.ENTHALPY_SLD_REF = 1.0                                       # reference enthalpy of solid at metal melting temp
        self.ENTHALPY_LIQ_REF = self.ENTHALPY_SLD_REF + self.LATENT_HEAT  # reference enthalpy of liquid at metal melting temp
        
        # Laser properties
        self.LASER_POWER = 1e0 * self.LATENT_HEAT
        self.LASER_FWHM = 0.10             # Laser full-width-half-maximum diameter
        self.EXTINCTION_COEF_GAS = 0.0
        self.EXTINCTION_COEF_SLD = 30.0
        self.EXTINCTION_COEF_LIQ = 30.0
        
        # Scanning setting
        self.laser_loc_x = 0.0
        self.laser_loc_y = 0.0
        self.scan_vel = 0.0
        
        cout("Physical properties have been initialized.")
    
    
    def getLameCoefs(self):
        """
            Return the Lame coefficients for linear elasticity of the solid
            """
        E = self.YOUNGS_MODULUS
        v = self.POISSON_RATIO
        lam = E * v / ((1 + v) * (1 - 2*v))
        mu = 0.5 * E / (1 + v)
        return lam, mu
    
    
    def getSolidDensity(self, T):
        """
            Return the temperature dependent temperature of solid.
            """
        #rs_T = self.DENSITY_SLD_REF + 1e-7 * numpy.tanh(0.1 * T)
        rs_T = self.DENSITY_SLD_REF
        return rs_T
    
    
    def getGasDensity(self,T):
        """
            Return the temperature dependent temperature of gas.
            """
        return self.DENSITY_GAS_REF
    
    
    def getLiquidDensity(self,T):
        """
            Return the temperature dependent temperature of liquid.
            """
        return self.DENSITY_LIQ_REF
    
    
    def getEnthalpy(self,sVF,gVF,r,T):
        """
            Calculate overall enthalpy from temperature and volume fractions.
            args:
            sVF - solid volumne fraction
            gVF - gas volume fraction
            r   - overall density of the mixture
            T   - temperature
            """
        rg = self.getGasDensity(T)
        rs = self.getSolidDensity(T)
        
        sMF = rs * sVF / r      # mass fraction of solid
        gMF = rg * gVF / r      # mass fraction of gas
        lMF = 1.0 - sMF - gMF   # mass fraction of liquid
        
        hg = self.ENTHALPY_GAS_REF + self.CP_GAS * (T - self.TEMP_REF) # specific enthalpy of gas
        hs = self.ENTHALPY_SLD_REF + self.CP_SLD * (T - self.TEMP_REF) # specific enthalpy of solid
        hl = self.ENTHALPY_LIQ_REF + self.CP_LIQ * (T - self.TEMP_REF) # specific enthalpy of liquid
        
        return sMF * hs + lMF * hl + gMF * hg  # specific enthalpy of the mixture
    
    
    def getTemperature(self,sMF,lMF,gMF,h):
        """
            Calculate the mixture temperature at equilibrium state from the enthalpy and volume fractions.
            args:
            sMF - mass fraction of solid
            lMF - mass fraction of liquid
            gMF - mass fraction of gas
            h   - overall enthalpy of the mixture
            """
        hs = self.ENTHALPY_SLD_REF
        hl = self.ENTHALPY_LIQ_REF
        hg = self.ENTHALPY_GAS_REF
        h0 = sMF * hs + lMF * hl + gMF * hg                            # mixture reference enthalpy
        cp = sMF * self.CP_SLD + lMF * self.CP_LIQ + gMF * self.CP_GAS # mixture specific heat
        T = (h - h0) / cp + self.TEMP_REF
        
        # Correct the temperature for mixture at phase non-equilibrium.
        hHigh = hl * (lMF + sMF) + hg * gMF     # high limit of two-phase enthalpy region
        hLow  = hs * (lMF + sMF) + hg * gMF     # low  limit of two-phase enthalpy region
        ind = numpy.where(numpy.logical_and(h<=hHigh,h>=hLow))
        #T[ind] = self.TEMP_REF
        
        return T, hHigh, hLow
    
    
    def getThermalConductivity(self, T, gVF, sVF, lVF):
        """
            Return the relative thermal conductivity at given temperature.
            """
        kg = self.CONDUCTIVITY_GAS
        kl = self.CONDUCTIVITY_LIQ
        ks = self.CONDUCTIVITY_SLD
        k = kg * gVF + ks * sVF + kl * lVF
        return k
    
    
    def getSpecificHeat(self, gMF, sMF, lMF):
        """
            Return specific heat for mixture.
            """
        cps, cpl, cpg = self.CP_SLD, self.CP_LIQ, self.CP_GAS
        cp = cps * sMF + cpl * lMF + cpg * gMF
        return cp
    
    
    
    def getViscosity(self, T, lVF, gVF):
        """
            Return the kinematic viscosity.
            """
        visc = lVF * 1.0 + gVF * 1e-2 + (1.0-lVF-gVF) * 0.
        return visc
    
    
    def getGasPhaseMobility(self, T):
        return 1.0
    
    
    def getMetalPhaseChangeMobility(self, T):
        return 1.0
    
    
    def getMetalInterfaceMobility(self, T):
        return 1.0

def cout(msg,indent=0,topspace=0,bottomspace=0):
    PROCESSOR_ID = MPI.COMM_WORLD.rank
    print ("\n"*topspace + "[PROCESSOR {}] " + " "*indent*4 + ">>> " + msg + "\n"*bottomspace).format(PROCESSOR_ID)
    return PROCESSOR_ID
