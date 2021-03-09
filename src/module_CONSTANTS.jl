module CONSTANTS

export p_WTOMJ,p_ETOM,p_CPRHO,p_GAMMA,p_CVLQ,p_CVICE,p_LF,p_LS,p_RHOWG,p_SIGMA,p_SC,p_K,p_PI,p_ThCrit

const p_WTOMJ  = 0.0864         # WTOMJ  convert Watt to MegaJoule, (MJ m-2 d-1)/(watt/m2) = 86400 s/d * .000001 MJ/J
const p_ETOM   = 0.4085         # ETOM   (mm water)/(MJ/m2) using Lv 2448 MJ/Mg and density of water 1 Mg/m3 = 1E3 mm/m / (2448 MJ/Mg * 1 Mg/m3)
const p_CPRHO  = 1240.          # CPRHO - volumetric heat capacity of air, J m-3 K-1)
const p_GAMMA  = 0.067          # GAMMA - psychrometer constant, kPa/K
const p_CVLQ   = 0.00418        # CVLQ  - volumetric heat capacity of water, MJ m-2 mm-1 K-1
const p_CVICE  = 0.00192        # CVICE - volumetric heat capacity of ice, MJ m-2 mm-1 K-1
const p_LF     = 0.335          # LF  heat of fusion of water, MJ m-2 mm-1
const p_LS     = 2.824          # LS  latent heat of sublimation of snow, MJ m-2 mm-1
const p_RHOWG  = 0.00981        # RHOWG  density of water times gravity acceleration, MPa/m or kPa/mm
const p_SIGMA  = 5.67E-08       # SIGMA  Stefan-Boltzmann constant, W m-2 K-4)
const p_SC     = 1367.          # SC  solar constant, value from Lean (1991), W/m2
const p_K      = 0.4            # K  vonKarman constant
const p_PI     = 3.1416         # PI  pi

# TODO(bernhard): In Julia we could use in-built π. However, for comparison of results we
# might still use same Pi as used by Fortran code.

# from LWFBrook90R
const p_ThCrit = 1.E-4  # minimal fraction of water content above residual water content to allow water supply for evapotranspiration
                        # p_PsiCrit is derived from p_ThCrit
const p_DTIMIN = 1.e-9  # minimum time step for iteration interval, d

# from LWFBrook90R: additionally:
# integer, dimension(12), parameter :: DAYMO = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/) ! day of the
#           real(kind=8), parameter :: DT = 1.d0           ! DT  time step for DFILE interval,  must be 1 d
#           real(kind=8), parameter :: ThCrit = 1.e-4       ! * minimal fraction of water content above residual water content to allow water supply for evapotranspiration

# from LWFBrook90R: heat flow:
#           real(kind=8), parameter :: zeroCurRange = 1.d0 ! near zero degree Celsius caused by the latent heat of fusion and thawing, see for example Boike et al. 1998
#           real(kind=8), parameter :: zeroCurTemp = 0.d0  ! zeroCurTemp temperature for the so-called zero curtain
end