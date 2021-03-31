"""
# Interception and transpiration
Text copied from Ecoshift on module EVP:

"
BROOK90 routines in this module relate to interception and actual transpiration.
Subroutines INTER and INTER24, which handle the interception of rain or snow by the plant
canopy are equivalent routines. INTER24 is used when parameter `p_NPINT` is 1 and
precipitation is input once a day; it assumes that the daily precipitation all occurs
within DURATN hours in the middle of the day. INTER is used when `p_NPINT > 1` and
precipitation is input more than once a day; it assumes that precipitation rate is
constant through the precipitation time step. INTER and INTER24 are used both for rain and
snow, with different calling parameters and variables. PLNTRES calculates parameters
related to rhizosphere, root, and xylem resistance; it is called once at the beginning of
each day. These parameters affect only soil water supply rate, not potential
transpiration. TBYLAYER calculates the daily transpiration from each soil layer from the
potential transpiration (PTRAN) and the total soil water potential in each layer (PSITI).
"
"""
module EVP # Interception and transpiration

using ..CONSTANTS: p_RHOWG, p_PI # https://discourse.julialang.org/t/large-programs-structuring-modules-include-such-that-to-increase-performance-and-readability/29102/5

export PLNTRES, TBYLAYER, INTER, INTER24

"""
    PLNTRES()

Allocates total plant resistance to xylem and root layers.

Ecoshift:
Subroutine PLNTRES is called at the beginning of each day to obtain resistivities to liquid
water flow: rhizosphere resistivity for each soil layer, root resistivity in each soil
layer, and xylem resistivity. These parameters, together with soil water potential in each
layer (PSITI) and critical plant water potential (PSICR) control the supply of water to
transpiring leaves and thus the reduction of actual transpiration below potential
transpiration. As defined by Hunt et al. (1991) the resistances used here are "potential
difference resistivities", because the transpiration flux rate is in units of mm/d and the
potential gradient is in MPa. The resistivities have units of MPa d mm-1.

In most of this subsection and the following subsection (TBYLAYER), algebraic notation is
used instead of variable names. The correspondence is:

rri RROOTI()    ψt  PSIT    P   PTR fi  RTFRAC rx  RXYLEM  ψti PSITI() T   ATR Di  D() rp
RPLANT  ψc  PSICR   Ti  ATRANI()    Li  RTDENI ri  RI          S   SUPPLY  δi  DELT rt  RT
αi  ALPHA() The following additional parameters or constants are input to the routines

fx  FXYLEM  di  RELDEN  ρwg RHOWG R1  RTRAD   d   DISPC   π   PI Lr  RTLEN   Ki  KK()

Several other algebraic variables occur in the derivations below, but are not needed in the
program.

If soil water potential is uniform through the root zone and rhizosphere resistivity is
negligible, a bulk plant resistivity, rp, can be defined by

rp = ( ψt - ψ - ρw gd ) / T'

where T' is the transpiration rate, ψt is the soil water potential, and ψ is the leaf water
potential. The ρwgd term is the gravity potential difference above the ground surface, where
ρwg is the density of water times the acceleration of gravity, and d is the effective height
of canopy evaporation, taken here as the closed canopy zero-plane displacement (DISPC).
Change in water storage within the plant is ignored in BROOK90, as Hunt et al. (1991) have
shown that it does not matter to total daily transpiration.

In BROOK90, rp is determined primarily by the maximum bulk plant conductivity (MXKPL), which
is an input parameter. MXKPL is the water uptake rate for a closed canopy per unit of soil
to leaf potential difference. When the soil is wet so ψt is close to 0, many plants have a ψ
of around -1.5 MPa when transpiration rate is about 0.5 mm/hr, a normal midday rate for a
sunny day in temperate regions. This gives a MXKPL of 8 mm d-1 MPa-1. Abdul-Jabbar et al.
(1984) found literature values to range only from 7 to 30 mm d-1 MPa-1, which seems a
surprisingly narrow range for plant canopies of any species. Actual plant resistivity, rp,
is determined in subroutine CANOPY.

Fig. EVP-1. Resistances and potentials in the liquid flow path for transpiration, for 3
layers with roots. ψ is the leaf potential, ψx is potential at gorund level, and ψti are
total soil water potentials. rₓ is xylem resistance, rᵣᵢ are root resistances, and rₛᵢ are
rhizosphere resistances.

Figure EVP-1 shows the resistance network used in BROOK90. Each soil layer is considered to
have a rhizosphere resistivity, rsi, and a root resistivity, rri, in series. Each layer is
considered in parallel with the others. An additional resistivity in series accounts for
resistance to flow through the xylem above ground level, rx. The total soil-water potential,
ψti, differs among layers. The leaf water potential, ψ, is assumed constant through the
canopy. The potential at the ground surface is ψx. This system and the parameterization of
it was developed by Federer (1979) and Hunt and others(1991).

The xylem resistance, rx, is

rx = fx rp ;

where fx (FXYLEM) is an input parameter, which is probably close to zero for short canopies,
but may be 0.5 or higher for forests (Hunt et al. 1991). The use of fx here specifies the
amount of rp (RPLANT) that should not be divided among soil layers, so it is effectively the
fraction of the plant resistance that is above ground level.

The plant resistance that is below ground level, rp - rx, is the parallel combination of the
individual layer resistances, rri (Fig. EVP-1). BROOK90 assumes that the root resistivity
per unit length of root is constant, so rri (RROOTI) is inversely proportional to the
fraction of total root length that is in each layer, fi

rri = ( rp - rx ) / fi .

BROOK90 parameterizes root distribution by the relative density of roots in each soil layer
di (RELDENi). RELDENi is obtained from ROOTDEN in subroutine RTDEN . Then fi is

fi = di Di / S ( di Di )

where Di is the stone-free layer thickness, THICK(i) * (1 - STONEF(i)).

When the soil is at uniform water potential and rhizosphere resistances, rsi , are
negligible and fx = 0, then the relative transpiration withdrawal from each layer is
proportional to fi (see subroutine TBYLAYER). Increasing fx makes the withdrawal less
dependent on fi .

Usually the rhizosphere resistance only becomes significant when the soil is dry. Following
Federer (1979) and Cowan (1965), the rhizosphere resistance, rsi , calculated in subroutine
TBYLAYER is rsi = αi / Ki , where Ki is the hydraulic conductivity of the rhizosphere. The
value of αi (ALPHAi) for a layer is

αi = Ai / ρwg Di

where Ai from Cowan (1965) is

Ai = ( 1 / 8 π Li ) [ δi - 3 - 2 ( ln δi ) / ( 1 - δi ) ]

where Li is the root density in the layer (mm/mm3) and δi is the root volume fraction in the
layer, obtained as

δi = π R12 Li

where R1 is the average radius of the absorbing roots (RTRAD) , which is an input parameter.
(Note that δi has no relation to the daylength, δ ). The root density is

Li = fi Lr / Di

where Lr is the total length of absorbing roots per unit ground area in mm/mm2. The value of
Lr in m/m2 is found in subroutine CANOPY as the parameter MXRTLN reduced by DENSEF and
RELHT. The dependence on the seasonal RELHT assumes that root length increases
proportionally with height growth. MXRTLN is an input parameter expressing the length of
absorbing roots per unit ground area in m/m2. Lr and R1 only affect rhizosphere resistance
and thus are only important for dry soil or when Lr is small.
"
"""
function PLNTRES(NLAYER, p_soil, p_RTLEN, p_fT_RELDEN, p_RTRAD, p_RPLANT, p_FXYLEM, p_PI, p_RHOWG) # TODO(bernhard) find out how to import p_PI and p_RHOWG from the moduel CONSTANTS
    p_THICK = p_soil.p_THICK
    p_STONEF = p_soil.p_STONEF
    # compute stone-free layer thickness D
    D = fill(NaN, NLAYER)
    for i = 1:NLAYER
        D[i] = p_THICK[i] * (1 - p_STONEF[i])
    end

    # compute fraction of total root length that is in each layer
    RTFRAC = p_fT_RELDEN[1:NLAYER].*D[1:NLAYER] ./ sum(p_fT_RELDEN[1:NLAYER].*D[1:NLAYER])

    # compute xylem resistance, MPa d/mm
    RXYLEM = p_FXYLEM * p_RPLANT
    # compute RROOTI and ALPHA
    RROOTI = zeros(NLAYER)
    ALPHA  = zeros(NLAYER)
    for i = 1:NLAYER
        if (p_fT_RELDEN[i] < 0.00001 || p_RTLEN < 0.1)
            # no roots in layer
            RROOTI[i] = 1E+20
            ALPHA[i]  = 1E+20
        else
            # root resistance for layer
            RROOTI[i] = (p_RPLANT - RXYLEM) / RTFRAC[i]
            # rhizosphere resistance for layer
            RTDENI = RTFRAC[i] * 0.001 * p_RTLEN / D[i] # .001 is (mm/mm2)/(m/m2) conversion
            DELT = p_PI * p_RTRAD ^ 2 * RTDENI
            ALPHA[i] = (1 / (8 * p_PI * RTDENI)) * (DELT - 3 - 2 * (log(DELT)) / (1 - DELT))
            ALPHA[i] = ALPHA[i] * 0.001 * p_RHOWG / D[i] # .001 is MPa/kPa conversion
        end
    end
    return (RXYLEM, RROOTI, ALPHA)
end


"""

    TBYLAYER()

Allocate transporation among soil layers.

Ecoshift:
TBYLAYER determines the rate at which liquid water can be supplied to transpiring leaves,
compares this rate with the potential transpiration rate, sets actual transpiration equal to
the lesser of the two, and then allocates the transpiration among soil layers. This routine
is based on the model of Federer (1979), which has been widely used, e.g. Wetzel and Chang
(1987), Levine and Salvucci (1999).

The routine requires iteration when outflow from roots is prevented for layers in which
xylem potential ψx is greater than soil water potential ψti(Fig. EVP-1).

The resistance to uptake from each layer, ri, is

ri = rri + rsi = rri + αi / Ki

(Fig. EVP-1) where rri (RROOTI) and αi (ALPHAi) are variables from subroutine PLNTRES and Ki
is the rhizosphere conductivity of the layer (Cowan 1965, Federer 1979). Estimating
rhizosphere conductivity requires an iterative solution as described in Federer (1979), so
BROOK90 uses the Ki (KKi) of the bulk soil instead. This tends to underestimate the
rhizosphere resistance, but the error is probably unimportant unless rp is quite small. The
total root resistance to uptake is

rt = 1 / Σ ( 1 / ri ) .

The transpiration rate for each layer, Ti , is

Ti = ( ψti - ψx ) / ri

where ψti is the total water potential for the layer (PSITI) , and ψx is the xylem
potential. The total transpiration from all layers is

T = ΣTi = Σ( ψti / ri ) - ψx Σ ( 1 / ri ) = ( ψt - ψx ) / rt

where ψt is a weighted mean soil water potential defined as

ψt = ψt Σ (ψti / ri ) .

Eliminating ψx from the Ti and T equations gives

Ti = ( ψti - ψt + rt T ) / ri

which is the equation used to distribute T among layers.

In BROOK90 actual transpiration rate, T, is the lesser of potential transpiration rate, P,
and the maximum possible rate of water uptake from the soil, S. Following Federer (1979) and
Lynn and Carlson (1990), BROOK90 assumes that the S obtains when the leaf potential, ψ, is
the critical potential at which stomates close, ψc , which is an input parameter (PSICR).
The assumption that the T is the lesser of the P and S is equivalent to an assumption that
stomatal opening is not limited by leaf water potential, ψ, until the critical potential,
ψc, is reached. The stomata are then assumed to close as much as necessary to maintain ψ =
ψc. This abrupt switch from a demand-limited to supply-limited regime is oversimplified, but
is convenient for modeling. This type of behavior induces a flat-topped diurnal
transpiration curve (Lynn and Carlson 1990). The critical potential, ψc, represents the
maximum suction that the plant can exert to get water from the soil, the minimum soil-water
potential that the plant can induce, and the lower limit of soil water availability.

The water supply rate, S, is the potential difference between leaf water potential, ψ, and
the xylem potential, ψx, divided by the xylem resistance, rx , when ψ is equal to the
critical potential, ψc (see Fig EVP-1). With allowance for the gravity potential between the
ground and the canopy zero-plane displacement height, d (DISP), this is

S = ( ψx - ψc - ρw g d ) / rx

S below ground must be T = ( ψt - ψx ) / rt from above, and equating the two S equations to
eliminating ψx, gives

S = ( ψt - ψc - ρw g d) / ( rt + rx )

which is assumed to be constant throughout a day.

 Fig. EVP-2. Transpiration for the day, T' (shaded area), as the lesser of a constant water
 supply rate, S, and a half-sine potential transpiration, P'.

Following Federer (1982), BROOK90 assumes that the daytime potential transpiration rate, P',
varies as a half sine wave. The actual transpiration rate, T' , is the lesser of S and P'
(Fig. EVP-2). The average value of T' over the daytime, T, is

T = P ( 1 + R cos-1 R - sin (cos-1 R) )  R < 1

T = P       R >= 1

where

R = 2 S / π P

where P is the average of P' over the daytime. (Notation here differs slightly from Federer
(1982) who used S' for S, D/d for P, and T/d for T where d is daylength.) At night, P is
assumed constant and T is the lesser of S and P.

Normally, all layers with roots are included in the above calculations. However, when some
layers are wet and others are dry, it is possible for ψti (PSITI) in one or more layers to
be lower than ψx so that the roots in those layers are releasing water to the soil. The
question of whether such outflow occurs or is somehow prevented by the roots is
controversial (Hunt et al. 1991). More recent work shows that outflow does occur (Dawson
1993). Richards and Caldwell (1987) name the process "hydraulic lift" because it moves water
upward from wetter deeper soil layers through the plant to shallow drier layers. In BROOK90
outflow from roots is prevented when the parameter NOOUTF is set to 1, and is allowed if
NOOUTF is 0. In general, BROOK90 usually will work better when NOOUTF = 1; this is the
default. When NOOUTF = 1 and any Ti is negative, the layer with the most negative Ti is
eliminated and new values of rt , ψt , T, and Ti are obtained. If any Ti are still negative,
the elimination process is repeated. This elimination procedure causes transpiration from a
layer to cease when its potential is still greater than PSICR.
"""
function TBYLAYER(J, p_fu_PTR, p_fu_DISPC, p_fT_ALPHA, p_fu_KK, p_fT_RROOTI, p_fT_RXYLEM, u_aux_PSITI, NLAYER, p_PSICR, NOOUTF)

    FLAG = zeros(NLAYER)
    for i = 1:NLAYER
        if p_fT_RROOTI[i] > 1E+15
            # This layer has no roots
            FLAG[i] = 1
        elseif (NOOUTF && u_aux_PSITI[i] / 1000. <= p_PSICR)
            # This layer has no outflow from roots to soil
            FLAG[i] = 1
        else
            # This layer has roots connected to soil
            FLAG[i]= 0
        end
    end

    # Compute ATR and ATRANI
    # top of loop for recalculation of transpiration if more layers get flagged
    RI = zeros(NLAYER)
    ATRANI=zeros(NLAYER)
    ATR = 0.0
    while true
        # start loop with NEGFLAG = 0
        NEGFLAG = false

        ###
        SUM = 0
        for i = 1:NLAYER
            if (FLAG[i] == 0)
                RI[i] = p_fT_RROOTI[i] + p_fT_ALPHA[i] / p_fu_KK[i]
                SUM = SUM + 1.0 / RI[i]
            else
                ATRANI[i] = 0.0
            end
        end
        if SUM < 1E-20
                ATR = 0.
                PSIT = -1e10

                break
        else
                RT = 1.0 / SUM
        end

        ###
        # weighted mean soil water potential
        PSIT = 0
        for  i = 1:NLAYER
                if FLAG[i] == 0
                    PSIT = PSIT + RT * u_aux_PSITI[i] / RI[i]
                end
        end

        # soil water supply rate, assumed constant over day
        SUPPLY = (PSIT / 1000 - p_PSICR - p_RHOWG * p_fu_DISPC) / (RT + p_fT_RXYLEM)
        # transpiration rate limited by either PTR or SUPPLY

        if J == 1
        # daytime, PTR is average of a half sine over daytime
            R = (2 / p_PI) * (SUPPLY / p_fu_PTR)
            if R <= 0
                ATR = 0.
            elseif R < 1
                ATR = p_fu_PTR * (1 + R * acos(R) - sin(acos(R)))
            else
                ATR = p_fu_PTR
            end
        else
        # nighttime, PTR assumed constant over nighttime
            if (SUPPLY <= 0) || (p_fu_PTR <= 0)
                ATR = 0.
            else
                ATR = min(SUPPLY, p_fu_PTR)
            end
        end

        # distribute total transpiration rate to layers
        for  i = 1:NLAYER
            if FLAG[i] == 1
                ATRANI[i] = 0
            else
                ATRANI[i] = ((u_aux_PSITI[i]- PSIT) / 1000 + RT * ATR) / RI[i]

                # check for any negative transpiration losses
                if ATRANI[i] < -0.000001
                    NEGFLAG = true
                end
            end
        end

        ###
        if NOOUTF && NEGFLAG
            # find layer with most negative transpiration and omit it
            IDEL = 0
            TRMIN = 0
            for i = 1:NLAYER
                if (ATRANI[i] < TRMIN)
                    TRMIN = ATRANI[i]
                    IDEL = i
                end
            end
            FLAG[IDEL] = 1
        # repeat main loop with flagged layers excluded
        else
            break
        end
    end

    return (ATR, ATRANI)
end


"""
    INTER(p_fT_RFAL, p_fu_PINT, p_fu_LAI, p_fu_SAI, p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS, p_DTP, u_INTR)

Compute rain catch rate (interception) and evaporation rate of intercepted rain in mm/d.

Rain interception, used when p_NPINT > 1.

# Arguments
- `p_fT_RFAL`: rainfall rate, mm/d
- `p_fu_PINT`: potential interception rate, mm/d
- `p_fu_LAI`: projected leaf area index, m2/m2
- `p_fu_SAI`: projected stem area index, m2/m2
- `p_FRINTL`: intercepted fraction of p_fT_RFAL per unit p_fu_LAI
- `p_FRINTS`: intercepted fraction of p_fT_RFAL per unit p_fu_SAI
- `p_CINTRL`: maximum interception storage of rain per unit p_fu_LAI, mm
- `p_CINTRS`: maximum interception storage of rain per unit p_fu_SAI, mm
- `p_DTP`: precipitation interval time step, d
- `u_INTR`: intercepted rain storage, mm,

# Ecoshift:
Older studies of rain and snow interception regressed throughfall on precipitation, but such
interpretation ignored the fact that energy supply rather than water supply may limit
interception and also ignores storm duration/intensity and interstorm interval. Detailed
simulation models of rain interception over time through a storm have been developed (Rutter
et al. 1972) but these are too complex to include in a hydrologic model like BROOK90. Much
less is known about the complicated process of evaporation of intercepted snow.

Subroutine INTER accounts in the simplest way for the concepts of catch rate, evaporation
rate, and canopy capacity. The same algorithm is applied to both rain and snow, which are
considered to behave independently with respect to their interception. INTER is used when
precipitation data are input more than once a day in a precip. interval file (PINT > 1). If
only daily precipitation is input, then the modified procedure of INTER24 is used.

The conservation of mass equation for rain interception can be written as

dS/dt = C - I - D

where S is the amount of water stored on the canopy (mm), C is the catch rate, or rate of
water input to the canopy, I is the rate of evaporation of intercepted water, and D is the
drip rate, or rate of transfer of liquid water to the ground. The same equation applies to
snow or mixed snow and rain, when any solid-liquid phase change is ignored and D includes
all rain or snow blowing or falling from canopy to ground.

BROOK90 ignores D by defining C as a net catch rate (C - D), or only the portion of the
catch that will sooner or later evaporate, so, from the Flow Chart,

d INTR / dt = RINT - IRVP

for rain, and

d INTS / dt = SINT - ISVP

for snow, where INTR and INTS are the canopy storages, RINT and SINT are the net catch
rates, and IRVP and ISVP are the evaporation rates, for rain and snow respectively.

BROOK90 assumes that interception catch rates, RINT and SINT, are a constant fraction of
rainfall or snowfall until the canopy reaches a storage capacity. Until the capacity is
reached, RINT and SINT are assumed to be linear functions of LAI and SAI, so that

RINT = (FRINTL * LAI + FRINTS * SAI) * RFAL

and

SINT = (FSINTL * LAI + FSINTS * SAI) * SFAL

where RFAL and SFAL are rainfall rate and snowfall rate as determined from subroutine
SNOFRAC, FRINTL and FSINTL are the catch fraction per unit LAI for rain and snow,
respectively, and FRINTS and FSINTS are the catch fraction per unit SAI for rain and snow,
respectively.

The canopy has capacities or maximum values of INTR and INTS that depend on LAI and SAI. In
BROOK90 these dependencies are assumed linear. The parameters CINTRL and CINTRS are the
capacities for intercepted rain per unit LAI and SAI respectively, so that INTRMX, the
capacity for rain, is

INTRMX = CINTRL * LAI + CINTRS * SAI.

For snow,

INTSMX = CINTSL * LAI + CINTSS * SAI,

and the capacity parameters are generally larger than for rain. The eight interception
parameters, FRINTL, FRINTS, FSINTL, FSINTS, CINTRL, CINTRS, CINTSL, and CINTSS, only control
interception loss in small storms; interception loss in large storms is controlled by the
evaporation rate of intercepted water (PINT) and the storm intensity and duration.

The rate at which intercepted water evaporates (PINT) is calculated from the
Shuttleworth-Wallace equations by calling subroutine SWPE (Section PET) with the canopy
resistance rc = 0. The soil surface resistance (RSS) is not reduced for the PINT
calculation. The MSBDAYNIGHT routine does this separately for daytime and for nighttime
weather variables and the results are weighted by daylength (DAYLEN) to produce PINT. PINT
is considered to be constant throughout the daily time step; its actual diurnal variation is
ignored.

The canopy is considered to be either completely wetted or completely dry. Partial canopy
wetting and drying is not treated in BROOK90, though it is a key component of specific
models of the interception process (Rutter et al. 1972). Subroutine INTER determines the
actual catch rate (RINT or SINT) and the actual evaporation rate (IRVP or ISVP) for the
precipitation time step in the three cases that the canopy dries during the timestep, the
canopy wets but does not reach capacity, and the canopy reaches capacity. The routine
appropriately handles the case of a wet canopy with decreasing capacity because of
decreasing LAI or SAI by allowing RINT or SINT to be negative.
"""
function INTER(p_fT_RFAL, p_fu_PINT, p_fu_LAI, p_fu_SAI, p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS, p_DTP, u_INTR)
    # p_fT_RFAL rainfall rate, mm/d
    # p_fu_PINT potential interception rate, mm/d
    # p_fu_LAI  projected leaf area index, m2/m2
    # p_fu_SAI  projected stem area index, m2/m2
    # p_FRINTL  intercepted fraction of p_fT_RFAL per unit p_fu_LAI
    # p_FRINTS  intercepted fraction of p_fT_RFAL per unit p_fu_SAI
    # p_CINTRL  maximum interception storage of rain per unit p_fu_LAI, mm
    # p_CINTRS  maximum interception storage of rain per unit p_fu_SAI, mm
    # p_DTP     precipitation interval time step, d
    # u_INTR    intercepted rain storage, mm,

    # maximum RINT, mm/d
    CATCH = min(1, p_FRINTL * p_fu_LAI + p_FRINTS * p_fu_SAI) * p_fT_RFAL
    # maximum canopy storage for rain, mm
    INTRMX = p_CINTRL * p_fu_LAI + p_CINTRS * p_fu_SAI
    # first approximation to new canopy storage (INTR)
    NEWINT = u_INTR + (CATCH - p_fu_PINT) * p_DTP

    if NEWINT > 0
        # canopy is wet throughout DTP
        aux_du_IRVP = p_fu_PINT
        if NEWINT > INTRMX
            # canopy capacity is reached
            aux_du_RINT = p_fu_PINT + (INTRMX - u_INTR) / p_DTP
            # RINT can be negative if INTR exists and LAI or SAI is decreasing over time
        else
            # canopy capacity is not reached
            aux_du_RINT = CATCH
        end
    else
        # canopy dries during interval or stays dry
        aux_du_RINT = CATCH
        aux_du_IRVP = (u_INTR / p_DTP) + CATCH

        # IRVP is < PINT
    end

    # RINT     - rain catch rate, mm/d
    # IRVP     - evaporation rate of intercepted rain, mm/d
  return (aux_du_RINT, aux_du_IRVP)
end

"""
    INTER(p_fT_RFAL, p_fu_PINT, p_fu_LAI, p_fu_SAI, p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS, p_DTP, u_INTR)

Compute rain catch rate (interception) and evaporation rate of intercepted rain in mm/d.

Rain interception with duration in hours, used when p_NPINT = 1. Same routine is used for
snow interception, with different calling variables.

# Arguments:
- `p_fT_RFAL`: 24-hour average rainfall rate, mm/d
- `p_fu_PINT`: potential interception rate, mm/d
- `p_fu_LAI`: projected leaf area index, m2/m2
- `p_fu_SAI`: projected stem area index, m2/m2
- `p_FRINTL`: intercepted fraction of p_fT_RFAL per unit p_fu_LAI
- `p_FRINTS`: intercepted fraction of p_fT_RFAL per unit p_fu_SAI
- `p_CINTRL`: maximum interception storage of rain per unit p_fu_LAI, mm
- `p_CINTRS`: maximum interception storage of rain per unit p_fu_SAI, mm
- `p_DURATN`: average storm duration, hr
- `u_INTR`: intercepted rain storage, mm,
- `MONTHN`: Month of the year

# Ecoshift:
"
Subroutine INTER24 - daily interception
Proper representation and integration of the interception process is a problem for
hydrologic models that use a daily interval for precipitation input (p_NPINT = 1), because the
storm duration is not known. For a brief, intense storm, the canopy wets once and the
interception loss is limited primarily by canopy capacity. For a low intensity, all day
storm, the canopy stays wet and the interception loss is limited primarily by the potential
interception, PINT. This problem is worst when only daily precipitation is known, and
decreases as precipitation is given at shorter intervals.

Subroutine INTER24 was developed because the use of subroutine INTER for daily precipitation
consistently produced too much interception. INTER24 is a modification of INTER that loops
through the procedure every hour,using the PINT rate for each hour. DURATN is a parameter
that specifies the average hourly duration of precipitation for each month of the year.
INTER24 truncates DURATN to the next lower even integer, and then centers the "storm" on
noon. Thus if DURATN is input as 7.5, the daily precipitation is assumed to occur at a
constant rate from time 0900 to 1500. Centering on noon is only used to see how much
interception carries over into the next day. The algorithm for each hourly loop is the same
as for INTER, except that rates are in mm/hr and amounts are summed over the day. The
interception catch rate (RINT or SINT), and the evaporation rate (IRVP or ISVP) are returned
to MSBPREINT as average rates over the day.

To determine appropriate values of DURATN I examined hourly precipitation data for 4 years
(one year at Hubbard Brook) from the SAMSON data set. Averaging the number of hours per day
of precipitation of 0.02 inch (0.5 mm) or greater over days with such precipitation gave the
following results after a little smoothing

                        J   F   M   A   M   J   J   A   S   O   N   D
    San Juan PR         3   2   2   2   2   2   2   3   3   3   3   3
    Atlanta GA          5   5   5   5   4   4   3   3   4   4   5   6
    Caribou ME          4   4   5   5   4   4   4   4   4   6   6   5
    Madison WI          4   4   5   3   3   2   3   3   4   4   5   5
    Lake Charles LA     5   4   3   3   3   3   2   2   3   3   4   5
    Phoenix AZ          4   4   4   4   4   2   2   2   2   2   4   4
    Rapid City SD       3   3   3   4   4   3   2   2   2   2   4   4
    Tacoma WA           6   6   5   4   4   4   4   4   4   4   6   6
    Fairbanks AK        3   3   4   4   4   4   3   3   4   4   4   3
    Hubbard Brook NH    5   5   5   4   4   4   4   4   4   5   5   5

Apparently a default DURATN of 4 hours is appropriate for all months anywhere in the U.S.
"
"""
function INTER24(p_fT_RFAL, p_fu_PINT, p_fu_LAI, p_fu_SAI, p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS, p_DURATN, u_INTR, MONTHN)
    # p_fT_RFAL  24-hour average rainfall rate, mm/d
    # p_fu_PINT  potential interception rate, mm/d
    # p_fu_LAI projected leaf area index, m2/m2
    # p_fu_SAI projected stem area index, m2/m2
    # p_FRINTL  intercepted fraction of p_fT_RFAL per unit p_fu_LAI
    # p_FRINTS  intercepted fraction of p_fT_RFAL per unit p_fu_SAI
    # p_CINTRL  maximum interception storage of rain per unit p_fu_LAI, mm
    # p_CINTRS  maximum interception storage of rain per unit p_fu_SAI, mm
    # p_DURATN  average storm duration, hr
    # u_INTR  intercepted rain storage, mm,
    # MONTHN  Month of the year

    IHD = Int(floor((p_DURATN[MONTHN] + 0.1) / 2.0))
    DTH = 1                     # time step, = 1 hr
    # Define catch along day
    CATCH = fill(NaN, 24)       # maximum RINTHR, mm/hra
    for i = 0:23                # hour, 0 to 23
        if i < (12 - IHD) || i >= (12 + IHD)
            # before or after rain
            CATCH[i+1] = 0
        else
            # during rain, mm/hr is rate in mm/d divided by hr of rain/d
            CATCH[i+1] = min(1, (p_FRINTL * p_fu_LAI + p_FRINTS * p_fu_SAI)) * p_fT_RFAL / (2 * IHD)
        end
    end

    # Define RINTHR and IRVPHR
    INTRMX = p_CINTRL * p_fu_LAI + p_CINTRS * p_fu_SAI  # maximum canopy storage for rain, mm
    INTRNU = fill(NaN, 24)                              # canopy storage at end of hour, mm

    RINTHR = fill(NaN, 24) # rain catch rate for hour, mm/hr
    IRVPHR = fill(NaN, 24) # evaporation rate for hour, mm/hr

    # Set INTRNU starting value for first loop
    INTRNU[1]= u_INTR

    for i = 1:24
        # INTRNU - canopy storage at end of hour (mm)
        NEWINT = INTRNU[i] + (CATCH[i] - p_fu_PINT / 24) * DTH # first approximation to INTRNU, mm

        if (NEWINT > 0.0001)
            # canopy is wet throughout hour, evap rate is PINT
            IRVPHR[i] = p_fu_PINT / 24
            if (NEWINT > INTRMX)
                # canopy capacity is reached
                RINTHR[i] = IRVPHR[i] + (INTRMX - INTRNU[i]) / DTH
                # INTRMX - INTRNU can be negative if LAI or SAI is decreasing over time
            else
                # canopy capacity is not reached
                RINTHR[i] = CATCH[i]
            end
        else
            # canopy dries during hour or stays dry
            IRVPHR[i] = INTRNU[i] / DTH + CATCH[i]
            RINTHR[i] = CATCH[i]
            # IRVPHR for hour is < PI/24
        end

        # Set INTRNU starting value for next loop
        if i < 24
            INTRNU[i+1] = INTRNU[i] + (RINTHR[i] - IRVPHR[i]) * DTH
        end
    end

    aux_du_RINT = sum(RINTHR .* DTH)# = SMINT    per 1 d # daily accumulated actual catch, mm
    aux_du_IRVP = sum(IRVPHR .* DTH)# = SMVP     per 1 d # daily accumulated actual evaporation, mm

    # RINT     -  rain catch rate, mm/d
    # IRVP     -  evaporation rate of intercepted rain, mm/d
    return (aux_du_RINT, aux_du_IRVP)
end



end