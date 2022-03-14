@doc raw"""
# Snow accumulation and melt

Text copied from Ecoshift on module SNO:

"
Simulation of snow accumulation and melt is a complex subject (U.S. Army Corps Eng. 1956,
Anderson and Crawford 1964, Colbeck and Ray 1979). Detailed expressions for the energy
balance of snow surfaces have been developed (Anderson 1976, Dingman 1994) but they are
not generalized to all cover types and are too complex for BROOK90. Leavesley and
Striffler (1979) give an energy balance model that includes radiation melt and the effects
of canopy cover on it, but ignores convection-condensation melt. More complex algorithms
could be developed, but Colbeck et al. (1979) say "the energy exchange processes between
snow and a forest cover are not well enough understood to allow detailed modelling of the
melt process through use of the energy equation." The energy balance is made complicated
because of the heat of fusion in the water-ice phase change, and because the surface vapor
pressure for melting snow is fixed. Application of the Shuttleworth-Wallace two layer
approach to snow under sparse canopies remains in the future.

BROOK90 therefore falls back on the classic degree-day method for estimating snow energy
balance. Anderson (1979, p. 342) states that "Air temperature (ambient) is an adequate
index to surface energy exchange in most cases." This is not a perfect solution as
Anderson points out three cases in which air temperature fails: 1) warm temperatures with
little wind causes overestimates of melt, 2) high dewpoint with high wind causes
underestimates, and 3) low temperatures with clear sky and ripe snow causes
underestimates. A modified degree-day method that incorporates solar radiation and wind
speed could be added but would require development for sparse canopies. Another
improvement would separate the day into daytime and nighttime as is done in BROOK90 for
evaporation. But BROOK90 currently uses only the mean daily temperature (TA) for the snow
energy balance.

The "water equivalent" of snow (SNOW, mm) is the depth of water a snowpack would produce
if it were all melted; this is the BROOK90 variable that represents the snowpack. The
actual depth of snow, assuming a constant snow density (SNODEN), is used only to calculate
the amount of the canopy above the snow in subroutine CANOPY. Variable snow density (mass
per unit volume) is not simulated in BROOK90. When the snow is colder than 0°C, it has a
"cold content" (CC, MJ/m2), which is the energy needed to warm the snow to 0°C. When the
snow is at 0°C, part of SNOW can be liquid water (SNOWLQ, mm). The maximum liquid water
that can be retained by the snowpack without draining is a constant fraction (MAXLQF) of
SNOW; CC and SNOWLQ are always initialized as zero in BROOK90, so any initial SNOW is
considered to be at 0°C with no liquid water.

Groundmelt is snowmelt at the bottom of a snowpack; it occurs because of heat conduction
from the soil whenever the soil is unfrozen. A constant groundmelt rate (GRDMLT, mm/d) is
an constant parameter in BROOK90 and is applied whenever there is snow on the ground. The
possibilities of frozen soil or variable groundmelt are not considered.

Snowmelt (SMLT, mm/d) is the sum of groundmelt and drainage of excess liquid water from
the snowpack. Drainage occurs only after the snowpack is both isothermal at 0°C and is
holding the maximum possible liquid water; the snowpack is then "ripe". The gains and
losses of liquid water by the snowpack, including the refreezing of rain on cold snow, are
handled in the somewhat complicated subroutine SNOWPACK.

BROOK90 assumes that the snowpack is always isothermal. In reality, large and variable
temperature gradients can exist in thick snowpacks; simulating these is beyond the scope
of BROOK90. The snowpack temperature (TSNOW) at the beginning of the day is calculated in
MSBSETVARS from the cold content

TSNOW = -CC / (CVICE * SNOW)

where CVICE is the heat capacity of ice (0.00192 MJ m-2 mm-1 K-1) (Leavesley and
Striffler, 1979). TSNOW is used both in calculating snow evaporation (SNVP) and snow
energy flux (SNOEN).
"
"""
module SNO # SNOW ACCUMULATION AND MELT

using ..CONSTANTS: p_CVICE, p_CVLQ, p_LF, p_WTOMJ, p_LS, p_CPRHO, p_GAMMA # https://discourse.julialang.org/t/large-programs-structuring-modules-include-such-that-to-increase-performance-and-readability/29102/5
using ..PET: ESAT, SWGRA

export SNOFRAC, SNOVAP, SNOENRGY


"""
    SNOFRAC(p_fT_TMAX, p_fT_TMIN, p_RSTEMP)

Separate RFAL from SFAL.

# Ecoshift
"
Separation of daily precipitation into snow or rain is a major problem in hydrologic
modeling. For instance, if the wrong precipitation form is chosen in December, simulated
streamflow from that day's precipitation could be shifted from December to April or vice
versa, a rather significant effect! BROOK90 uses both daily maximum (TMAX) and daily minimum
(TMIN) temperatures to allow days on which mixed rain and snow falls. This reduces the
potential error from making the wrong choice. The algorithm seems to have been stated first
by Willen and Shumway (1971). When TMAX for the day is greater than the parameter RSTEMP and
TMIN is less than RSTEMP, the fraction of precipitation as snow, SNOFRC, is

SNOFRC = (RSTEMP - TMIN) / (TMAX - TMIN)

where RSTEMP is the "base" temperature for the rain-snow transition. When TMAX < RSTEMP,
SNOFRC = 1; when TMIN > RSTEMP, SNOFRC = 0. The default value of RSTEMP is -0.5°C because
that seems to work best at Hubbard Brook. If precipitation is input more than once a day,
the same SNOFRC is used for all precipitation intervals.
"
"""
function SNOFRAC(p_fT_TMAX, p_fT_TMIN, p_RSTEMP)
    if (p_fT_TMIN >= p_RSTEMP)
        p_fT_SNOFRC = 0.0
    elseif (p_fT_TMAX < p_RSTEMP)
        p_fT_SNOFRC = 1.0
    else
        p_fT_SNOFRC = 1.0 - (p_fT_TMAX - p_RSTEMP) / (p_fT_TMAX - p_fT_TMIN)
    end
    return p_fT_SNOFRC
end

"""
    SNOVAP(p_fu_TSNOW, p_fu_TA, p_fT_EA, UA, p_fu_ZA, p_fu_HEIGHT, p_fu_Z0, p_fu_DISP,
    p_fu_Z0C,p_fu_DISPC, p_fu_Z0GS, p_LWIDTH, p_RHOTP, p_NN, p_fu_LAI, p_fu_SAI, p_KSNVP)

Compute potential snow evaporation (mm/d).


# Ecoshift
"
Evaporation rate from the snowpack, and its negative, condensation, are evaluated
using the aerodynamic flux equation

E = (cp ρ / g Ls rw ) (e0 - ea) / (raa + ras)

where ea is the vapor pressure of the air, e0 is the surface vapor pressure, and raa and ras
are the Shuttleworth-Wallace resistances described in section PET. The constants cp ρ
(CPRHO), γ (GAMMA), and the latent heat of sublimation Ls ρw (LS) are constant. BROOK90
assumes that the snowpack is always isothermal and that its temperature does not change
diurnally. When the snowpack temperature (TSNOW) is less than 0°C, the surface vapor
pressure, e0, is the saturated vapor pressure at TSNOW and is obtained by calling subroutine
ESAT. When TSNOW is 0°C, e0 is 0.61 kPa; use of Ls instead of the latent heat of
vaporization, Lv, in this case is slightly wrong. The vapor pressure ea (EA) is the input
vapor pressure at reference height za . The resistances raa and ras are obtained from
subroutine SWGRA using the daily average wind speed (UA). The value of E in mm/d returned
from SNOVAP is called PSNVP because it can be reduced in SNOWPACK if snow disappears.

For evaporation from snow in the open, U.S. Army Corps of Engineers (1956) gives

E = 1.9 ua (e0 - ea)

for evaporation in the open. This yields E = 0.6 mm when ua = 3 m/s and the vapor pressure
difference is 0.1 kPa. The two E equations are the same when ras = 0, za - d = 2 m, and z0 =
1 mm. Subroutine SWGRA thus gives the appropriate raa for snow in the open.

Colbeck et al. (1979) state "Evaporation from the snow in a forest has received a great deal
of attention, with many investigators concluding that it is small.... Although there are
many reports of high evaporative losses from forests, these have not been verified from heat
balance considerations." Generally, literature values are around 0.5 mm/d in the open, with
monthly values around 10 mm and annual values of 20 mm or more. Anderson (1976) gives 15-20
mm annually for Sleepers River, VT.

The modified Shuttleworth and Gurney resistance formulations in subroutine SWGRA for a
leafless tall canopy give raa about 3 s/m and ras about 40 s/m for a weather station wind
speed of 3 m/s. A common vapor pressure difference of 0.1 kPa then gives a very high
evaporation of 1.6 mm/d or 50 mm/ month. The problem may be that the resistances calculated
by SWGRA are too small, either because of too large a roughness parameter for leafless
deciduous forests, or because stability effects are ignored. To fix the problem BROOK90
includes KSNVP, which is an arbitrary constant multiplier of PSNVP. Use KSNVP = 1.0 for open
open areas or short canopies normally covered by snow; but for forest, KSNVP of 0.3 gives
more reasonable values of SNVP. More work is obviously needed on amount and theory of snow
evaporation under forests.

Note that although evaporation and condensation of water are simulated in SNOVAP, the
accompanying latent transfer is not simulated. The snow energy balance in subroutine
SNOENRGY is (unfortunately) decoupled from the snow evaporation-condensation process.
"
"""
function SNOVAP(p_fu_TSNOW, p_fu_TA, p_fT_EA, UA, p_fu_ZA, p_fu_HEIGHT, p_fu_Z0, p_fu_DISP,
    p_fu_Z0C, p_fu_DISPC, p_fu_Z0GS, p_LWIDTH, p_RHOTP, p_NN, p_fu_LAI, p_fu_SAI, p_KSNVP)
    # ignores effect of interception on p_fu_PSNVP or of p_fu_PSNVP on p_fu_PTRAN
    if (p_fu_TSNOW > -0.1)
        ESNOW = 0.61
    else
        # snow surface vapor pressure saturated at lower of p_fu_TA and p_fu_TSNOW
        ESNOW, dummy = ESAT(min(p_fu_TA, p_fu_TSNOW))
    end
    RAA, RAC, RAS = SWGRA(UA, p_fu_ZA, p_fu_HEIGHT, p_fu_Z0, p_fu_DISP, p_fu_Z0C,
                          p_fu_DISPC, p_fu_Z0GS, p_LWIDTH, p_RHOTP, p_NN, p_fu_LAI, p_fu_SAI)

    p_fu_PSNVP = (p_WTOMJ / p_LS) * (p_CPRHO / p_GAMMA) * (ESNOW - p_fT_EA) / (RAA + RAS)
    # fix for PSNVP overestimate
    p_fu_PSNVP = p_KSNVP * p_fu_PSNVP

    return p_fu_PSNVP # potential snow evaporation (mm/d)
end

"""
    SNOENRGY(p_fu_TSNOW, p_fu_TA, p_fT_DAYLEN, p_CCFAC, p_MELFAC, p_fT_SLFDAY, p_fu_LAI,
    p_fu_SAI, p_LAIMLT, p_SAIMLT)

Compute energy flux density to snow surface, MJ m-2 d-1.

# Ecoshift
"
The energy flux density to the snow surface (SNOEN, MJ m-2 d-1) is calculated in subroutine
SNOENRGY independently of precipitation for the day.

When TA is <= 0°C, SNOEN is the energy used to heat or cool the snowpack

SNOEN = CCFAC * 2 * DAYLEN * (TA - TSNOW)

where CCFAC is an input parameter, and DAYLEN is the daytime fraction of a day. CCFAC is the
below-zero degree-day factor for a day with a daylength of 0.5 d. Anderson (1973) pointed
out that this degree-day factor appears to vary seasonally. In BROOK90, this seasonality is
incorporated by using 2 * DAYLEN as a multiplier in the above equation. BROOK90 is not very
sensitive to CCFAC unless it is close to 0; larger values of CCFAC make the snow melt later
because the snowpack cools more. When CCFAC = 0, snow temperature is always 0°C and there is
never any cold content.

When TA is greater than 0°C, energy is provided that can melt snow, and SNOEN is calculated
differently. The energy supply rate, SNOEN, is then

SNOEN = MELFAC * 2 * DAYLEN * SLFDAY * TA * exp(-LAIMLT * LAI - SAIMLT * SAI)

where MELFAC is the melting degree-day factor for a day with a daylength of 0.5 d and no
canopy, SLFDAY is the ratio of potential insolation on the slope to that on a horizontal
surface, and the input parameters LAIMLT and SAIMLT express the dependence of SNOEN on leaf
area index (LAI) and stem area index (SAI). MELFAC uses 0°C as a base and is zero for TA
below this. Inclusion of SLFDAY in the MELT equation arises from an assumption that
radiation melt plays an important role. If this is not so, then SLFDAY multiplier could be
omitted, and snowmelt would not depend on slope-aspect. The functional forms of the SAI and
LAI dependencies are based on the somewhat arbitrary curves used by Federer and Lash (1978).
"
"""
function SNOENRGY(p_fu_TSNOW, p_fu_TA, p_fT_DAYLEN, p_CCFAC, p_MELFAC, p_fT_SLFDAY, p_fu_LAI, p_fu_SAI, p_LAIMLT, p_SAIMLT)
    # snow surface energy balance
    if (p_fu_TA <= 0)
        # snow warms or cools proportional to snow-air temperature difference
        p_fu_SNOEN = p_CCFAC * 2.0 * p_fT_DAYLEN * (p_fu_TA - p_fu_TSNOW)
    else
        # energy input proportional to TA, modified by cover, slope-aspect, and daylength
        p_fu_SNOEN = p_MELFAC * 2.0 * p_fT_DAYLEN * p_fu_TA * exp(-p_SAIMLT * p_fu_SAI) * exp(-p_LAIMLT * p_fu_LAI) * p_fT_SLFDAY
    end

    return p_fu_SNOEN # energy flux density to snow surface, MJ m-2 d-1
end

"""
    SNOWPACK(p_fu_RTHR, p_fu_STHR, p_fu_PSNVP, p_fu_SNOEN, u_CC, u_SNOW,
    u_SNOWLQ, p_DTP, p_fu_TA, p_MAXLQF, p_GRDMLT,
    p_CVICE, p_LF, p_CVLQ)

Update status of snowpack.

# Ecoshift
"
In each precipitation interval, throughfall of rain (RTHR) and throughfall of snow (STHR)
are calculated and subroutine SNOWPACK is entered if there is STHR or if there is SNOW. This
subroutine adds throughfall to the snowpack, subtracts groundmelt, snow evaporation, and
drainage of liquid water, and calculates the new cold content (CC) and liquid water content
(SNOWLQ) of the snowpack. There are a number of different ways to program all this adding
and subtracting of energy and water. The program flow of SNOWPACK has been selected to make
the many realizable situations as clear as possible. Much shorter algorithms could have been
used, but at the expense of clarity. Any alterations to this routine must be made carefully,
keeping all the possibilities in mind. However, unlike routine SNOENRGY, the content of the
SNOWPACK routine is essentially standard for all snow models and alterations generally
should not be necessary.

In SNOWPACK, snow throughfall (STHR) is first added to SNOW. If TA is <0°C, the cold content
of the new snow is increased by

CC = CC - CVICE * TA * STHR * DTP

where DTP is the precipitation interval time step and CVICE is the volumetric heat capacity
of ice. If this addition of cold snow causes both cold content (CC) and liquid water
(SNOWLQ) to coexist, then liquid water is refrozen; CC is used to refreeze part or all of
the liquid.

Groundmelt (GRDMLT) and snow evaporation-condensation are dealt with next. In the
precipitation time interval (DTP), the fraction of the snowpack that melts by groundmelt and
evaporates is

FRAC = (GRDMLT + PSNVP) * DTP / SNOW.

If FRAC is > 1 then all the snow melts and evaporates at rates proportional to GRDMLT and
PSNVP. If FRAC is < 1, SNVP is equal to potential snow evaporation (PSNVP) and GRDMLT drains
from the snowpack as snowmelt (SMLT). An assumption is made that evaporation and groundmelt
remove any liquid water and cold content that is associated with the snow removed, so SNOW,
SNOWLQ, and CC are all reduced by 1 - FRAC. If PSNVP is negative, condensation is occurring;
SNOW, SNOWLQ, and CC are correspondingly increased. This is simple, but not quite accurate.
If no snow is left the routine ends.

The amount of snowpack warming or cooling is calculated next. The equivalent amount of ice
melted by the energy input from SNOEN and the heat included in any warm rain is

EQEN = DTP * (SNOEN + RTHR * RMAX(TA, 0) * CVLQ) / LF

where CVLQ is the specific heat of water, and LF is the latent heat of fusion of water. Both
CVLQ and LF are constant. When EQEN is less than 0, the snow is cooling; first any SNOWLQ is
refrozen, then CC is increased. CC is not allowed to be reduced so that TSNOW is below the
mean daily air temperature, TA, although it may remain colder, i.e. if TA < 0°, CC can only
be reduced to TSNOW = TA. When EQEN is greater than 0, the snow is warming; first CC is
reduced, then SNOWLQ is produced, and finally melt (SMLT) occurs.

Finally, any rain throughfall (RTHR) is added to the snowpack. If any CC exists it refreezes
rain until the CC is "used up". Any additional rain then increases SNOWLQ; when the maximum
SNOWLQ is reached, the input of rain to the snow (RSNO) has also reached its maximum. In all
cases the final results are a new SNOW, new SNOLQ and CC, and a value of RSNO.

MSBPREINT then calculates the rain passing through the snow (RNET) as

RNET = RTHR - RSNO.

When SNOW exists at the beginning of the day, soil evaporation (SLVP) is zero.
"
"""
function SNOWPACK(p_fu_RTHR, p_fu_STHR, p_fu_PSNVP, p_fu_SNOEN, u_CC, u_SNOW,
                  u_SNOWLQ, p_DTP, p_fu_TA, p_MAXLQF, p_GRDMLT,
                  p_CVICE, p_LF, p_CVLQ)

    ####
    # Operator splitting: step 1: update u_SNOW and u_CC using snow throughfall:
    # Update snow amount and its cold content, while u_SNOWLQ remains unchanged
    du_SNOW_step1 = p_fu_STHR
    du_CC_step1   = p_CVICE * max(-p_fu_TA, 0.) * p_fu_STHR
    # Update u_SNOW, u_CC
    u_SNOW = u_SNOW + du_SNOW_step1 * p_DTP
    u_CC   = u_CC   + du_CC_step1   * p_DTP

    ####
    # Operator splitting: step 2: update u_SNOWLQ and u_CC using cold content:
    # Variant 1:
    if (u_CC > 0. && u_SNOWLQ > 0.)
        if (u_CC > u_SNOWLQ * p_LF)
            # refreeze all liquid
            u_CC = u_CC - u_SNOWLQ * p_LF
            u_SNOWLQ = 0.
        else
            # refreeze part
            u_SNOWLQ = u_SNOWLQ - u_CC / p_LF
            u_CC = 0.
        end
    end
    # Variant 2: deleted

    ####
    # Operator splitting: step 3: update u_SNOW, u_SNOWLQ, and u_CC as well as
    #                             define aux_du_SMLT and aux_du_SNVP using p_GRDMLT and PSNVP
    # GRDMLT - rate of groundmelt of snowpack (mm/d)
    # PSNVP  - potential snow evaporation (mm/d)

    # groundmelt and evaporation loss as fraction of u_SNOW
    # variant1:
    FRAC = (p_GRDMLT + p_fu_PSNVP) * p_DTP / u_SNOW
    # FRAC can be negative if condensation (p_fu_PSNVP<0) exceeds groundmelt
    if (FRAC < 1)
        # entire potential_SNOW_decrease can be realized
        aux_du_SMLT = p_GRDMLT
        aux_du_SNVP = p_fu_PSNVP
        # reduce u_CC, u_SNOWLQ, and u_SNOW proportionally for groundmelt and evaporation
        # increase them proportionally if condensation exceeds groundmelt
        u_SNOW   = u_SNOW   * (1. - FRAC)
        u_SNOWLQ = u_SNOWLQ * (1. - FRAC)
        u_CC     = u_CC     * (1. - FRAC)
    else
        # all u_SNOW disappears from groundmelt and/or evaporation
        # potential rates need to be adapted to actual ones
        aux_du_SMLT = p_GRDMLT / FRAC
        aux_du_SNVP = p_fu_PSNVP / FRAC
        aux_du_RSNO = 0. # FB: unused
        u_SNOW      = 0.
        u_SNOWLQ    = 0.
        u_CC        = 0.
    end
    # variant2: deleted


    ####
    # Operator splitting: step 4: update u_SNOWLQ, and u_CC by considering
    #                             snowpack cooling or warming
    if (u_SNOW > 0)
        # equivalent ice melted by energy input including warm rain (mm)
        EQEN = p_DTP * (p_fu_SNOEN + p_fu_RTHR * max(p_fu_TA, 0) * p_CVLQ) / p_LF
        if (EQEN <= 0)
            # case 1) snowpack cooling
            NMLT = -EQEN
            if (NMLT < u_SNOWLQ)
                # 1a) only part of u_SNOWLQ refreezes
                u_CC = 0.
                # should be 0 already because u_SNOWLQ is positive
                u_SNOWLQ = u_SNOWLQ - NMLT
            else
                # 1b) all u_SNOWLQ (if any) refreezes, remaining NMLT increases u_CC
                NMLT = NMLT - u_SNOWLQ
                u_SNOWLQ = 0.
                u_CC = u_CC + NMLT * p_LF
                # do not allow p_fu_TSNOW to cool below p_fu_TA
                u_CC = min(u_CC, -p_fu_TA * u_SNOW * p_CVICE)
            end
        else
            # case 2) snowpack warming  (cant have both u_CC and u_SNOWLQ)
            if (EQEN * p_LF < u_CC || p_fu_TA < 0)
                # 2a) reduce but dont eliminate u_CC
                if (p_fu_TA < 0)
                    # do not allow p_fu_TSNOW to warm above p_fu_TA when p_fu_TA < 0
                    u_CC = max(u_CC - EQEN * p_LF, -p_fu_TA * u_SNOW * p_CVICE)
                    u_SNOWLQ = 0.
                else
                    u_CC = u_CC - EQEN * p_LF
                    u_SNOWLQ = 0.
                end
            else
                # 2b) u_CC eliminated
                EQEN = EQEN - u_CC / p_LF
                if (EQEN <= p_MAXLQF * u_SNOW - u_SNOWLQ)
                    # 2b-i) remaining energy increases liquid water
                    u_SNOWLQ = u_SNOWLQ + EQEN
                    u_CC = 0.
                    # aux_du_SMLT and u_SNOW unchanged
                else
                    # 2b-ii) liquid water capacity reached, u_SNOW melt produced
                    EQEN = EQEN - (p_MAXLQF * u_SNOW - u_SNOWLQ)
                    if (u_SNOW * (1. - p_MAXLQF) > EQEN)
                        # 2b-ii-A) melt is ice plus the liquid included in it
                        aux_du_SMLT = aux_du_SMLT + (EQEN / p_DTP) / (1. - p_MAXLQF)
                        u_SNOW = u_SNOW - EQEN / (1. - p_MAXLQF)
                        u_SNOWLQ = p_MAXLQF * u_SNOW
                        u_CC = 0.
                    else
                        # 2b-ii-B) all u_SNOW melts
                        aux_du_SMLT = aux_du_SMLT + u_SNOW / p_DTP
                        u_SNOW = 0.
                        u_SNOWLQ = 0.
                        u_CC = 0.
                    end
                end
            end
        end

        # add rain to snowpack,
        if (p_fu_RTHR == 0 || u_SNOW == 0)
            # FB: note u_SNOW will never == 0 at this point, as this is inside: "if (u_SNOW > 0)"
            aux_du_RSNO = 0.
        else
            # rain on u_SNOW
            RIN = p_fu_RTHR * p_DTP
            if (u_CC > 0)
                # use u_CC to refreeze rain
                if (u_CC > RIN * p_LF)
                    # refreezes all rain
                    u_CC = u_CC - RIN * p_LF
                    aux_du_RSNO = p_fu_RTHR
                    u_SNOW = u_SNOW + RIN
                else
                    # u_CC refreezes part of rain
                    u_SNOW = u_SNOW + u_CC / p_LF
                    aux_du_RSNO = (u_CC / p_LF) / p_DTP
                    u_CC = 0.
                    # remaining rain
                    RIN = RIN - aux_du_RSNO * p_DTP
                    # increase liquid water, u_SNOWLQ initially zero
                    if (RIN < p_MAXLQF * u_SNOW / (1. - p_MAXLQF))
                        # remaining RIN all to u_SNOWLQ
                        u_SNOWLQ = RIN
                        aux_du_RSNO = aux_du_RSNO + RIN / p_DTP
                        u_SNOW = u_SNOW + RIN
                    else
                        u_SNOWLQ = p_MAXLQF * u_SNOW / (1. - p_MAXLQF)
                        aux_du_RSNO = aux_du_RSNO + u_SNOWLQ / p_DTP
                        u_SNOW = u_SNOW + u_SNOWLQ
                    end
                end
            else
                # u_CC = 0..
                if (u_SNOWLQ >= p_MAXLQF * u_SNOW)
                    # u_SNOW already holding maximum liquid
                    aux_du_RSNO = 0.
                else
                    ALQ = p_MAXLQF * u_SNOW - u_SNOWLQ
                    if (RIN < ALQ)
                        # all RIN to u_SNOW
                        aux_du_RSNO = p_fu_RTHR
                        u_SNOWLQ = u_SNOWLQ + RIN
                        u_SNOW = u_SNOW + RIN
                    else
                        # maximum liquid reached
                        aux_du_RSNO = (ALQ / (1. - p_MAXLQF)) / p_DTP
                        u_SNOW = u_SNOW + aux_du_RSNO * p_DTP
                        u_SNOWLQ = p_MAXLQF * u_SNOW
                    end
                end
            end
        end
    end
    return (u_CC, u_SNOW, u_SNOWLQ, aux_du_RSNO, aux_du_SNVP, aux_du_SMLT)
end



end
