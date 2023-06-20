# TODO(bernhard): think about where to put function definition of MSBSETVARS(), MSBDAYNIGHT()
"""
    MSBSETVARS()

Compute state dependent parameters for updating states INTS, INTR, SNOW, CC, SNOWLQ in
callback function.

A) get: 1) sunshine durations, 2) SFAL, 3) plant resistance (TODO: could be done before simulation) 
B) get_windspeed_from_canopy_and_snowpack: (p_fu_UADTM, p_fu_UANTM) = f(...) 
C) get_potential_snowpack_and_soil_evaporation_from_canopy,snowpack,atmosphere,: (p_fu_UADTM, p_fu_UANTM) = f(...)

# Arguments
- many

" From ecoshift: 
Subroutine MSUBSETVARS contains subroutines that calulate derived variables for the day.
SUNDS, CANOPY, ROUGH, and PLNTRES are called to get solar, canopy structure, roughness, and
plant resistance variables that depend on day of the year. Subsurface heat flux (SHEAT) is
always set to zero. Subroutine WEATHER estimates missing values, modifies input data as
necessary, corrects weather station wind speed to wind speed above the canopy, and
determines daytime and nighttime temperatures and wind speeds. Subroutine SNOFRAC determines
the fraction of daily precipitation that is snow (SNOFRC). If there is no snow on the
ground, the soil evaporation resistance (RSS) is obtained as function FRSS. When there is
snow on the ground, the snowpack temperature (TSNOW) is calculated from cold content (CC).
Subroutine SNOVAP then estimates the snow evaporation rate. Subroutine SNOENRGY obtains the
energy available for snowmelt (SNOEN) from mean daily temperature. The factor is modified
for canopy cover as determined by LAI and SAI. Snow evaporation or condensation depends on
the aerodynamic resistances and the vapor gradient; however, an arbitrary reduction factor
is required. "
"""
function MSBSETVARS(# arguments
                    FLAG_MualVanGen, NLAYER, p_soil,
                    # for SUNDS
                    p_LAT, p_ESLOPE, DOY, p_L1, p_L2,
                    # for CANOPY
                    p_fT_HEIGHT, p_fT_LAI, p_fT_SAI, u_SNOW, p_SNODEN, p_MXRTLN, p_MXKPL, p_fT_DENSEF,
                    #
                    p_Z0S, p_Z0G,
                    # for ROUGH
                    p_ZMINH, p_CZS, p_CZR, p_HS, p_HR, p_LPC, #p_CS,
                    # for PLNTRES
                    p_fT_RELDEN, p_RTRAD, p_FXYLEM,
                    # for WEATHER
                    p_fT_TMAX, p_fT_TMIN, p_fT_VAPPRES, p_fT_UW, p_WNDRAT, p_FETCH, p_Z0W, p_ZW, p_fT_SOLRAD,
                    # for SNOFRAC
                    p_RSTEMP,
                    #
                    u_CC, p_CVICE,
                    # for SNOVAP
                    p_LWIDTH, p_RHOTP, p_NN, p_KSNVP,
                    #
                    p_ALBSN, p_ALB,
                    # for FRSS
                    p_RSSA, p_RSSB, u_aux_PSIM, #u_aux_PSIM[1]
                    # for SNOENRGY
                    p_CCFAC, p_MELFAC, p_LAIMLT, p_SAIMLT)

    # start A: get: 1) sunshine durations, 2) SFAL, 3) plant resistance
    # TODO(bernhard): a) Do this outside of integration loop in define_LWFB90_p() p_fT_DAYLEN
    # solar parameters depending on only on DOY
    p_fT_DAYLEN, p_fT_I0HDAY, p_fT_SLFDAY = LWFBrook90.SUN.SUNDS(p_LAT, p_ESLOPE, DOY, p_L1, p_L2, LWFBrook90.CONSTANTS.p_SC, LWFBrook90.CONSTANTS.p_PI, LWFBrook90.CONSTANTS.p_WTOMJ)
    # canopy evolution depending only on DOY
    p_fT_LAIeff, p_fT_SAIeff, p_fT_RTLENeff, p_fT_RPLANT = LWFBrook90.PET.CANOPY_timeEvolution(
        p_fT_LAI, # leaf area index, m2/m2, minimum of 0.00001
        p_fT_SAI,  # stem area index, m2/m2
        p_MXRTLN,  # maximum root length per unit land area, m/m2
        p_MXKPL,   # maximum plant conductivity, (mm/d)/MPa
        #p_CS,     # ratio of projected SAI to canopy height, m-1, not needed in this version
        p_fT_DENSEF)  # density factor
    # fraction of precipitation as p_fT_SFAL
    p_fT_SNOFRC= LWFBrook90.SNO.SNOFRAC(p_fT_TMAX, p_fT_TMIN, p_RSTEMP)
    # plant resistance components
    p_fT_RXYLEM, p_fT_RROOTI, p_fT_ALPHA = LWFBrook90.EVP.PLNTRES(NLAYER, p_soil, p_fT_RTLENeff, p_fT_RELDEN, p_RTRAD, p_fT_RPLANT, p_FXYLEM, LWFBrook90.CONSTANTS.p_PI, LWFBrook90.CONSTANTS.p_RHOWG)
    # end A: get: 1) sunshine durations, 2) SFAL, 3) plant resistance

    # start B: get_windspeed_from_canopy_and_snowpack: (p_fu_UADTM, p_fu_UANTM) = f(...)
    # Compute influence of snow on canopy height, LAI, and roughness and how this influences wind speed
    # Modify canopy parameters (depending on DOY -> fT) due to dependency on snow depth (state parameter -> fu)
    p_fu_HEIGHTeff, p_fu_LAIeff = LWFBrook90.PET.CANOPY_snowCover( # height and LAI reduced by snowcover
        p_fT_HEIGHT,
        p_fT_LAIeff,  # leaf area index, m2/m2, not yet reduced by snowcover
        u_SNOW,       # water equivalent of snow on the ground, mm SWE
        p_SNODEN)     # snow density, mm SWE/mm depth
    # roughness parameters depending on u_SNOW
    p_fu_Z0GS, p_fu_Z0C, p_fu_DISPC, p_fu_Z0, p_fu_DISP, p_fu_ZA =
            LWFBrook90.PET.ROUGH(p_fu_HEIGHTeff, p_ZMINH, p_fu_LAIeff, p_fT_SAIeff,
                                 p_CZS, p_CZR, p_HS, p_HR, p_LPC, ifelse(u_SNOW > 0, p_Z0S, p_Z0G))
    # calculated weather data
    p_fu_SHEAT = 0.
    (p_fT_SOLRADC, p_fT_TA, p_fT_TADTM, p_fT_TANTM, UA, p_fu_UADTM, p_fu_UANTM) =
        LWFBrook90.PET.WEATHER(p_fT_TMAX, p_fT_TMIN, p_fT_DAYLEN, p_fT_I0HDAY, p_fT_VAPPRES, p_fT_UW, p_fu_ZA, p_fu_DISP, p_fu_Z0, p_WNDRAT, p_FETCH, p_Z0W, p_ZW, p_fT_SOLRAD)
    # end B: get_windspeed_from_canopy_and_snowpack: (p_fu_UADTM, p_fu_UANTM) = f(...)

    # start C: get_potential_snowpack_and_soil_evaporation_from_canopy,snowpack,atmosphere,: (p_fu_UADTM, p_fu_UANTM) = f(...)
    if (u_SNOW > 0)
        # snowpack temperature at beginning of day
        p_fu_TSNOW = -u_CC / (p_CVICE * u_SNOW)
        # potential snow evaporation PSNVP
        p_fu_PSNVP=LWFBrook90.SNO.SNOVAP(p_fu_TSNOW, p_fT_TA, p_fT_VAPPRES, UA, p_fu_ZA, p_fu_HEIGHTeff, p_fu_Z0, p_fu_DISP, p_fu_Z0C, p_fu_DISPC, p_fu_Z0GS, p_LWIDTH, p_RHOTP, p_NN, p_fu_LAIeff, p_fT_SAIeff, p_KSNVP)
        p_fu_ALBEDO = p_ALBSN
    else
        p_fu_TSNOW = 0.
        p_fu_PSNVP = 0.
        p_fu_ALBEDO = p_ALB
    end
    # snow surface energy balance (is performed even when SNOW=0 in case snow is added during day)
    p_fu_SNOEN = LWFBrook90.SNO.SNOENRGY(p_fu_TSNOW, p_fT_TA, p_fT_DAYLEN, p_CCFAC, p_MELFAC, p_fT_SLFDAY, p_fu_LAIeff, p_fT_SAIeff, p_LAIMLT, p_SAIMLT)
    # soil evaporation resistance
    if (u_SNOW > 0)
        p_fu_RSS = 0.
    else
        p_fu_RSS = LWFBrook90.PET.FRSS(p_RSSA, p_RSSB, u_aux_PSIM[1], p_soil) # TODO: with smaller discretizations, maybe use more than one layer? (e.g. damping depth see p.115/119 or Table9.2 in Campbell-1998- Env.Biophys.)
        # check for zero or negative p_fu_RSS (TODO: not done in LWFBrook90)
        if (p_fu_RSS < 0.000001)
           error("""
            p_fu_RSS is very small or negative ($p_fu_RSS). Run ends.
            Happended at PSIM: $(u_aux_PSIM[1])kPa.
            Check p_RSSA ($p_RSSA) and p_RSSB ($p_RSSB) values.""")
        end
    end
    # end C: get_potential_snowpack_and_soil_evaporation_from_canopy,snowpack,atmosphere,: (p_fu_UADTM, p_fu_UANTM) = f(...)

    return (p_fT_DAYLEN, p_fT_I0HDAY, p_fT_SLFDAY,
            p_fu_HEIGHTeff, p_fu_LAIeff, p_fT_SAIeff,
            p_fu_Z0GS, p_fu_Z0C, p_fu_DISPC, p_fu_Z0, p_fu_DISP, p_fu_ZA,
            p_fT_RXYLEM, p_fT_RROOTI, p_fT_ALPHA,
            p_fu_SHEAT,
            p_fT_SOLRADC, p_fT_TA, p_fT_TADTM, p_fT_TANTM, p_fu_UADTM, p_fu_UANTM,
            p_fT_SNOFRC,
            p_fu_TSNOW,p_fu_PSNVP, p_fu_ALBEDO,p_fu_RSS,
            p_fu_SNOEN)
end

"""
    MSBDAYNIGHT()

Calculate average daily rate of potential and actual interception, evaporation, and
transpiration by considering weighted average of rate during day and rate during night:
Subroutine MSBDAYNIGHT - day-night loop: Compute day and night rates

# Arguments
- many

http://www.ecoshift.net/brook/b90.html:
"
Subroutine MSBDAYNIGHT contains routines that calculate the five components of evaporation
(see Flow Chart):

- evaporation of intercepted rain (IRVP) 
- evaporation of intercepted snow (ISVP) 
- evaporation from snow (SNVP) 
- soil evaporation (SLVP) from the top soil layer 
- transpiration (TRANI) from each soil layer that contains roots 

These evaporation values are obtained separately for daytime and nightime, then combined 
into a daily values.
Note that evaporation of intercepted storage and snow evaporation (IRVP, ISVP, SNVP) are 
reduced if their sources disappear. (This happens later in MSBPREINT.)

Potential evaporation rates are obtained using the Shuttleworth and Wallace (1985)
modification of the Penman-Monteith approach. Daily solar radiation is corrected for slope,
is allocated to the daytime, and is converted to average daytime rate. Subroutine AVAILEN
calculates available energy (net radiation minus SHEAT=0) at the top (AA) and at the bottom
(ASUBS) of the canopy, using a Beers Law extinction coefficient. The three aerodynamic
resistances (RAA, RAC, RAS) needed by the Shuttleworth-Wallace method are obtained in
subroutine SWGRA, using algorithms of Shuttleworth and Gurney (1990). These resistances
depend on leaf area index (LAI), which can vary seasonally, and on canopy height, which
determines stem area index (SAI). The canopy surface resistance to transpiration (RSC) for
the daytime is obtained in subroutine SRSC; it depends on maximum leaf conductance, reduced
for humidity, temperature, and light penetration. At night RSC is the reciprocal of leaf
area index (LAI) times minimum leaf conductance (GLMIN). Soil evaporation resistance (RSS)
depends on soil water potential in the top soil layer. Subroutine SWPE uses AA, ASUBS, RSC,
RSS, RAA, RAC, RAS and the vapor pressure deficit (VPD) to calculate potential transpiration
(PTR) and the associated ground or soil evaporation (GER) as given by Shuttleworth and
Wallace (1985). Subroutine SWPE is then called again with RSC = 0 to give the intercepted
evaporation rate and its associated soil evaporation (PIR and GIR). Subroutine TBYLAYER
obtains actual transpiration by layer (TRANI). Actual transpiration is the lesser of
potential transpiration and a soil water supply rate determined by the resistance to liquid
water flow in the plants and on root distribution and soil water potential in the soil
layers. If the actual transpiration is less than the potential, a new, higher GER is
calculated by subroutine SWGE. After the MSBDAYNIGHT day-night loop, these evaporation rates
are weighted for daytime and nighttime according to daylength (DAYLEN), and the daily
average rates are then used in later calculations.

The "precipitation interval" is equal to one day when daily precipitation is provided along
with other daily weather data in the Data File; parameter NPINT is then set to 1 and the
precipitation loop is passed through once a day. Alternatively, precipitation data at fixed
intervals of less than a day can be read from a Precip Interval File. The Precip Interval
File can have one line per day, allowing easy reruns with different daily precipitation.
Then NPINT is set to the number of precipitation intervals per day; the precipitation loop
is passed through NPINT times per day, and a line of the Precip Interval File is read each
time. If available, measured flow for the interval can also be input.
"
http://www.ecoshift.net/brook/pet.html
"
BROOK90 obtains evaporation rates separately for daytime and nighttime within a day-night
evaporation loop. All solar radiation (SOLRAD) is assigned to the daytime. The atmospheric
humidity (EA) is assumed constant through the day ("day" refers to 24 hours). The daytime
and nighttime values of air temperature and wind speed are obtained in subroutine WEATHER
using function WNDADJ. Vapor pressure deficit (VPD) is obtained using subroutine ESAT.
Subroutine CANOPY uses function INTERP to obtain canopy structure variables for the day.
Subroutine ROUGH gets canopy roughness parameters. Within a day-night loop, the three
aerodynamic resistances needed by the Shuttleworth-Wallace method are calculated in
subroutine SWGRA. The canopy surface resistance (RSC) for the daytime is obtained from
subroutine SRSC, and the soil surface resistance (RSS) in function FRSS. Subroutine SWPE
uses function PM along with the various resistances to obtain potential transpiration rate
(PTR) and the associated ground or soil evaporation rate (GER) by the Shuttleworth-Wallace
equations. Subroutine SWPE is called again with RSC = 0 to give potential interception rate
(PIR) and its associated soil evaporation rate (GIR). Subroutine TBYLAYER obtains actual
transpiration by layer (ATRANI). If the actual transpiration is less than the potential, a
new, higher GER is calculated by subroutine SWGE. BROOK90 then weights the daytime and
nighttime rates by the solar daylength (DAYLEN) to obtain average rates for the day, PTRAN,
GEVP, PINT, GIVP, and TRANI, which are used in later calculations.
"
"""
function MSBDAYNIGHT(p_fT_SLFDAY, p_fT_SOLRADC, p_WTOMJ, p_fT_DAYLEN, p_fT_TADTM, p_fu_UADTM, p_fT_TANTM, p_fu_UANTM,
                     p_fT_I0HDAY,
                     # for AVAILEN:
                     p_fu_ALBEDO, p_C1, p_C2, p_C3, p_fT_VAPPRES, p_fu_SHEAT, p_CR, p_fu_LAIeff, p_fT_SAIeff,
                     # for SWGRA:
                     p_fu_ZA, p_fu_HEIGHTeff, p_fu_Z0, p_fu_DISP, p_fu_Z0C, p_fu_DISPC, p_fu_Z0GS, p_LWIDTH, p_RHOTP, p_NN,
                     # for SRSC:
                     p_fT_TA, p_GLMIN, p_GLMAX, p_R5, p_CVPD, p_RM, p_TL, p_T1, p_T2, p_TH,
                     # for SWPE:
                     p_fu_RSS,
                     # for TBYLAYER:
                     p_fT_ALPHA, p_fu_KK, p_fT_RROOTI, p_fT_RXYLEM, u_aux_PSITI, NLAYER, p_PSICR, NOOUTF)
    # MSBDAYNIGHT() computes the five components of evaporation:
    # - aux_du_ISVP: evaporation of intercepted snow
    # - aux_du_IRVP: evaporation of intercepted rain
    # - aux_du_SNVP: evaporation from snow
    # - aux_du_SLVP: soil evaporation from the top soil layer
    # - aux_du_TRANI: transpiration from each soil layer that contains roots

    p_fu_PTR = fill(NaN, 2)
    p_fu_GER = fill(NaN, 2)
    p_fu_PIR = fill(NaN, 2)
    p_fu_GIR = fill(NaN, 2)
    p_fu_ATRI=fill(NaN,2,NLAYER)

    ATR = fill(NaN, 2)
    SLRAD=fill(NaN,2)
    for J = 1:2 # 1 for daytime, 2 for nighttime

        # net radiation
        if (J ==1)
            SLRAD[J] = p_fT_SLFDAY * p_fT_SOLRADC / (p_WTOMJ * p_fT_DAYLEN)
            TAJ = p_fT_TADTM
            UAJ = p_fu_UADTM
        else
            SLRAD[J] = 0
            TAJ = p_fT_TANTM
            UAJ = p_fu_UANTM
        end

        # if (p_fT_I0HDAY <= 0.01)
        #     # TODO(bernhard): Brook90 did treat this case specially, LWFBrook90 did not
        #
        #     # no sunrise, assume 50% clouds for longwave
        #     cloud_fraction = 0.5
        # else
        cloud_fraction = p_fT_SOLRADC / p_fT_I0HDAY
        # end
        # Available energy above canopy (AA) and at soil/substrate (ASUBS)
        AA, ASUBS =
            LWFBrook90.SUN.AVAILEN(SLRAD[J], p_fu_ALBEDO, p_C1, p_C2, p_C3, TAJ, p_fT_VAPPRES,
                    cloud_fraction,
                    p_fu_SHEAT, p_CR, p_fu_LAIeff, p_fT_SAIeff)

        # vapor pressure deficit
        ES, DELTA = LWFBrook90.PET.ESAT(TAJ)
        VPD = ES - p_fT_VAPPRES
        # S.-W. resistances
        RAA, RAC, RAS = LWFBrook90.PET.SWGRA(UAJ, p_fu_ZA, p_fu_HEIGHTeff, p_fu_Z0, p_fu_DISP, p_fu_Z0C, p_fu_DISPC, p_fu_Z0GS, p_LWIDTH, p_RHOTP, p_NN, p_fu_LAIeff, p_fT_SAIeff)
        if (J == 1)
            RSC=LWFBrook90.PET.SRSC(SLRAD[J], p_fT_TA, VPD, p_fu_LAIeff, p_fT_SAIeff, p_GLMIN, p_GLMAX, p_R5, p_CVPD, p_RM, p_CR, p_TL, p_T1, p_T2, p_TH)
        else
            RSC = 1 / (p_GLMIN * p_fu_LAIeff)
        end
        # RSC: canopy surface resistance

        #print("\nIDAY:$(@sprintf("% 3d", IDAY)), J = $J     AA:$(@sprintf("% 8.4f", AA)), ASUBS:$(@sprintf("% 8.4f", ASUBS)), VPD:$(@sprintf("% 8.4f", VPD)), RAA:$(@sprintf("% 8.4f", RAA)), RAC:$(@sprintf("% 8.4f", RAC)), RAS:$(@sprintf("% 8.4f", RAS)), RSC:$(@sprintf("% 8.4f", RSC)), p_fu_RSS:$(@sprintf("% 8.4f", p_fu_RSS)), DELTA:$(@sprintf("% 8.4f", DELTA))")

        # S.-W. potential transpiration and ground evaporation rates
        p_fu_PTR[J], p_fu_GER[J] =  LWFBrook90.PET.SWPE(AA, ASUBS, VPD, RAA, RAC, RAS, RSC, p_fu_RSS, DELTA)
        # S.-W. potential interception evaporation and ground evap. rates
        # RSC = 0, p_fu_RSS not changed
        p_fu_PIR[J], p_fu_GIR[J] =  LWFBrook90.PET.SWPE(AA, ASUBS, VPD, RAA, RAC, RAS, 0,   p_fu_RSS, DELTA)

        # actual transpiration and ground evaporation rates
        if (p_fu_PTR[J] > 0.001)
            ATR[J], ATRANI = LWFBrook90.EVP.TBYLAYER(J, p_fu_PTR[J], p_fu_DISPC, p_fT_ALPHA, p_fu_KK, p_fT_RROOTI, p_fT_RXYLEM, u_aux_PSITI, NLAYER, p_PSICR, NOOUTF)
            for i = 1:NLAYER
                p_fu_ATRI[J,i] = ATRANI[i]
            end
            if (ATR[J] < p_fu_PTR[J])
                # soil water limitation, new GER
                p_fu_GER[J]=LWFBrook90.PET.SWGE(AA, ASUBS, VPD, RAA, RAS, p_fu_RSS, DELTA, ATR[J])
            end
        else
            # no transpiration, condensation ignored
            p_fu_PTR[J] = 0
            ATR[J] = 0
            for i = 1:NLAYER
                p_fu_ATRI[J,i] = 0
            end
            p_fu_GER[J]=LWFBrook90.PET.SWGE(AA, ASUBS, VPD, RAA, RAS, p_fu_RSS, DELTA, 0)
        end
    end
    #print("\nIDAY:$(@sprintf("% 3d", IDAY)), p_fu_GER[1]: $(@sprintf("% 8.4f",p_fu_GER[1])), p_fu_GIR[1]: $(@sprintf("% 8.4f",p_fu_GIR[1]))")
    # print(", p_fu_GER[2]: $(@sprintf("% 8.4f",p_fu_GER[2])), p_fu_GIR[2]: $(@sprintf("% 8.4f",p_fu_GIR[2]))")

    # print(", AA:$(@sprintf("% 8.4f", AA)), ASUBS:$(@sprintf("% 8.4f", ASUBS)), VPD:$(@sprintf("% 8.4f", VPD)), RAA:$(@sprintf("% 8.4f", RAA)), RAC:$(@sprintf("% 8.4f", RAC)), RAS:$(@sprintf("% 8.4f", RAS)), p_fu_RSS:$(@sprintf("% 8.4f", p_fu_RSS)), DELTA:$(@sprintf("% 8.4f", DELTA))")
    # print(", J:$(@sprintf("% 8.4f", J)), p_fu_PIR[J]:$(@sprintf("% 8.4f", p_fu_PIR[J])))")
    # print("\n        ")
    return (p_fu_PTR, # potential transpiration rate for daytime or night (mm/d)
            p_fu_GER, # ground evaporation rate for daytime or night (mm/d)
            p_fu_PIR, # potential interception rate for daytime or night (mm/d)
            p_fu_GIR, # ground evap. rate with intercep. for daytime or night (mm/d)
            p_fu_ATRI)# actual transp.rate from layer for daytime and night (mm/d)

    # return (#SLRAD[2],SLRAD[1],TAJ,UAJ, SOVERI, AA, ASUBS
    #         #ES, DELTA, VPD, RAA, RAC, RAS, RSC,
    #         p_fu_PTR, p_fu_GER, p_fu_PIR, p_fu_GIR,
    #         #ATR, ATRANI,
    #         p_fu_ATRI)
end

"""
    MSBDAYNIGHT_postprocess()

Calculate average daily rate of potential and actual interception, evaporation, and
transpiration by considering weighted average of rate during day and rate during night:
Subroutine MSBDAYNIGHT_postprocess - Combine day and night rates to average daily rate

# Arguments
- many

http://www.ecoshift.net/brook/pet.html
"""
function MSBDAYNIGHT_postprocess(NLAYER,
                                 p_fu_PTR, # potential transpiration rate for daytime or night (mm/d)
                                 p_fu_GER, #      ground evaporation rate for daytime or night (mm/d)
                                 p_fu_PIR, # potential   evaporation rate of   interception storage for daytime or night (mm/d)
                                 p_fu_GIR, # ground evaporation      rate with interception storage for daytime or night (mm/d)
                                 p_fu_ATRI,# actual transp.rate from layer for daytime and night (mm/d)
                                 p_fT_DAYLEN)

    # average rates over day (mm/day)
    # mm/day   =  mm/day      * day_fraction+ mm/day      *  day_fraction
    p_fu_PTRAN = (p_fu_PTR[1] * p_fT_DAYLEN + p_fu_PTR[2] * (1 - p_fT_DAYLEN))
    p_fu_GEVP  = (p_fu_GER[1] * p_fT_DAYLEN + p_fu_GER[2] * (1 - p_fT_DAYLEN))
    p_fu_PINT  = (p_fu_PIR[1] * p_fT_DAYLEN + p_fu_PIR[2] * (1 - p_fT_DAYLEN))
    p_fu_GIVP  = (p_fu_GIR[1] * p_fT_DAYLEN + p_fu_GIR[2] * (1 - p_fT_DAYLEN))

    aux_du_TRANI=zeros(NLAYER)

    for i = 1:NLAYER
        aux_du_TRANI[i] = (p_fu_ATRI[1, i] * p_fT_DAYLEN + p_fu_ATRI[2, i] * (1 - p_fT_DAYLEN))
    end

    return (p_fu_PTRAN, # average potential transpiration rate for day (mm/d)
            p_fu_GEVP,  # average      ground evaporation rate for day (mm/d)
            p_fu_PINT,  # average potential   evaporation rate of   interception storage for day (mm/d)
            p_fu_GIVP,  # average ground evaporation      rate with interception storage for day (mm/d)
            aux_du_TRANI) # average transpiration rate for day from layer (mm/d)
end

"""
    MSBPREINT()

http://www.ecoshift.net/brook/b90.html
In subroutine MSBPREINT, precipitation is separated into rain and snow using SNOFRC. If
NPINT = 1, subroutine INTER24 is called twice, once for snow interception and once for rain
interception; this routine uses the monthly parameter DURATN, which is the average storm
duration in hours. If NPINT > 1, subroutine INTER is used instead, and the precipitation is
assumed to occur over the whole precipitation interval. Transpiration (TRAN) for the
interval is summed over layers and reduced by the fraction of time the canopy is wet; soil
evaporation (SLVP) is GIR when the canopy is wet and GER when it is dry. If a snowpack
exists, subroutine SNOWPACK is called to use SNOEN and the rain and snow throughfall to
calculate snowmelt (SMLT), cold content (CC), and liquid water content of the snow. Net rain
to the soil surface (RNET) is rain throughfall (RTHR) minus rain absorbed by the snowpack
(RSNO). Water reaching the ground surface is RNET + SMLT.
"""
function MSBPREINT(#arguments:
                   p_fT_PREC, p_DTP, p_fT_SNOFRC, p_NPINT, p_fu_PINT, p_fT_TA,
                   # for INTER (snow)
                   u_INTS, p_fu_LAI, p_fu_SAI, p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS,
                   # for INTER (rain)
                   u_INTR, p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS,
                   # for INTER24 (snow + rain)
                   p_DURATN, MONTHN,
                   #
                   u_SNOW, p_fu_PTRAN, NLAYER, aux_du_TRANI, p_fu_GIVP, p_fu_GEVP,
                   # for SNOWPACK
                   u_CC, u_SNOWLQ, p_fu_PSNVP, p_fu_SNOEN, p_MAXLQF, p_GRDMLT,
                   p_CVICE, p_LF, p_CVLQ)

    p_fT_SFAL = p_fT_SNOFRC * p_fT_PREC # rate in mm/day
    p_fT_RFAL = p_fT_PREC - p_fT_SFAL   # rate in mm/day

    # A) Interception by canopy ################
    if (p_NPINT > 1.0)
        # more than one precip interval in day
        error("Case with multiple precipitation intervals (using PRECDAT and precip_interval != 1) is not implemented.")
        # # snow interception
        # if (p_fu_PINT < 0 && p_fT_TA > 0)
        #     # prevent frost when too warm, carry negative p_fu_PINT to rain  # NOTE: this case is not done in LWFBRook90R
        #     aux_du_SINT, aux_du_ISVP = LWFBrook90.EVP.INTER(p_fT_SFAL, 0, p_fu_LAI, p_fu_SAI, p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS, p_DTP, u_INTS)
        # else
        #     aux_du_SINT, aux_du_ISVP = LWFBrook90.EVP.INTER(p_fT_SFAL, p_fu_PINT, p_fu_LAI, p_fu_SAI, p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS, p_DTP, u_INTS)
        # end
        # # rain interception,  note potential interception rate is p_fu_PINT-aux_du_ISVP (mm/day)
        # aux_du_RINT, aux_du_IRVP = LWFBrook90.EVP.INTER(p_fT_RFAL, p_fu_PINT - aux_du_ISVP, p_fu_LAI, p_fu_SAI, p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS, p_DTP, u_INTR)
    else
        # one precip interval in day, use storm p_DURATN and INTER24
        # a) snow interception storage (adding SFAL and potential condensation/evaporation p_fu_PINT)
        if (p_fu_PINT < 0 && p_fT_TA > 0)
            # if potential evaporation rate of interception storage is negative (i.e. condensation),
            # atmospheric humidity will be deposited in the interception storage. However, when too warm TA > 0, deposit
            # a negative p_fu_PINT into the rain interception storage.

            # prevent frost when too warm, carry negative p_fu_PINT to rain interception storage # NOTE: this case is not done in LWFBRook90R
            aux_du_SINT, aux_du_ISVP = LWFBrook90.EVP.INTER24(p_fT_SFAL, 0, p_fu_LAI, p_fu_SAI, p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS, p_DURATN, u_INTS, MONTHN)
        else
            # else if PINT>0 (potential evaporation of snow intercep. storage)
            # or   if PINT<0 (condensation) and it is cold enough that it deposits as frost:
            aux_du_SINT, aux_du_ISVP = LWFBrook90.EVP.INTER24(p_fT_SFAL, p_fu_PINT, p_fu_LAI, p_fu_SAI, p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS, p_DURATN, u_INTS, MONTHN)
        end
        # b) rain interception,  note potential interception rate is p_fu_PINT-aux_du_ISVP (mm/day)
        # interception catch rate, interception evaporation rate (both in mm/day):
        aux_du_RINT, aux_du_IRVP = LWFBrook90.EVP.INTER24(p_fT_RFAL, p_fu_PINT - aux_du_ISVP, p_fu_LAI, p_fu_SAI, p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS, p_DURATN, u_INTR, MONTHN)
    end

    # fraction of time (precip interval) that canopy is wet
    p_fu_WETFR = min(1.0, (aux_du_IRVP + aux_du_ISVP) / p_fu_PINT)
    # end A) Interception by canopy ################

    # throughfall arriving onto snowpack (or soil layer)
    p_fu_RTHR = p_fT_RFAL - aux_du_RINT # rates in mm/h
    p_fu_STHR = p_fT_SFAL - aux_du_SINT # rates in mm/h

    # B) Effect of wet/snow-covered canopy on fluxes ################
    # reduce potential and actual transpiration for fraction of precip interval that canopy is wet
    p_fu_PTRAN = (1.0 - p_fu_WETFR) * p_fu_PTRAN         # only used for cum_d_ptran
    for i = 1:NLAYER
        aux_du_TRANI[i] = (1.0 - p_fu_WETFR) * aux_du_TRANI[i]
        # LWFBrook90 additionally: if(u_aux_PSIM[i] < PsiCrit[i]) FRSS=1.e+20 end # where PsiCrit (! not PSICR) is a cutoff pressure for evaporation but also transpiration
    end
    # end B) Effect of wet/snow-covered canopy on fluxes ################
    
    # C) Effect of snowpack on fluxes ################
    # compute soil evaporation as average of wet and non-wet canopy computations
    if (u_SNOW <= 0 && p_fu_STHR <= 0)
        # no previous snow, no added snow: soil evaporation weighted for p_fu_WETFR
        aux_du_SLVP = p_fu_WETFR * p_fu_GIVP + (1.0 - p_fu_WETFR) * p_fu_GEVP
    else
        # either previous snow, or added snow, or both => prevents soil evaporation
        aux_du_SLVP = 0.0
    end
    # end C) Effect of snowpack on fluxes ################

    # D) Interception by snowpack (RSNO), snowmelt (SMLT) and snow evaporation (SNVP) ################
    # Compute evolution of snowpack (evaporation, acumulation, etc...)
    if (u_SNOW <= 0 && p_fu_STHR <= 0)
        # no previous snow, no added snow
        aux_du_RSNO = 0.0
        aux_du_SNVP = 0.0
        aux_du_SMLT = 0.0
        # u_CC, u_SNOW, u_SNOWLQ remain unchanged
    elseif (u_SNOW <= 0 && p_fu_STHR > 0)
        # no previous snow, some added snow: initialize zero CC and SNOWLQ
        u_CC     = 0.0
        u_SNOWLQ = 0.0
        # snow accumulation and melt
        (u_SNOW, u_CC, u_SNOWLQ, aux_du_RSNO, aux_du_SNVP, aux_du_SMLT) =
                LWFBrook90.SNO.SNOWPACK(
                    p_fu_RTHR, p_fu_STHR, p_fu_PSNVP, p_fu_SNOEN,
                    # States that are overwritten:
                    u_CC, u_SNOW, u_SNOWLQ,
                    p_DTP, p_fT_TA, p_MAXLQF, p_GRDMLT,
                    p_CVICE, p_LF, p_CVLQ)
    else
        # previous snow with/without added snow: evolve snowpack (accumulation & melt)
        # snow accumulation and melt
        (u_SNOW, u_CC, u_SNOWLQ, aux_du_RSNO, aux_du_SNVP, aux_du_SMLT) =
                LWFBrook90.SNO.SNOWPACK(
                    p_fu_RTHR, p_fu_STHR, p_fu_PSNVP, p_fu_SNOEN,
                    # States that are overwritten:
                    u_CC, u_SNOW, u_SNOWLQ,
                    p_DTP, p_fT_TA, p_MAXLQF, p_GRDMLT,
                    p_CVICE, p_LF, p_CVLQ)
    end
    # end D) Interception by snowpack (RSNO), snowmelt (SMLT) and snow evaporation (SNVP) ################

    # Compute net rain arriving at soil surface (together with snow melt SMLT)
    p_fu_RNET = p_fu_RTHR - aux_du_RSNO 
    # aux_du_SMLT: snowmelt arriving at soil surface (together with snow melt RNET)

    return (# compute some fluxes as intermediate results:
            p_fT_SFAL, p_fT_RFAL, p_fu_RNET, p_fu_PTRAN,
            # compute changes in soil water storage:
            aux_du_TRANI, aux_du_SLVP,
            # compute change in interception storage:
            aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP,
            # compute change in snow storage:
            aux_du_RSNO, aux_du_SNVP, aux_du_SMLT, p_fu_STHR,
            # compute updated states:
            u_SNOW, u_CC, u_SNOWLQ)
end

"""
    MSBITERATE()

http://www.ecoshift.net/brook/b90.html 

Subsurface water movement is determined in several to many iterations per precipitation
time-step. Remaining model calculations are done within subroutine MSBITERATE for each
iteration loop.

Net throughfall (RNET) plus snowmelt (SMLT) may:

1. infiltrate into the soil matrix of the surface horizon (INFLI(1)), 
2. infiltrate directly to deeper horizons via vertical macropore flow (INFLI), 
3. go immediately to streamflow via vertical macropore flow followed by downslope pipe flow (BYFLI), 
4. or go immediately to streamflow via impaction on a variable saturated source area (SRFL). 

The fraction of area acting as a saturated source area (SAFRAC) is obtained in subroutine
SRFLFR. Source area flow (SRFL) is obtained as SAFRAC plus impervious area (IMPERV) times
RNET + SMLT. Infiltration rate (SLFL) is RNET + SMLT - SRFL. The fraction of infiltration to
each layer that bypasses the layer and becomes an output via bypass flow (BYFLI) is
calculated in subroutine BYFLFR. For each layer, the downslope flow rate by matrix flow
(DSFLI) is obtained from subroutine DSLOP. In general, one or more of SRFL, BYFL, and DSFL
will be set to zero by the user.

If the water potential difference between layers is less than the parameter DPSIMX, vertical
flow (VRFLI) is zero; otherwise subroutine VERT obtains VRFLI between layers from a weighted
hydraulic conductivity and the water potential difference between the layers. For the bottom
layer, outflow to groundwater is the hydraulic conductivity of the layer times a parameter
(DRAIN), which can vary from 0 to 1. This assumes a gravity potential gradient.

Subroutine INFLOW is called to get net inflow into each layer (NTFLI) using parameter DTIMAX
as a first approximation for iteration time step. The rate of change of matric potential
with water content (DPSIDW) from function FDPSIDW is used with NTFLI in subroutine ITER to
obtain the maximum iteration time step (DTI) allowed by two parameters. The parameters are
DPSIMX and the maximum allowed change in soil water content (DSWMAX). INFLOW is called again
with the new DTI to get the final NTFLI, VRFLI, BYFLI, and matric uptake (INFLI).

Groundwater discharge to streamflow (GWFL) and deep seepage (SEEP) are obtained from
subroutine GWATER. GWFL is simulated as a fixed fraction of groundwater each day and SEEP is
a fixed fraction of GWFL.

Simulated streamflow is the sum of SRFL, BYFL, DSFL, and GWFL. This can be compared with
measured streamflow if that is available.

At the end of each iteration time-step, soil water content of each layer (SWATI) is updated
by adding NTFLI * DTI. Groundwater storage is also updated. Then new soil water variables
are calculated using function FPSIM and subroutine SOILVAR. Available water (AWAT) is
calculated for output as water held in the root zone between field capacity and PSICR.

SWCHEK tests that the SWATI remain between zero and saturation; if not, the program stops.
At the end of each day the water balance is checked and the program stops if there is an
error greater than 0.003 mm. These crashes should only occur if there are parameter or
programming errors.


"""
function MSBITERATE(FLAG_MualVanGen, NLAYER, p_QLAYER, p_soil,
                    # for SRFLFR:
                    u_SWATI, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC,
                    #
                    p_IMPERV, p_fu_RNET, aux_du_SMLT,
                    p_LENGTH_SLOPE, p_DSLOPE,
                    # for DSLOP:
                    p_RHOWG, u_aux_PSIM, p_fu_KK,
                    #
                    u_aux_PSITI, p_DPSIMAX,
                    #
                    p_DRAIN, p_DTP, t, p_DTIMAX,
                    # for INFLOW:
                    p_INFRAC, p_fu_BYFRAC, aux_du_TRANI, aux_du_SLVP,
                    # for FDPSIDW:
                    u_aux_WETNES,
                    # for ITER:
                    p_DSWMAX, u_aux_θ)

    ##########################
    ## On soil surface, partition incoming rain (RNET) and melt water (SMLT)
    # into either above ground source area flow (streamflow, SRFL) or
    # below ground ("infiltrated") input to soil (SLFL)
    # source area flow rate
    if (p_QLAYER > 0)
        SAFRAC=LWFBrook90.WAT.SRFLFR(p_QLAYER, u_SWATI, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC)
    else
        SAFRAC = 0.
    end
    p_fu_SRFL = min(1., (p_IMPERV + SAFRAC)) * (p_fu_RNET + aux_du_SMLT)

    # Derive water supply rate to soil surface:
    p_fu_SLFL = p_fu_RNET + aux_du_SMLT - p_fu_SRFL

    ##########################
    ## Within soil compute flows from layers:
    #  a) downslope, lateral flow from the layers: DSFLI
    aux_du_DSFLI = fill(NaN, NLAYER) # 0.000001 seconds (1 allocation: 144 bytes)
    if (p_LENGTH_SLOPE == 0 || p_DSLOPE == 0)
        # added in Version 4
        aux_du_DSFLI .= 0
    else
        aux_du_DSFLI .= LWFBrook90.WAT.DSLOP(p_DSLOPE, p_LENGTH_SLOPE, p_RHOWG, p_soil, u_aux_PSIM, p_fu_KK)
    end

    #  b) vertical flux between layers: VRFLI (Richards equation)
    aux_du_VRFLI = fill(NaN, NLAYER) # 0.000001 seconds (1 allocation: 144 bytes)
    for i = NLAYER:-1:1
        # vertical flow rates
        if (i < NLAYER)
            if (abs(u_aux_PSITI[i] - u_aux_PSITI[i+1]) < p_DPSIMAX)
                aux_du_VRFLI[i] = 0
            else
                aux_du_VRFLI[i] =
                    LWFBrook90.WAT.VERT(p_fu_KK[i],     p_fu_KK[i+1],
                                            p_soil.p_KSAT[i],      p_soil.p_KSAT[i+1],
                                            p_soil.p_THICK[i],     p_soil.p_THICK[i+1],
                                            u_aux_PSITI[i], u_aux_PSITI[i+1],
                                            p_soil.p_STONEF[i],    p_soil.p_STONEF[i+1],
                                            p_RHOWG)
            end
        else
        # bottom layer i == NLAYER
            if (p_DRAIN > 0.00001)
            # gravity drainage only
                aux_du_VRFLI[NLAYER] = p_DRAIN * p_fu_KK[NLAYER] * (1 - p_soil.p_STONEF[NLAYER])
            else
            # bottom of profile sealed
                aux_du_VRFLI[NLAYER] = 0
            end
        end
    end

    # first approximation on aux_du_VRFLI
    aux_du_VRFLI_1st_approx = aux_du_VRFLI
    #@debug "u_aux_PSITI[1]: $(u_aux_PSITI[1]), sum(u_aux_PSITI): $(sum(u_aux_PSITI))"
    #@debug "a) aux_du_VRFLI_1st_approx[1]: $(aux_du_VRFLI_1st_approx[1]), sum(aux_du_VRFLI_1st_approx): $(sum(aux_du_VRFLI_1st_approx))"

    # first approximation for iteration time step,time remaining or DTIMAX
    DTRI = p_DTP - (t % p_DTP) # Time remaining
    # DTRI = 1.0 - (t % 1.0)   # as p_DTP is 1.0 days in a default simulation
    DTI = min(DTRI, p_DTIMAX)

    # net inflow to each layer including E and T withdrawal adjusted for interception
    aux_du_VRFLI, aux_du_INFLI, aux_du_BYFLI, du_NTFLI =
        LWFBrook90.WAT.INFLOW(NLAYER, DTI, p_INFRAC, p_fu_BYFRAC, p_fu_SLFL, aux_du_DSFLI, aux_du_TRANI,
                                    aux_du_SLVP, p_soil.p_SWATMAX, u_SWATI,
                                    aux_du_VRFLI_1st_approx) # 0.000003 seconds (4 allocations: 576 bytes)

    # limit step size
    #       TODO(bernhard): This could alternatively be achieved with a stepsize limiter callback
    #                       see: https://diffeq.sciml.ai/stable/features/callback_library/#Stepsize-Limiters
    #       TODO(bernhard): Or it could alternatively be set manually using set_proposed_dt!
    #                       see: https://diffeq.sciml.ai/stable/basics/integrator/#SciMLBase.set_proposed_dt!
    #   ITER computes DTI so that the potential difference (due to aux_du_VRFLI)
    #   between adjacent layers does not change sign during the iteration time step
    if (true) # NOTE: when using DiffEq.jl the integrator time step is determined by solve().
        # One might think, that therefore the adaptive time step control of LWFBrook can be deactivated.
        # However, this is not the case. DTI is used in INFLOW() to compute the fluxes aux_du_VRFLI, ... etc.

        DPSIDW = LWFBrook90.KPT.FDPSIDWF(u_aux_WETNES, p_soil) # 0.000004 seconds (3 allocations: 432 bytes)

        DTINEW=LWFBrook90.WAT.ITER(NLAYER, FLAG_MualVanGen, DTI, LWFBrook90.CONSTANTS.p_DTIMIN, DPSIDW,
                                        du_NTFLI, u_aux_PSITI, u_aux_θ, p_DSWMAX, p_DPSIMAX,
                        p_soil)

        # recompute step with updated DTI if needed
        if (DTINEW < DTI)
            # recalculate flow rates with new DTI
            DTI = DTINEW
            aux_du_VRFLI, aux_du_INFLI, aux_du_BYFLI, du_NTFLI =
                LWFBrook90.WAT.INFLOW(NLAYER, DTI, p_INFRAC, p_fu_BYFRAC, p_fu_SLFL, aux_du_DSFLI, aux_du_TRANI,
                                            aux_du_SLVP, p_soil.p_SWATMAX, u_SWATI,
                                            aux_du_VRFLI_1st_approx)
        end
    end

    return (p_fu_SRFL, p_fu_SLFL, aux_du_DSFLI, aux_du_VRFLI, aux_du_INFLI, aux_du_BYFLI, du_NTFLI, DTI)
end


"""
    compute_isotope_U_of_INTS_INTR_SNOW_and_SLFL(
        p_δ2H_PREC, p_δ18O_PREC, p_fT_TADTM, p_fT_VAPPRES,
        # for INTS (in: SINT; out: ISVP):
        u_INTS, aux_du_SINT, aux_du_ISVP, p_DTP, u_δ2H_INTS, u_δ18O_INTS,
        # for INTR (in: RINT; out: IRVP):
        u_INTR, aux_du_RINT, aux_du_IRVP, u_δ2H_INTR, u_δ18O_INTR,
        # for SNOW (in: STHR, RSNO (both δ_PREC); out: SMLT, SNVP (δ_SNOW and fractionated)):
        p_fu_STHR, aux_du_RSNO, aux_du_SMLT, aux_du_SNVP, u_δ2H_SNOW, u_δ18O_SNOW,
        # to compute isotopic signature of soil infiltration: SLFL
        p_fu_RNET)

Computes updated values of states INTS, INTR, and SNOW as well as their isotopic composition.
Compute mixing and evaporative fractionation, and also compute the isotopic composotion of
the resulting flux that infiltrates into the soil: p_fu_δ18O_SLFL, p_fu_δ2H_SLFL

The function is called in the daily callback.

The function returns: p_fu_δ18O_SLFL, p_fu_δ2H_SLFL, as well as:
    u_INTS
    u_δ18O_INTS
    u_δ2H_INTS
    u_INTR
    u_δ18O_INTR
    u_δ2H_INTR
    u_SNOW
    u_δ18O_SNOW
    u_δ2H_SNOW
"""
function compute_isotope_U_of_INTS_INTR_SNOW_and_SLFL(
    p_δ2H_PREC, p_δ18O_PREC, p_fT_TADTM, p_fT_VAPPRES,
    # for INTS (in: SINT; out: ISVP):
    u_INTS, aux_du_SINT, aux_du_ISVP, p_DTP, u_δ2H_INTS, u_δ18O_INTS,
    # for INTR (in: RINT; out: IRVP):
    u_INTR, aux_du_RINT, aux_du_IRVP, u_δ2H_INTR, u_δ18O_INTR,
    # for SNOW (in: STHR, RSNO (both δ_PREC); out: SMLT, SNVP (δ_SNOW and fractionated)):
    u_SNOW, u_SNOW_MSBupdate, p_fu_STHR, aux_du_RSNO, aux_du_SMLT, aux_du_SNVP, u_δ2H_SNOW, u_δ18O_SNOW,
    # to compute isotopic signature of soil infiltration: SLFL
    p_fu_RNET)
    # check chart "../docs/src/assets/b90flow.gif"

    ##################
    ##################
    ##################
    # Define conditions for isotope calculations:
    Tc = p_fT_TADTM  # °C, average daytime air temperature
    h = min(1.0, p_fT_VAPPRES / LWFBrook90.PET.ESAT(Tc)[1]) # -, relative humidity of the atmosphere (vappress_atm/1 atm)
    γ = 1.0          # -, thermodynamic activity coefficient of evaporating water
    X_INTS = 0.5  # -, turbulence incex of the atmosphere above the evaporating water
    X_INTR = 0.5  # -, turbulence incex of the atmosphere above the evaporating water
    X_SNOW = 1.0  # -, turbulence incex of the atmosphere above the evaporating water

    # 1c) Atmospheric vapor composition assumed to be in equilibrium with precipitation
    #     Benettin-2018-Hydrol_Earh_Syst_Sci citing Gibson et al. 2008
    δₐ(δ_p, α_eq) = (δ_p - 1000 * (α_eq - 1))/α_eq
    δ²H_a  = δₐ(p_δ2H_PREC,  LWFBrook90.ISO.α²H_eq(Tc))
    δ¹⁸O_a = δₐ(p_δ18O_PREC, LWFBrook90.ISO.α¹⁸O_eq(Tc))

    # # For soil:
    # Xa = 0.5 # between molecular and turbulent
    # Xs = 1.0 # molecular only
    # X_SOIL = ((θ_s-θ_res)*Xa + (θ_sat-θ_s)*Xs)/(θ_sat-θ_res)
    # # Taken from Zhou-2021-Environ_Model_Softw (X is called n_k there)
    ##################
    ##################
    ##################


                # # 2a) INTS (in: SINT*δ_SINT; out: ISVP*δ_ISVP)
                # #          with δ_SINT = δ_PREC; δ_ISVP = f(f, α, ...)
                # # Operator step 1
                # @assert aux_du_SINT >= 0 "aux_du_SINT should not be negative"
                # if ((u_INTS == 0) & (aux_du_SINT == 0)) # initially no intercepted snow and no new is added
                #     u_δ18O_INTS_final = 0
                #     u_δ2H_INTS_final  = 0
                #     u_INTS_final      = 0
                # else
                #     if ((u_INTS == 0) & (aux_du_SINT > 0)) # initially no intercepted snow but some is added
                #         u_δ18O_INTS_first = p_δ18O_PREC
                #         u_δ2H_INTS_first  = p_δ2H_PREC
                #     else # u_INTS is not zero, and some/none is added (aux_du_SINT=0 or aux_du_SINT>0)
                #         u_δ18O_INTS_first = u_δ18O_INTS + aux_du_SINT*p_DTP/u_INTS * (p_δ18O_PREC-u_δ18O_INTS)
                #         u_δ2H_INTS_first  = u_δ2H_INTS  + aux_du_SINT*p_DTP/u_INTS * (p_δ2H_PREC -u_δ2H_INTS )
                #     end
                #     u_INTS_first = u_INTS + aux_du_SINT*p_DTP
                #     # Operator step 2: now taking care of aux_du_ISVP
                #     u_INTS_final = u_INTS + (aux_du_SINT - aux_du_ISVP) * p_DTP
                #     f_INTS = max(0, u_INTS_final/u_INTS_first) # fraction remaining after evaporation
                #     # ε_δ18O = 1/(LWFBrook90.ISO.α¹⁸O_eq(Tc) * LWFBrook90.ISO.α¹⁸O_dif^X_INTS) - 1
                #     # ε_δ2H  = 1/(LWFBrook90.ISO.α²H_eq(Tc)  * LWFBrook90.ISO.α²H_dif^X_INTS)  - 1
                #     # u_δ18O_INTS_final = 1000 * ( (1 + u_δ18O_INTS_first/1000)(f_INTS)^ε_δ18O - 1 )
                #     # u_δ2H_INTS_final  = 1000 * ( (1 + u_δ2H_INTS_first /1000)(f_INTS)^ε_δ2H  - 1 )
                #     u_δ18O_INTS_final = u_δ18O_INTS_first#1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ18O_INTS_first/1000, δ¹⁸O_a/1000, f_INTS, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_INTS)
                #     u_δ2H_INTS_final  = u_δ2H_INTS_first#1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ2H_INTS_first /1000,  δ²H_a/1000, f_INTS, h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_INTS)
                # end

                # # 2b) INTR (in: RINT*δ_RINT; out: IRVP*δ_IRVP)
                # #          with δ_RINT = δ_PREC; δ_IRVP = f(f, α, ...)
                # # Operator step 1
                # @assert aux_du_RINT >= 0 "aux_du_RINT should not be negative"
                # if ((u_INTR == 0) & (aux_du_RINT == 0)) # initially no intercepted rain and no new is added
                #     u_δ18O_INTR_final = 0
                #     u_δ2H_INTR_final  = 0
                #     u_INTR_final      = 0
                # else
                #     if ((u_INTR == 0) & (aux_du_RINT > 0)) # initially no intercepted rain but some is added
                #         u_δ18O_INTR_first = p_δ18O_PREC
                #         u_δ2H_INTR_first  = p_δ2H_PREC
                #     else # u_INTR is not zero, and some/none is added (aux_du_RINT=0 or aux_du_RINT>0)
                #         u_δ18O_INTR_first = u_δ18O_INTR + aux_du_RINT*p_DTP/u_INTR * (p_δ18O_PREC-u_δ18O_INTR)
                #         u_δ2H_INTR_first  = u_δ2H_INTR  + aux_du_RINT*p_DTP/u_INTR * (p_δ2H_PREC -u_δ2H_INTR )
                #     end
                #     u_INTR_first = u_INTR + aux_du_RINT*p_DTP
                #     # Operator step 2: now taking care of aux_du_IRVP
                #     u_INTR_final = u_INTR + (aux_du_RINT - aux_du_IRVP) * p_DTP
                #     f_INTR = max(0, u_INTR_final/u_INTR_first) # fraction remaining after evaporation
                #     # ε_δ18O = 1/(LWFBrook90.ISO.α¹⁸O_eq(Tc) * LWFBrook90.ISO.α¹⁸O_dif^X_INTR) - 1
                #     # ε_δ2H  = 1/(LWFBrook90.ISO.α²H_eq(Tc)  * LWFBrook90.ISO.α²H_dif^X_INTR)  - 1
                #     # u_δ18O_INTR_final = 1000 * ( (1 + u_δ18O_INTR_first/1000)(f_INTR)^ε_δ18O - 1 )
                #     # u_δ2H_INTR_final  = 1000 * ( (1 + u_δ2H_INTR_first /1000)(f_INTR)^ε_δ2H  - 1 )
                #     u_δ18O_INTR_final = u_δ18O_INTR_first# 1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ18O_INTR_first/1000, δ¹⁸O_a/1000, f_INTR, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_INTR)
                #     u_δ2H_INTR_final  = u_δ2H_INTR_first# 1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ2H_INTR_first /1000,  δ²H_a/1000, f_INTR, h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_INTR)
                # end

                # # 2c) SNOW (in: STHR*δ_STHR, RSNO*δ_RSNO; out: SMLT*δ_SMLT, SNVP*δ_SNVP)
                # #          with δ_STHR = δ_PREC, δ_RSNO = δ_PREC; δ_SMLT = δ_SNOW, δ_SNVP = f(f, α, ...)

                # # NOTE: for SNOW the isotope balance is greatly simplified. The most precise
                # #       approach would be to define mixing concentrattion within `SNOWPACK()` in `module_SNO.jl`
                # # Operator step 1
                # @assert p_fu_STHR + aux_du_RSNO >= 0 "p_fu_STHR + aux_du_RSNO should not be negative"
                # if ((u_SNOW == 0) & (p_fu_STHR + aux_du_RSNO == 0)) # initially no snowpack and no new is added
                #     u_δ18O_SNOW_final = 0
                #     u_δ2H_SNOW_final  = 0
                #     u_SNOW_final      = 0
                # else
                #     if ((u_SNOW == 0) & (p_fu_STHR + aux_du_RSNO > 0)) # initially no snowpack but some is added
                #         u_δ18O_SNOW_first = p_δ18O_PREC
                #         u_δ2H_SNOW_first  = p_δ2H_PREC
                #     else # u_SNOW is not zero, and some/none is added ((p_fu_STHR + aux_du_RSNO)=0 or (p_fu_STHR + aux_du_RSNO)>0)
                #         # p_fu_STHR, aux_du_RSNO, aux_du_SMLT, aux_du_SNVP
                #         u_δ18O_SNOW_first = u_δ18O_SNOW + (p_fu_STHR + aux_du_RSNO)*p_DTP/u_SNOW * (p_δ18O_PREC - u_δ18O_SNOW)
                #                                     # NOTE: because the outflow term for aux_du_SMLT has
                #                                     #       an isotope concentration of u_δ18O_SNOW and it is thus not needed:
                #                                     # - (aux_du_SMLT)*p_DTP/u_SNOW * (u_δ18O_SNOW-u_δ18O_SNOW)
                #         u_δ2H_SNOW_first = u_δ2H_SNOW + (p_fu_STHR + aux_du_RSNO)*p_DTP/u_SNOW * (p_δ2H_PREC - u_δ2H_SNOW)
                #                                     # NOTE: because the outflow term for aux_du_SMLT has
                #                                     #       an isotope concentration of u_δ2H_SNOW and it is thus not needed:
                #                                     # - (aux_du_SMLT)*p_DTP/u_SNOW * (u_δ2H_SNOW-u_δ2H_SNOW)
                #     end
                #     u_SNOW_first = u_SNOW + p_DTP * (p_fu_STHR + aux_du_RSNO - aux_du_SMLT)
                #     # Operator step 2: now taking care of aux_du_SNVP
                #     u_SNOW_final = u_SNOW + p_DTP * (p_fu_STHR + aux_du_RSNO - aux_du_SMLT - aux_du_SNVP)
                #     f_SNOW = max(0, u_SNOW_final/u_SNOW_first) # fraction remaining after evaporation
                #     # ε_δ18O = 1/(LWFBrook90.ISO.α¹⁸O_eq(Tc) * LWFBrook90.ISO.α¹⁸O_dif^X_SNOW) - 1
                #     # ε_δ2H  = 1/(LWFBrook90.ISO.α²H_eq(Tc)  * LWFBrook90.ISO.α²H_dif^X_SNOW)  - 1
                #     # u_δ18O_SNOW_final = 1000 * ( (1 + u_δ18O_SNOW_first/1000)(f_SNOW)^ε_δ18O - 1 )
                #     # u_δ2H_SNOW_final  = 1000 * ( (1 + u_δ2H_SNOW_first /1000)(f_SNOW)^ε_δ2H  - 1 )
                #     u_δ18O_SNOW_final = u_δ18O_SNOW_first#1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ18O_SNOW_first/1000, δ¹⁸O_a/1000, f_SNOW, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_SNOW)
                #     u_δ2H_SNOW_final  = u_δ2H_SNOW_first#1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ2H_SNOW_first /1000,  δ²H_a/1000, f_SNOW, h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_SNOW)
                # end

    # 2a) u_INTS (in: aux_du_SINT*δ_SINT; out: aux_du_ISVP*δ_ISVP)
    #            with δ_SINT = δ_PREC; δ_ISVP = f(f, α, ...)
    # 2b) u_INTR (in: aux_du_RINT*δ_RINT; out: aux_du_IRVP*δ_IRVP)
    #            with δ_RINT = δ_PREC; δ_IRVP = f(f, α, ...)
    # 2c) SNOW (in: p_fu_STHR*δ_STHR, aux_du_RSNO*δ_RSNO; out: aux_du_SMLT*δ_SMLT, aux_du_SNVP*δ_SNVP)
    #          with δ_STHR = δ_PREC, δ_RSNO = δ_PREC; δ_SMLT = δ_SNOW, δ_SNVP = f(f, α, ...)





    # update_δ_with_mixing_and_evaporation(dt, u₀, δ₀, inflow, δin, outflow, E, δₐ, h, α_eq, α_dif, γ, X)
    u_INTS_final, u_δ2H_INTS_final  = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_INTS, u_δ2H_INTS,  aux_du_SINT,               p_δ2H_PREC,                0,           LWFBrook90.ISO.R_VSMOW²H,  aux_du_ISVP, δ²H_a,  h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_INTS)
    _,            u_δ18O_INTS_final = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_INTS, u_δ18O_INTS, aux_du_SINT,               p_δ18O_PREC,               0,           LWFBrook90.ISO.R_VSMOW¹⁸O, aux_du_ISVP, δ¹⁸O_a, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_INTS)
    u_INTR_final, u_δ2H_INTR_final  = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_INTR, u_δ2H_INTR,  aux_du_RINT,               p_δ2H_PREC,                0,           LWFBrook90.ISO.R_VSMOW²H,  aux_du_IRVP, δ²H_a,  h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_INTR)
    _,            u_δ18O_INTR_final = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_INTR, u_δ18O_INTR, aux_du_RINT,               p_δ18O_PREC,               0,           LWFBrook90.ISO.R_VSMOW¹⁸O, aux_du_IRVP, δ¹⁸O_a, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_INTR)
    # NOTE: for SNOW the isotope balance is greatly simplified. The most precise
    #       approach would be to define mixing concentrattion within `SNOWPACK()` in `module_SNO.jl`
    u_SNOW_final, u_δ2H_SNOW_final  = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_SNOW, u_δ2H_SNOW,  (p_fu_STHR, aux_du_RSNO), (p_δ2H_PREC,  p_δ2H_PREC),  aux_du_SMLT, LWFBrook90.ISO.R_VSMOW²H,  aux_du_SNVP, δ²H_a,  h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_SNOW)
    _,            u_δ18O_SNOW_final = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_SNOW, u_δ18O_SNOW, (p_fu_STHR, aux_du_RSNO), (p_δ18O_PREC, p_δ18O_PREC), aux_du_SMLT, LWFBrook90.ISO.R_VSMOW¹⁸O, aux_du_SNVP, δ¹⁸O_a, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_SNOW)



    # 3) also compute δ_SLFL as mix of δ_SMLT with δ_RNET (i.e. water that infiltrates)
    # if (aux_du_SMLT + p_fu_RNET == 0)
    #     p_fu_δ18O_SLFL = δ18O_empty
    #     p_fu_δ2H_SLFL  = δ2H_empty
    # elseif (aux_du_SMLT == 0)
    #     p_fu_δ18O_SLFL = p_δ18O_PREC
    #     p_fu_δ2H_SLFL  = p_δ2H_PREC
    # elseif (p_fu_RNET == 0)
    #     p_fu_δ18O_SLFL = u_δ18O_SNOW_final  # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
    #     p_fu_δ2H_SLFL  = u_δ2H_SNOW_final  # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
    # else # both fluxes are non-null and we need to compute their mix
    #     p_fu_δ18O_SLFL = (u_δ18O_SNOW_final * aux_du_SMLT + p_δ18O_PREC * p_fu_RNET) / (aux_du_SMLT + p_fu_RNET)  # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
    #     p_fu_δ2H_SLFL  = (u_δ2H_SNOW_final * aux_du_SMLT  + p_δ2H_PREC * p_fu_RNET)  / (aux_du_SMLT + p_fu_RNET)  # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
    # end
    # TODO(bernhard): deactivate the following workaround and activate above code
    p_fu_δ18O_SLFL = p_δ18O_PREC
    p_fu_δ2H_SLFL  = p_δ2H_PREC

    return (p_fu_δ18O_SLFL, p_fu_δ2H_SLFL,
        u_INTS_final, u_δ18O_INTS_final, u_δ2H_INTS_final,
        u_INTR_final, u_δ18O_INTR_final, u_δ2H_INTR_final,
        u_SNOW_final, u_δ18O_SNOW_final, u_δ2H_SNOW_final)

    # In-place modify
    # u_INTS      = u_INTS_final
    # u_δ18O_INTS = u_δ18O_INTS_final
    # u_δ2H_INTS  = u_δ2H_INTS_final
    # u_INTR      = u_INTR_final
    # u_δ18O_INTR = u_δ18O_INTR_final
    # u_δ2H_INTR  = u_δ2H_INTR_final
    # u_SNOW      = u_SNOW_final
    # u_δ18O_SNOW = u_δ18O_SNOW_final
    # u_δ2H_SNOW  = u_δ2H_SNOW_final
end

"""
    compute_isotope_du_GWAT_SWATI(
        # for GWAT:
        u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT,
        # for SWATI:
        du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI, # (non-fractionating)
        aux_du_SLVP, p_fT_TADTM, p_fT_VAPPRES, p_δ2H_PREC, p_δ18O_PREC, u_aux_WETNES, # (fractionating)
        u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, EffectiveDiffusivity_18O, EffectiveDiffusivity_2H,
        )

Computes rate of change (du) of isotopic composition of states:
δ18O_GWAT, δ2H_GWAT, δ18O_SWATI, δ2H_SWATI
By computing the rates of change this function is needed when we track the isotopic composition
by numerically solving an ODE. Alternatively we can use a callback and updated the isotopic
compositions separately from solving the flow equation (by using `compute_isotope_u_GWAT_SWATI`).

The function returns: du_δ18O_GWAT, du_δ2H_GWAT, du_δ18O_SWATI, du_δ2H_SWATI that can be used
for numerically solving the evolution of state (either in function f or in cb by using a simple
time stepping algorithm).
"""
function compute_isotope_du_GWAT_SWATI(
    # for GWAT:
    u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT,
    # for SWATI:
    du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI, # (non-fractionating)
    aux_du_SLVP, p_fT_TADTM, p_fT_VAPPRES, p_δ2H_PREC, p_δ18O_PREC, u_aux_WETNES, # (fractionating)
    u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, EffectiveDiffusivity_18O, EffectiveDiffusivity_2H,
    )

    return compute_isotope_GWAT_SWATI("numerical-du", NaN,
        # for GWAT:
        u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT,NaN,NaN,
        # for SWATI:
        du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI, # (non-fractionating)
        aux_du_SLVP, p_fT_TADTM, p_fT_VAPPRES, p_δ2H_PREC, p_δ18O_PREC, u_aux_WETNES, # (fractionating)
        u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, EffectiveDiffusivity_18O, EffectiveDiffusivity_2H,
        )
end

"""
    compute_isotope_u_GWAT_SWATI(integrator,
        # for GWAT:
        u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT, du_GWFL, du_SEEP,
        # for SWATI:
        du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI, # (non-fractionating)
        aux_du_SLVP, p_fT_TADTM, p_fT_VAPPRES, p_δ2H_PREC, p_δ18O_PREC, u_aux_WETNES, # (fractionating)
        u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, EffectiveDiffusivity_18O, EffectiveDiffusivity_2H,
        )

Computes new states of isotopic composition of GWAT and SWATI by using an operator splitting
approach between the non-fractionating fluxes (mixing) and fractionating fluxes (evaporation)
in the water balance.
This allows to use analytical solutions by considering the effects to be separated. It computes
first an auxiliary state vector after all the mixing fluxes are considered and then for each
water pool evolution of the isotopic composition due to fractionating fluxes is considered
by assuming the fractionating flux is the only outgoing flux during the entire time step.

When no operator splitting is desired, the user could use the numerical approach that
considers all of the fluxes to modify the isotopic compositions simultaneously. This is done
by using the alternative function `compute_isotope_du_GWAT_SWATI()`

The function returns: u_δ18O_GWAT, u_δ2H_GWAT, u_δ18O_SWATI, u_δ2H_SWATI
"""
function compute_isotope_u_GWAT_SWATI(integrator,
    # TODO(bernhard): make non-allocating and rename as compute_isotope_u_GWAT_SWATI!()
    # for GWAT:
    u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT, du_GWFL, du_SEEP,
    # for SWATI:
    du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI, # (non-fractionating)
    aux_du_SLVP, p_fT_TADTM, p_fT_VAPPRES, p_δ2H_PREC, p_δ18O_PREC, u_aux_WETNES, # (fractionating)
    u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, EffectiveDiffusivity_18O, EffectiveDiffusivity_2H,
    )

    # Code below is very condensed. Refer to implementation notes in programming comments
    # noted within the function `compute_isotope_du_GWAT_SWATI()`.

    return compute_isotope_GWAT_SWATI("analytical-u", integrator,
        # for GWAT:
        u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT, du_GWFL, du_SEEP,
        # for SWATI:
        du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI, # (non-fractionating)
        aux_du_SLVP, p_fT_TADTM, p_fT_VAPPRES, p_δ2H_PREC, p_δ18O_PREC, u_aux_WETNES, # (fractionating)
        u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, EffectiveDiffusivity_18O, EffectiveDiffusivity_2H,
        )
end


function compute_isotope_GWAT_SWATI(solution_type, integrator,
        # for GWAT:
        u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT, du_GWFL, du_SEEP,
        # for SWATI:
        du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI, # (non-fractionating)
        aux_du_SLVP, p_fT_TADTM, p_fT_VAPPRES, p_δ2H_PREC, p_δ18O_PREC, u_aux_WETNES, # (fractionating)
        u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, EffectiveDiffusivity_18O, EffectiveDiffusivity_2H,
        )

    @assert solution_type == "numerical-du" || solution_type == "analytical-u" """
        Unknown solution_type ($solution_type) provided to function compute_isotope_GWAT_SWATI.
        Should be either "analytical-u" or "numerical-du".
        """

    NLAYER = length(aux_du_VRFLI)
    # Theory:
        # MASS BALANCE WATER
        # Actual water mass balance equation over computational cell between interfaces z_upper and z_lower:
        # ∂/∂t ∫θ dz = ∫∂/∂z[K_h(∂h/∂z + 1)] dz - ∫S dz  # (where S are source/sink terms)
        # ∂/∂t SWATI = [K_h(∂h/∂z + 1)]_(z_lower)^(z_upper) - ∫S dz
        #            = [K_h(∂h/∂z + 1)]_(z_lower)^(z_upper) + (INFLI-TRANI-DSFL-SLVP)
        #            = q(z_upper) - q(z_lower) + INFLI - TRANI - DSFL - SLVP

        # MASS BALANCE SOLUTE (i.e. HEAVY ISOTOPE)
        # Actual heavy isotope mass balance equation over computational cell:
        # ∂/∂t ∫θ Cᵢ dz = Dᵢ(z_upper)∂Cᵢ/∂z(z_upper) - Dᵢ(z_lower)∂Cᵢ/∂z(z_lower)
        #                 - q(z_upper) Cᵢ(z_upper) + q(z_lower) Cᵢ(z_lower)
        #                 + (INFLI*Cᵢ_{INFLI} - TRANI*Cᵢ_{TRANI} - DSFL*Cᵢ_{DSFL} - SLVP*Cᵢ_{SLVP})
        #               = r-h-s
        # with constant θ and Cᵢ over cell:
        # ∂/∂t θ Cᵢ ∫dz = ∂/∂t θ Cᵢ * (thickness) = thickness * [ θ * ∂/∂t Cᵢ + Cᵢ * ∂/∂t θ ] # where thickness = z_upper - z_lower
        #                                         = [ SWATI * ∂/∂t Cᵢ + Cᵢ * ∂/∂t SWATI ]
        #                                         = r-h-s
        # ==>
        # ∂/∂t Cᵢ = [- Cᵢ/SWATI * ∂/∂t SWATI] + 1/SWATI * [
        #                                                  diff(z_upper) - diff(z_lower) - qCᵢ(z_upper) + qCᵢ(z_lower) +
        #                                                  INFLI*Cᵢ_{INFLI} - TRANI*Cᵢ_{TRANI} - DSFL*Cᵢ_{DSFL} - SLVP*Cᵢ_{SLVP}
        #                                                 ]
        # For mixing fluxes the compositions Cᵢ are like the originating storage.
        # For diffusion fluxes the compositions Cᵢ are computed based on the diffusivity.
        # For evaporation fluxes the compositions Cᵢ are computed using the Craig-Gordon model:
        #   Cᵢ = [ (δᵢ - ε_eq)/α_eq - h*δₐ - ε_kin ] / [ 1 - h + ε_kin/1000] # [‰]
        # source: Benettin 2018 HESS eq. 1 and Gibson 2016 QSR eq. 3
        #         and Craig Gordon 1965 eq. 23 (and text above)
        #
        # Above equation of the Craig-Gordon model is not used, instead one is derived from Gonfiantini 2018:
        #   4) R_esc = R_w / ( α * (α_dif)^X )
        #   5) R_cap = R_A / (α_dif)^X
        #   6) Evaporation E = (γ - h) φ_vap
        #   7a) Upward isotopic flux: E_isotop = φ_esc * R_esc - φ_cap * R_cap
        #                                      = γ*φ_vap * R_esc - h*φ_vap * R_cap
        #                                      = (γ*R_esc - h*R_cap) * φ_vap
        #       Net evaporating flux is thus: R_E = E_isotop/E
        #                                         = (γ*R_esc - h*R_cap) / (γ - h)
        #                                         = (γ*R_w / ( α * (α_dif)^X ) - h*R_A / (α_dif)^X) / (γ - h)
        #                                         = (γ*R_w/α - h*R_A) / ((γ - h) * (α_dif)^X)
        #       And in delta notation this is:
        #                                     δ_E = 1000*(1 + 1/((γ - h)*(α_dif)^X) * (γ/α*(1+δ_w/1000) - h*(1+δ_A/1000)))
        #
        # This can be used either:
        #  a) to solve the evolution of δ numerically considering all the in and outflows at the same time
        #  b) solve the evolution of δ analytically (assuming a single outflow due to evaporation) (Note that we could use
        #     operator splitting to compute first the new δ due to non-fractionating in-/outflows and then the fractionating
        #     evaporation). ==> R_w = R_{w,0} + (A/B*R_A)*f^B -A/B*R_A,
        #                   where: B = γ / ( α * (α_dif)^X * (γ - h) ) -1
        #                   where: A/B = - h * α / (γ - ( α * (α_dif)^X * (γ - h) ))
        #                   and unused: A = - h / ((α_dif)^X * (γ - h))

        # Using approach a): the mass balance equation above stated for ∂/∂t Cᵢ is transormed to ∂/∂t δᵢ and solved numerically
        # Using approach b): the mass balance equation above is split into fractionating and non-fractionating fluxes:
        #                    non-fractionating: ∂/∂t Cᵢ = [- Cᵢ/SWATI * ∂/∂t SWATI] + 1/SWATI * [ diff(z_upper) - diff(z_lower) - qCᵢ(z_upper) + qCᵢ(z_lower) + INFLI*Cᵢ_{INFLI} - TRANI*Cᵢ_{TRANI} - DSFL*Cᵢ_{DSFL}                  ]
        #                    fractionating    :         = [- Cᵢ/SWATI * ∂/∂t SWATI] + 1/SWATI * [                                                                                                                    - SLVP*Cᵢ_{SLVP} ]
        #       non-fractionating fluxes are solved with a simple forward Euler scheme (note that this is an approximation for the diffusive fluxes)
        #       fractionating fluxes are solved with analytical integration of a reservoir with a single (fractionating) sink: R_w = ...




    # Note that above mass balance equation is expressed with C, which is exact if we take for C
    # the isotope-amount fraction x, a.k.a. isotopic abundance, a.k.a. atom fraction
    # (see Coplen-2011-Rapid_Commun_Mass_Spectrom). Using δ would introduce an error.
    #
    # δ and x are related in the following way:
    #   x = 1 / (1 + 1/(R_std*(δ+1)) )
    #   δ = 1 / ( R_std*(1/x - 1) ) - 1, respectively

    δ_to_x(δ_permil,R_std) = 1 ./ (1 .+ 1 ./ (R_std .* ( δ_permil./1000 .+ 1 )) )
    x_to_δ(x,R_std) = (1 ./ ( R_std .* (1 ./ x .- 1) ) .- 1) .* 1000 # in permil

    dxdt_to_dδdt(dxdt, x, R_std) = dxdt .* 1 ./ R_std .* 1 ./ (x .- 1).^2 .* 1000  # using dδ/dt = dδ/dx * dx/dt

    # where R_std are:
    R_VSMOW¹⁸O = 2005.2e-6 # (source: Baertschi-1976-Earth_Planet_Sci_Lett)
    R_VSMOW²H  = 155.76e-6 # (source: Hagemann-1970-Tellus)

    ##################
    ##################
    ##################
    # 0) Define conditions for isotope calculations:
    Tc = p_fT_TADTM  # °C, average daytime air temperature
    h = min(1.0, p_fT_VAPPRES / LWFBrook90.PET.ESAT(Tc)[1]) # -, relative humidity of the atmosphere (vappress_atm/1 atm)
    γ = 1.0          # -, thermodynamic activity coefficient of evaporating water
    # X_INTS = 0.5  # -, turbulence incex of the atmosphere above the evaporating water
    # X_INTR = 0.5  # -, turbulence incex of the atmosphere above the evaporating water
    # X_SNOW = 1.0  # -, turbulence incex of the atmosphere above the evaporating water

    # Atmospheric vapor composition assumed to be in equilibrium with precipitation
    #     Benettin-2018-Hydrol_Earh_Syst_Sci citing Gibson et al. 2008
    δₐ(δ_p, α_eq) = (δ_p - 1000 * (α_eq - 1))/α_eq
    δ²H_a  = δₐ(p_δ2H_PREC,  LWFBrook90.ISO.α²H_eq(Tc))
    δ¹⁸O_a = δₐ(p_δ18O_PREC, LWFBrook90.ISO.α¹⁸O_eq(Tc))

    # For soil:
    Xa = 0.5 # between molecular and turbulent
    Xs = 1.0 # molecular only
    X_SOIL = u_aux_WETNES[1] * Xa + (1 - u_aux_WETNES[1]) * Xs
    # when fully saturated, WETNES == 1, water vapor leaves quickly to atmosphere
    # when hardly saturated, WETNES -> 0, water vapor crosses soil pores until it reaches atmosphere
    # Equivalent:
    # X_SOIL = ((u_aux_θ[1]-θ_res)*Xa + (θ_sat-u_aux_θ[1])*Xs)/(θ_sat-θ_res)
    # Taken from Zhou-2021-Environ_Model_Softw (X is called n_k there)
    ##################
    ##################
    ##################


    ##################
    ##################
    ##################
    # 1) Compute isotopic composition (or solute concentrations) of evaporating soil water
    ### a) Fix constants
    α¹⁸O_eq = LWFBrook90.ISO.α¹⁸O_eq(Tc)
    α²H_eq  = LWFBrook90.ISO.α²H_eq(Tc)
    # ε¹⁸O_dif = (α¹⁸O_dif-1)*1000  # diffusive (i.e. kinetic fractionation)
    # ε²H_dif  = (α²H_dif-1)*1000   # diffusive (i.e. kinetic fractionation)
    # ε¹⁸O_eq = (α¹⁸O_eq-1)*1000    # equilibrium fractionation
    # ε²H_eq  = (α²H_eq-1)*1000     # equilibrium fractionation

    if solution_type == "numerical-du"
        ### b1) Compute δ signature of evaporating flux (for numerical solution)
        # Equation derived based on Gonfiantini (see 60 lines above in comment)
        # δ_E = 1000*(1 + 1/((γ - h)*(α_dif)^X) * (γ/α*(1+δ_w/1000) - h*(1+δ_A/1000)))
        δ¹⁸O_SLVP = 1000*( 1 + (γ/α¹⁸O_eq*(1 + u_δ18O_SWATI[1] / 1000) - h*(1 + δ¹⁸O_a/1000)) /
                                ((γ - h)*(LWFBrook90.ISO.α¹⁸O_dif)^X_SOIL)
        )
        δ²H_SLVP  = 1000*( 1 + (γ/α²H_eq*(1 + u_δ2H_SWATI[1] / 1000) - h*(1 + δ²H_a/1000)) /
                                ((γ - h)*(LWFBrook90.ISO.α²H_dif)^X_SOIL)
        )
        # TODO(bernhard): for debugging:
        δ¹⁸O_SLVP = u_δ18O_SWATI[1] # disabling evaporation fractionation
        δ²H_SLVP  = u_δ2H_SWATI[1]  # disabling evaporation fractionation
        # (Above is an alternative to formulation in Benettin 2018 HESS eq. 1 and Gibson 2016)
        # Cᵢ_SLVP = ( (Cᵢ - ε¹⁸O_eq)/α¹⁸O_eq - h*δ¹⁸O_a - ε¹⁸O_dif ) /
        #             (1 - h + ε¹⁸O_dif/1000) # [‰]

    else
        # TODO(bernhard): current analytical code doesn't take diffusive fluxes into account
        dt = integrator.dt # should be equal to integrator.tprev - t

        ### b2) Compute δ signature taking into account effect of evaporating flux (for analytical solution)
        # Equation derived based on Gonfiantini (see 60 lines above in comment)
        #   R_w = R_{w,0} + (A/B*R_A)*f^B -A/B*R_A,
        #                   where: B = γ / ( α * (α_dif)^X * (γ - h) ) -1
        #                   where: A/B = - h * α / (γ - ( α * (α_dif)^X * (γ - h) ))
        # B¹⁸O = γ / ( α¹⁸O_eq * (LWFBrook90.ISO.α¹⁸O_dif)^X_SOIL * (γ - h) ) -1
        # B²H  = γ / ( α²H_eq  * (LWFBrook90.ISO.α²H_dif )^X_SOIL * (γ - h) ) -1
        # AdivB¹⁸O = - h * α¹⁸O_eq / (γ - ( α¹⁸O_eq * (LWFBrook90.ISO.α¹⁸O_dif)^X_SOIL * (γ - h) ))
        # AdivB²H  = - h * α¹⁸O_eq / (γ - ( α¹⁸O_eq * (LWFBrook90.ISO.α²H_dif )^X_SOIL * (γ - h) ))

        # Precise: R_w = R_{w,0} + (A/B*R_A)*f^B -A/B*R_A,
        # Approximative (?): δ = (δ0 + 1 + A/B*(δA + 1)) * f^B - (1 + A/B*(δA + 1))

        # # u_SWATI1_before_evaporation_flux = u_SWATI[1] + aux_du_SLVP * dt
        # # f_SWATI1 = aux_du_SLVP / u_SWATI1_before_evaporation_flux
        # f_SWATI1 = aux_du_SLVP / (u_SWATI[1] + aux_du_SLVP * dt) # fraction remaining due to fractionating flux

        # u_δ18O_SWATI[1] = ( u_δ18O_SWATI[1] + 1 + AdivB¹⁸O * (δ¹⁸O_a + 1) ) * f_SWATI1^B¹⁸O - (1 + AdivB¹⁸O * (δ¹⁸O_a + 1) )
        # u_δ2H_SWATI[1]  = ( u_δ2H_SWATI[1]  + 1 + AdivB²H  * (δ²H_a  + 1) ) * f_SWATI1^B²H  - (1 + AdivB²H  * (δ²H_a  + 1) )

        #                    non-fractionating: ∂/∂t Cᵢ = [- Cᵢ/SWATI * ∂/∂t SWATI] + 1/SWATI * [ diff(z_upper) - diff(z_lower) - qCᵢ(z_upper) + qCᵢ(z_lower) + INFLI*Cᵢ_{INFLI} - TRANI*Cᵢ_{TRANI} - DSFL*Cᵢ_{DSFL}                  ]
        #                    fractionating    :         = [- Cᵢ/SWATI * ∂/∂t SWATI] + 1/SWATI * [                                                                                                                    - SLVP*Cᵢ_{SLVP} ]

        u_δ18O_SWATI_final = u_δ18O_SWATI
        u_δ2H_SWATI_final  = u_δ2H_SWATI
        u_δ18O_GWAT_final = NaN
        u_δ2H_GWAT_final  = NaN

        for i in 1:NLAYER
            # u_SWATI1 (in: INFL[1]; out: TRANI[1], VRFL[1], DSFL[1], SLVP)
            #          with δ_INFL[1] = (δ_RNET*RNET + δ_SMLT*SMLT)/(RNET+SMLT);
            #          with δ_SLVP = f(f, α, ...)
            #          with INFL[1] = RNET + SMLT - SRFL - sum(BYFL) - sum(INFL[2:end])
            # u_SWATIi (in: VRFL[i-1], INFL[i]; out: TRANI[i], VRFL[i], DSFL[i])
            #          with δ_INFL[i] = (δ_RNET*RNET + δ_SMLT*SMLT)/(RNET+SMLT);
            #          with δ_VRFL[i-1] = δ_SWATI[i-1]
            # u_SWATIN (in: VRFL[N-1], INFL[N]; out: TRANI[N], VRFL[N], DSFL[N])
            #          with δ_INFL[N] = (δ_RNET*RNET + δ_SMLT*SMLT)/(RNET+SMLT);
            #          with δ_VRFL[N-1] = δ_SWATI[N-1]

            # TODO(bernhard): using mixing with δ-values is only approximative
            if i == 1
                δin18 = δ18O_INFLI
                δin2  = δ2H_INFLI
                inflow  = aux_du_INFLI[1]
                outflow = [aux_du_TRANI[1], aux_du_VRFLI[1], aux_du_DSFLI[1]]
                E       = aux_du_SLVP
                # update_δ_with_mixing_and_evaporation(dt, u₀, δ₀, inflow, δin, outflow, E, δₐ, h, α_eq, α_dif, γ, X)
                _, u_δ18O_SWATI_final[1] = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_SWATI[1], u_δ18O_SWATI[1], inflow, δin18, outflow, LWFBrook90.ISO.R_VSMOW²H,  E, δ¹⁸O_a, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif,  γ, X_SOIL)
                _, u_δ2H_SWATI_final[1]  = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_SWATI[1], u_δ2H_SWATI[1],  inflow, δin2,  outflow, LWFBrook90.ISO.R_VSMOW¹⁸O, E, δ²H_a,  h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,   γ, X_SOIL)
            else
                δin18 = [δ18O_INFLI, u_δ18O_SWATI[i-1]]
                δin2  = [δ2H_INFLI , u_δ2H_SWATI[i-1] ]
                inflow  = [aux_du_INFLI[i], aux_du_VRFLI[i-1]]
                outflow = [aux_du_TRANI[i], aux_du_VRFLI[i], aux_du_DSFLI[i]]
                E       = 0
                # update_δ_with_mixing_and_evaporation(dt, u₀, δ₀, inflow, δin, outflow, E, δₐ, h, α_eq, α_dif, γ, X)
                _, u_δ18O_SWATI_final[i] = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_SWATI[i], u_δ18O_SWATI[i], inflow, δin18, outflow, LWFBrook90.ISO.R_VSMOW²H,  E, δ¹⁸O_a, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif,  γ, X_SOIL)
                _, u_δ2H_SWATI_final[i]  = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_SWATI[i], u_δ2H_SWATI[i],  inflow, δin2,  outflow, LWFBrook90.ISO.R_VSMOW¹⁸O, E, δ²H_a,  h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,   γ, X_SOIL)
            end
            # # update_δ_with_mixing_and_evaporation(dt, u₀, δ₀, inflow, δin, outflow, E, δₐ, h, α_eq, α_dif, γ, X)
            # u_INTS_final, u_δ2H_INTS_final  = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_INTS, u_δ2H_INTS,  aux_du_SINT,               p_δ2H_PREC,                0,           aux_du_ISVP, δ²H_a,  h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_INTS)
            # _,            u_δ18O_INTS_final = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_INTS, u_δ18O_INTS, aux_du_SINT,               p_δ18O_PREC,               0,           aux_du_ISVP, δ¹⁸O_a, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_INTS)
            # u_INTR_final, u_δ2H_INTR_final  = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_INTR, u_δ2H_INTR,  aux_du_RINT,               p_δ2H_PREC,                0,           aux_du_IRVP, δ²H_a,  h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_INTR)
            # _,            u_δ18O_INTR_final = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_INTR, u_δ18O_INTR, aux_du_RINT,               p_δ18O_PREC,               0,           aux_du_IRVP, δ¹⁸O_a, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_INTR)
            # # NOTE: for SNOW the isotope balance is greatly simplified. The most precise
            # #       approach would be to define mixing concentrattion within `SNOWPACK()` in `module_SNO.jl`
            # u_SNOW_final, u_δ2H_SNOW_final  = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_SNOW, u_δ2H_SNOW,  (p_fu_STHR, aux_du_RSNO), (p_δ2H_PREC,  p_δ2H_PREC),  aux_du_SMLT, aux_du_SNVP, δ²H_a,  h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_SNOW)
            # _,            u_δ18O_SNOW_final = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_SNOW, u_δ18O_SNOW, (p_fu_STHR, aux_du_RSNO), (p_δ18O_PREC, p_δ18O_PREC), aux_du_SMLT, aux_du_SNVP, δ¹⁸O_a, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_SNOW)
        end

        δin18 = u_δ18O_SWATI[NLAYER]
        δin2  = u_δ2H_SWATI[NLAYER]
        inflow  = aux_du_VRFLI[NLAYER]
        outflow = [du_GWFL, du_SEEP]
        E       = 0
        # update_δ_with_mixing_and_evaporation(dt, u₀, δ₀, inflow, δin, outflow, E, δₐ, h, α_eq, α_dif, γ, X)
        _, u_δ18O_GWAT_final = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_GWAT, u_δ18O_GWAT, inflow, δin18, outflow, LWFBrook90.ISO.R_VSMOW²H,  E, δ¹⁸O_a, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif,  γ, X_SOIL)
        _, u_δ2H_GWAT_final  = LWFBrook90.ISO.update_δ_with_mixing_and_evaporation(dt, u_GWAT, u_δ2H_GWAT,  inflow, δin2,  outflow, LWFBrook90.ISO.R_VSMOW¹⁸O, E, δ²H_a,  h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,   γ, X_SOIL)
    end


    ##################
    ##################
    ##################


    ##################
    ##################
    ##################
    # 2) Compute isotopic composition (or solute concentrations) of all other water fluxes
    if solution_type == "numerical-du"
        # Note: transform the fluxes into isotope-amount fraction x, a.k.a. isotopic abundance, a.k.a. atom fraction
        Cᵢ¹⁸O = δ_to_x.(u_δ18O_SWATI, R_VSMOW¹⁸O)
        Cᵢ²H = δ_to_x.(u_δ2H_SWATI,  R_VSMOW²H)

        diff_upp = 0 #.* similar(aux_du_VRFLI) #TODO(bernhard): change this to include diffusive flux of isotopes, (note it generated an earlier bug when I included it as vector of 0's instead of a scalar of 0)
        diff_low = 0 #.* similar(aux_du_VRFLI) #TODO(bernhard): change this to include diffusive flux of isotopes, (note it generated an earlier bug when I included it as vector of 0's instead of a scalar of 0)
        qCᵢ¹⁸O_upp  = [0; aux_du_VRFLI[1:(NLAYER-1)]] .* [0; δ_to_x.(u_δ18O_SWATI[1:(NLAYER-1)], R_VSMOW¹⁸O)]
        qCᵢ¹⁸O_low  =     aux_du_VRFLI[1:(NLAYER)]    .*     δ_to_x.(u_δ18O_SWATI[1:(NLAYER)], R_VSMOW¹⁸O)
        qCᵢ²H_upp   = [0; aux_du_VRFLI[1:(NLAYER-1)]] .* [0; δ_to_x.(u_δ2H_SWATI[1:(NLAYER-1)],  R_VSMOW²H )]
        qCᵢ²H_low   =     aux_du_VRFLI[1:(NLAYER)]    .*     δ_to_x.(u_δ2H_SWATI[1:(NLAYER)],  R_VSMOW²H )

        Cᵢ¹⁸O_INFLI = δ_to_x.(p_δ18O_PREC, R_VSMOW¹⁸O)  # TODO(bernhard): for debugging, remove this again and replace with δ18O_INFLI
        # Cᵢ¹⁸O_INFLI = δ18O_INFLI
        # Cᵢ¹⁸O_INFLI = ifelse(sum(aux_du_INFLI) == 0, 0, δ18O_INFLI) # in case there is no inflow δ18O_INFLI was set to NaN, set it to zero for below equation
        Cᵢ¹⁸O_TRANI = Cᵢ¹⁸O # no fractionation occurring, i.e. outflux composition equal to storage composition
        Cᵢ¹⁸O_DSFL  = Cᵢ¹⁸O # no fractionation occurring, i.e. outflux composition equal to storage composition
        Cᵢ¹⁸O_SLVP  = [δ_to_x.(δ¹⁸O_SLVP, R_VSMOW¹⁸O); fill(0, NLAYER-1)]
        Cᵢ²H_INFLI = δ_to_x.(p_δ2H_PREC,  R_VSMOW²H )   # TODO(bernhard): for debugging, remove this again and replace with δ2H_INFLI
        # Cᵢ²H_INFLI = δ2H_INFLI
        # Cᵢ²H_INFLI = ifelse(sum(aux_du_INFLI) == 0, 0, δ2H_INFLI) # in case there is no inflow δ2H_INFLI was set to NaN, set it to zero for below equation
        Cᵢ²H_TRANI = Cᵢ²H # no fractionation occurring, i.e. outflux composition equal to storage composition
        Cᵢ²H_DSFL  = Cᵢ²H # no fractionation occurring, i.e. outflux composition equal to storage composition
        Cᵢ²H_SLVP  = [δ_to_x.(δ²H_SLVP,  R_VSMOW²H );  fill(0, NLAYER-1)]
    end
    ##################
    ##################
    ##################


    ##################
    ##################
    ##################
    if solution_type == "numerical-du"
        # 3) Some other terms in the isotope balance equation
        dVdt = du_NTFLI # [mm/day]
        # θ = u_aux_θ
        # dθdt = aux_du_VRFLI[NLAYER] - du_GWFL - du_SEEP # dV/dt [-/day] # TODO(bernhard)

                        # du_NTFL was computed in MSBITERATE as following:
                                    #                      = ∫∂/∂z[K_h(∂h/∂z + 1)] dz   - ∫S dz
                                    # for i=1:    NTFLI[i] = 0          - VRFLI[i]      + INFLI[i] - aux_du_TRANI[i] - aux_du_DSFLI[i] - aux_du_SLVP
                                    # for i=x:    NTFLI[i] = VRFLI[i-1] - VRFLI[i]      + INFLI[i] - aux_du_TRANI[i] - aux_du_DSFLI[i]
    end
    ##################
    ##################
    ##################


    ##################
    ##################
    ##################
    if solution_type == "numerical-du"
        # 4a) Isotope balance equation:

        # SOIL WATER
        ### THEORY:
        ### ∂/∂t Cᵢ     = [- Cᵢ/SWATI * ∂/∂t SWATI]           + 1/SWATI *
        ###                      [diff(z_upper) - diff(z_lower) - qCᵢ(z_upper) + qCᵢ(z_lower) +
        ###                       INFLI*Cᵢ_{INFLI} - TRANI*Cᵢ_{TRANI} - DSFL*Cᵢ_{DSFL} - SLVP*Cᵢ_{SLVP}]

        du_Cᵢ¹⁸_SWATI = -Cᵢ¹⁸O./u_SWATI .* dVdt .+ 1 ./ u_SWATI .* (
                                - diff_upp .+ diff_low .+ qCᵢ¹⁸O_upp .- qCᵢ¹⁸O_low .+ # TODO(bernhard): check sign for qC
                                aux_du_INFLI.*Cᵢ¹⁸O_INFLI .- aux_du_TRANI.*Cᵢ¹⁸O_TRANI .- aux_du_DSFLI.*Cᵢ¹⁸O_DSFL .- aux_du_SLVP.*Cᵢ¹⁸O_SLVP
                            )
        du_Cᵢ²H_SWATI = -Cᵢ²H./u_SWATI .* dVdt .+ 1 ./ u_SWATI .* (
                                -diff_upp .+ diff_low .+ qCᵢ²H_upp .- qCᵢ²H_low .+ # TODO(bernhard): check sign for qC
                                aux_du_INFLI.*Cᵢ²H_INFLI .- aux_du_TRANI.*Cᵢ²H_TRANI .- aux_du_DSFLI.*Cᵢ²H_DSFL .- aux_du_SLVP.*Cᵢ²H_SLVP
                            )
        # # NOTE: below max(0.001,u_SWATI) makes the code more robust
        # du_Cᵢ¹⁸_SWATI = -Cᵢ¹⁸O./max.(0.001,u_SWATI) .* dVdt .+ 1 ./ max.(0.001,u_SWATI) .* (
        #                         diff_upp .- diff_low .+ qCᵢ¹⁸O_upp .- qCᵢ¹⁸O_low .+ # TODO(bernhard): check sign for qC
        #                         aux_du_INFLI.*Cᵢ¹⁸O_INFLI .- aux_du_TRANI.*Cᵢ¹⁸O_TRANI .- aux_du_DSFLI.*Cᵢ¹⁸O_DSFL .- aux_du_SLVP.*Cᵢ¹⁸O_SLVP
        #                     )
        # du_Cᵢ²H_SWATI = -Cᵢ²H./max.(0.001,u_SWATI) .* dVdt .+ 1 ./ max.(0.001,u_SWATI) .* (
        #                         diff_upp .- diff_low .+ qCᵢ²H_upp .- qCᵢ²H_low .+ # TODO(bernhard): check sign for qC
        #                         aux_du_INFLI.*Cᵢ²H_INFLI .- aux_du_TRANI.*Cᵢ²H_TRANI .- aux_du_DSFLI.*Cᵢ²H_DSFL .- aux_du_SLVP.*Cᵢ²H_SLVP
        #                     )

        # GROUND WATER (GWAT)
        δ18O_empty = NaN
        δ2H_empty  = NaN
        @assert aux_du_VRFLI[NLAYER] >= 0 "aux_du_VRFLI[NLAYER] should not be negative"

        if ((u_GWAT == 0) & (aux_du_VRFLI[NLAYER] == 0)) # initially no groundwater and no new is added
            du_δ18O_GWAT = δ18O_empty
            du_δ2H_GWAT  = δ2H_empty
        elseif ((u_GWAT == 0) & (aux_du_VRFLI[NLAYER] > 0)) # initially no groundwater but some is added
            # If no groundwater initially (u_δ_GWAT were NaN)
            # In that case override u to fixed value and set du to zero
            #TODO(bernhard): this is not working in f, we'd need to do it in the callback to set u: idx_u_scalar_isotopes_d18O = integrator.p[1][4][8  ]
            #TODO(bernhard): this is not working in f, we'd need to do it in the callback to set u: idx_u_scalar_isotopes_d2H  = integrator.p[1][4][10  ]
            #TODO(bernhard): this is not working in f, we'd need to do it in the callback to set u: integrator.u[idx_u_scalar_isotopes_d18O[1]] = u_δ18O_SWATI[end] # u_δ18O_GWAT
            #TODO(bernhard): this is not working in f, we'd need to do it in the callback to set u: integrator.u[idx_u_scalar_isotopes_d2H[1] ] = u_δ2H_SWATI[end]  # u_δ2H_GWAT

            du_δ18O_GWAT = 0
            du_δ2H_GWAT  = 0
        else
            # If there is groundwater initially set du to required value
            # composition of fluxes
            δ18O_in_GWAT  = u_δ18O_SWATI[end]
            δ2H_in_GWAT   = u_δ2H_SWATI[end]
            # δ18O_out_GWAT = u_δ18O_GWAT # not needed, because it does not change composition
            # δ2H_out_GWAT  = u_δ2H_GWAT  # not needed, because it does not change composition

            # mass balance:
            # dVdt = aux_du_VRFLI[NLAYER] - du_GWFL - du_SEEP # dV/dt [mm/day]
            V_GWAT = u_GWAT
            in_GWAT = aux_du_VRFLI[NLAYER]
            # out_GWAT = (du_GWFL + du_SEEP)

            # isotope balance
            Cᵢ¹⁸O_GWAT = δ_to_x.(u_δ18O_GWAT,R_VSMOW¹⁸O)
            Cᵢ²H_GWAT  = δ_to_x.(u_δ2H_GWAT,R_VSMOW²H )
            du_Cᵢ¹⁸_GWAT = in_GWAT/V_GWAT*(δ_to_x.(δ18O_in_GWAT, R_VSMOW¹⁸O) -  Cᵢ¹⁸O_GWAT) # - out_GWAT*(δ18O_out_GWAT - u_δ18O_GWAT) # <- last part == 0
            du_Cᵢ²H_GWAT = in_GWAT/V_GWAT*(δ_to_x.(δ2H_in_GWAT,  R_VSMOW²H ) -  Cᵢ²H_GWAT) # - out_GWAT*(δ2H_out_GWAT  - u_δ2H_GWAT)  # <- last part == 0

            # go back from atom fraction to delta values
            du_δ18O_GWAT  = dxdt_to_dδdt.(du_Cᵢ¹⁸_GWAT, Cᵢ¹⁸O_GWAT, R_VSMOW¹⁸O)
            du_δ2H_GWAT   = dxdt_to_dδdt.(du_Cᵢ²H_GWAT, Cᵢ²H_GWAT, R_VSMOW²H)
        end
        # go back from atom fraction to delta values
        du_δ18O_SWATI = dxdt_to_dδdt.(du_Cᵢ¹⁸_SWATI, Cᵢ¹⁸O, R_VSMOW¹⁸O)
        du_δ2H_SWATI  = dxdt_to_dδdt.(du_Cᵢ²H_SWATI, Cᵢ²H, R_VSMOW²H)
    else
        # 4a) Isotope balance equation: i) mixing, ii) fractionation
        # was all computed above with function `update_δ_with_mixing_and_evaporation`

        u_δ18O_GWAT = u_δ18O_GWAT_final
        u_δ2H_GWAT = u_δ2H_GWAT_final
        u_δ18O_SWATI = u_δ18O_SWATI_final
        u_δ2H_SWATI = u_δ2H_SWATI_final
    end

    ##################
    ##################
    ##################

    if solution_type == "numerical-du"
        return du_δ18O_GWAT, du_δ2H_GWAT, du_δ18O_SWATI, du_δ2H_SWATI
    else
        return u_δ18O_GWAT, u_δ2H_GWAT, u_δ18O_SWATI, u_δ2H_SWATI
    end
end
