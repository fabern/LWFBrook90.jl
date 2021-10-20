# TODO(bernhard): think about where to put function definition of MSBSETVARS(), MSBDAYNIGHT()
"""
    MSBSETVARS()

Compute state dependent parameters for updating states INTS, INTR, SNOW, CC, SNOWLQ in
callback function.

# Arguments
- many
"""
function MSBSETVARS(IDAY, #TODO(bernhard) just for debug... remove again!
                    # arguments
                    FLAG_MualVanGen, NLAYER, p_soil,
                    # for SUNDS
                    p_LAT, p_ESLOPE, DOY, p_L1, p_L2,
                    # for CANOPY
                    p_fT_HEIGHT, p_fT_LAI, p_fT_SAI, u_SNOW, p_SNODEN, p_MXRTLN, p_MXKPL, p_fT_DENSEF,
                    #
                    p_Z0S, p_Z0G,
                    # for ROUGH
                    p_ZMINH, p_CZS, p_CZR, p_HS, p_HR, p_LPC, p_CS,
                    # for PLNTRES
                    p_fT_RELDEN, p_RTRAD, p_FXYLEM,
                    # for WEATHER
                    p_fT_TMAX, p_fT_TMIN, p_fT_EA, p_fT_UW, p_WNDRAT, p_FETCH, p_Z0W, p_ZW, p_fT_SOLRAD,
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
    #
    # solar parameters depending on only on DOY
    # TODO(bernhard): a) Do this outside of integration loop in define_LWFB90_p() p_fT_DAYLEN
    p_fT_DAYLEN, p_fT_I0HDAY, p_fT_SLFDAY = LWFBrook90.SUN.SUNDS(p_LAT, p_ESLOPE, DOY, p_L1, p_L2, LWFBrook90.CONSTANTS.p_SC, LWFBrook90.CONSTANTS.p_PI, LWFBrook90.CONSTANTS.p_WTOMJ)

    # canopy parameters depending on DOY as well as different state parameters of snow depth
    p_fu_HEIGHTeff, p_fu_LAIeff, p_fT_SAIeff, p_fT_RTLEN, p_fT_RPLANT =
        LWFBrook90.PET.LWFBrook90_CANOPY(p_fT_HEIGHT,
                          p_fT_LAI,  # leaf area index, m2/m2, minimum of 0.00001
                          p_fT_SAI,  # stem area index, m2/m2
                          u_SNOW,    # water equivalent of snow on the ground, mm SWE
                          p_SNODEN,  # snow density, mm SWE/mm depth
                          p_MXRTLN,  # maximum root length per unit land area, m/m2
                          p_MXKPL,   # maximum plant conductivity, (mm/d)/MPa
                          p_fT_DENSEF)

    # roughness parameters depending on u_SNOW
    if (u_SNOW > 0)
        p_fu_Z0GS = p_Z0S
    else
        p_fu_Z0GS = p_Z0G
    end
    p_fu_Z0GS, p_fu_Z0C, p_fu_DISPC, p_fu_Z0, p_fu_DISP, p_fu_ZA =
            LWFBrook90.PET.ROUGH(p_fu_HEIGHTeff, p_ZMINH, p_fu_LAIeff, p_fT_SAIeff,
                                      p_CZS, p_CZR, p_HS, p_HR, p_LPC, p_CS, p_fu_Z0GS)

    # plant resistance components
    p_fT_RXYLEM, p_fT_RROOTI, p_fT_ALPHA = LWFBrook90.EVP.PLNTRES(NLAYER, p_soil, p_fT_RTLEN, p_fT_RELDEN, p_RTRAD, p_fT_RPLANT, p_FXYLEM, LWFBrook90.CONSTANTS.p_PI, LWFBrook90.CONSTANTS.p_RHOWG)

    # calculated weather data
    p_fu_SHEAT = 0.
    (p_fu_SOLRADC, p_fu_TA, p_fu_TADTM, p_fu_TANTM, UA, p_fu_UADTM, p_fu_UANTM) =
        LWFBrook90.PET.WEATHER(p_fT_TMAX, p_fT_TMIN, p_fT_DAYLEN, p_fT_I0HDAY, p_fT_EA, p_fT_UW, p_fu_ZA, p_fu_DISP, p_fu_Z0, p_WNDRAT, p_FETCH, p_Z0W, p_ZW, p_fT_SOLRAD)
    # fraction of precipitation as p_fT_SFAL
    p_fT_SNOFRC= LWFBrook90.SNO.SNOFRAC(p_fT_TMAX, p_fT_TMIN, p_RSTEMP)

    if (u_SNOW > 0)
        # snowpack temperature at beginning of day
        p_fu_TSNOW = -u_CC / (p_CVICE * u_SNOW)
        # potential snow evaporation PSNVP
        p_fu_PSNVP=LWFBrook90.SNO.SNOVAP(p_fu_TSNOW, p_fu_TA, p_fT_EA, UA, p_fu_ZA, p_fu_HEIGHTeff, p_fu_Z0, p_fu_DISP, p_fu_Z0C, p_fu_DISPC, p_fu_Z0GS, p_LWIDTH, p_RHOTP, p_NN, p_fu_LAIeff, p_fT_SAIeff, p_KSNVP)
        p_fu_ALBEDO = p_ALBSN
        p_fu_RSS = 0.
    else
        p_fu_TSNOW = 0.
        p_fu_PSNVP = 0.
        p_fu_ALBEDO = p_ALB
        # soil evaporation resistance
        p_fu_RSS = LWFBrook90.PET.FRSS(p_RSSA, p_RSSB, u_aux_PSIM[1], p_soil)

        # check for zero or negative p_fu_RSS (TODO: not done in LWFBrook90)
        #if (p_fu_RSS < 0.000001)
        #    error("p_fu_RSS is very small or negative. Run ends. Check p_RSSA and p_RSSB values.")
        #end
    end

    # snow surface energy balance (is performed even when SNOW=0 in case snow is added during day)
    p_fu_SNOEN = LWFBrook90.SNO.SNOENRGY(p_fu_TSNOW, p_fu_TA, p_fT_DAYLEN, p_CCFAC, p_MELFAC, p_fT_SLFDAY, p_fu_LAIeff, p_fT_SAIeff, p_LAIMLT, p_SAIMLT)

    return (p_fT_DAYLEN, p_fT_I0HDAY, p_fT_SLFDAY,
            p_fu_HEIGHTeff, p_fu_LAIeff, p_fT_SAIeff, p_fT_RTLEN, p_fT_RPLANT,
            p_fu_Z0GS, p_fu_Z0C, p_fu_DISPC, p_fu_Z0, p_fu_DISP, p_fu_ZA,
            p_fT_RXYLEM, p_fT_RROOTI, p_fT_ALPHA,
            p_fu_SHEAT,
            p_fu_SOLRADC, p_fu_TA, p_fu_TADTM, p_fu_TANTM, p_fu_UADTM, p_fu_UANTM,
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

http://www.ecoshift.net/brook/pet.html

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
"""
function MSBDAYNIGHT(IDAY, #TODO(bernhard) just for debug... remove again!
                     FLAG_MualVanGen,
                     # arguments
                     p_fT_SLFDAY, p_fu_SOLRADC, p_WTOMJ, p_fT_DAYLEN, p_fu_TADTM, p_fu_UADTM, p_fu_TANTM, p_fu_UANTM,
                     p_fT_I0HDAY,
                     # for AVAILEN:
                     p_fu_ALBEDO, p_C1, p_C2, p_C3, p_fT_EA, p_fu_SHEAT, p_CR, p_fu_LAIeff, p_fT_SAIeff,
                     # for SWGRA:
                     p_fu_ZA, p_fu_HEIGHTeff, p_fu_Z0, p_fu_DISP, p_fu_Z0C, p_fu_DISPC, p_fu_Z0GS, p_LWIDTH, p_RHOTP, p_NN,
                     # for SRSC:
                     p_fu_TA, p_GLMIN, p_GLMAX, p_R5, p_CVPD, p_RM, p_TL, p_T1, p_T2, p_TH,
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

    if FLAG_MualVanGen == 1
        p_fu_PGER = fill(NaN, 2) # fill in values further down
    else
        p_fu_PGER = fill(NaN, 2) # don't fill in values, simply return NaN
    end

    ATR = fill(NaN, 2)
    SLRAD=fill(NaN,2)
    for J = 1:2 # 1 for daytime, 2 for nighttime

        # net radiation
        if (J ==1)
            SLRAD[J] = p_fT_SLFDAY * p_fu_SOLRADC / (p_WTOMJ * p_fT_DAYLEN)
            TAJ = p_fu_TADTM
            UAJ = p_fu_UADTM
        else
            SLRAD[J] = 0
            TAJ = p_fu_TANTM
            UAJ = p_fu_UANTM
        end

        # if (p_fT_I0HDAY <= 0.01)
        #     # TODO(bernhard): Brook90 did treat this case specially, LWFBrook90 did not
        #
        #     # no sunrise, assume 50% clouds for longwave
        #     cloud_fraction = 0.5
        # else
        cloud_fraction = p_fu_SOLRADC / p_fT_I0HDAY
        # end
        # Available energy above canopy (AA) and at soil/substrate (ASUBS)
        AA, ASUBS =
            LWFBrook90.SUN.AVAILEN(SLRAD[J], p_fu_ALBEDO, p_C1, p_C2, p_C3, TAJ, p_fT_EA,
                    cloud_fraction,
                    p_fu_SHEAT, p_CR, p_fu_LAIeff, p_fT_SAIeff)

        # vapor pressure deficit
        ES, DELTA = LWFBrook90.PET.ESAT(TAJ)
        VPD = ES - p_fT_EA
        # S.-W. resistances
        RAA, RAC, RAS = LWFBrook90.PET.SWGRA(UAJ, p_fu_ZA, p_fu_HEIGHTeff, p_fu_Z0, p_fu_DISP, p_fu_Z0C, p_fu_DISPC, p_fu_Z0GS, p_LWIDTH, p_RHOTP, p_NN, p_fu_LAIeff, p_fT_SAIeff)
        if (J == 1)
            RSC=LWFBrook90.PET.SRSC(SLRAD[J], p_fu_TA, VPD, p_fu_LAIeff, p_fT_SAIeff, p_GLMIN, p_GLMAX, p_R5, p_CVPD, p_RM, p_CR, p_TL, p_T1, p_T2, p_TH)
        else
            RSC = 1 / (p_GLMIN * p_fu_LAIeff)
        end
        # RSC: canopy surface resistance

        #print("\nIDAY:$(@sprintf("% 3d", IDAY)), J = $J     AA:$(@sprintf("% 8.4f", AA)), ASUBS:$(@sprintf("% 8.4f", ASUBS)), VPD:$(@sprintf("% 8.4f", VPD)), RAA:$(@sprintf("% 8.4f", RAA)), RAC:$(@sprintf("% 8.4f", RAC)), RAS:$(@sprintf("% 8.4f", RAS)), RSC:$(@sprintf("% 8.4f", RSC)), p_fu_RSS:$(@sprintf("% 8.4f", p_fu_RSS)), DELTA:$(@sprintf("% 8.4f", DELTA))")

        # S.-W. potential transpiration and ground evaporation rates
        p_fu_PTR[J], p_fu_GER[J] =  LWFBrook90.PET.SWPE(AA, ASUBS, VPD, RAA, RAC, RAS, RSC, p_fu_RSS, DELTA)
        # S.-W. potential interception and ground evap. rates
        # RSC = 0, p_fu_RSS not changed
        p_fu_PIR[J], p_fu_GIR[J] =  LWFBrook90.PET.SWPE(AA, ASUBS, VPD, RAA, RAC, RAS, 0, p_fu_RSS, DELTA)

        if FLAG_MualVanGen == 1
            # S.-W. potential interception and ground evap. rates
            # RSC not changed, p_fu_RSS = 0
            _, p_fu_PGER[J] =  LWFBrook90.PET.SWPE(AA, ASUBS, VPD, RAA, RAC, RAS, RSC, 0, DELTA)
        end
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
            p_fu_ATRI,# actual transp.rate from layer for daytime and night (mm/d)
            p_fu_PGER)

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
function MSBDAYNIGHT_postprocess(FLAG_MualVanGen, NLAYER,
                                 p_fu_PTR, # potential transpiration rate for daytime or night (mm/d)
                                 p_fu_GER, # ground evaporation rate for daytime or night (mm/d)
                                 p_fu_PIR, # potential interception rate for daytime or night (mm/d)
                                 p_fu_GIR, # ground evap. rate with intercep. for daytime or night (mm/d)
                                 p_fu_ATRI,# actual transp.rate from layer for daytime and night (mm/d)
                                 p_fT_DAYLEN,
                                 p_fu_PGER)

    # average rates over day (mm/day)
    # mm/day   =  mm/day      * day_fraction+ mm/day      *  day_fraction
    p_fu_PTRAN = (p_fu_PTR[1] * p_fT_DAYLEN + p_fu_PTR[2] * (1 - p_fT_DAYLEN))
    p_fu_GEVP  = (p_fu_GER[1] * p_fT_DAYLEN + p_fu_GER[2] * (1 - p_fT_DAYLEN))
    p_fu_PINT  = (p_fu_PIR[1] * p_fT_DAYLEN + p_fu_PIR[2] * (1 - p_fT_DAYLEN))
    p_fu_GIVP  = (p_fu_GIR[1] * p_fT_DAYLEN + p_fu_GIR[2] * (1 - p_fT_DAYLEN))

    if FLAG_MualVanGen==1
        p_fu_PSLVP = (p_fu_PGER[1] * p_fT_DAYLEN + p_fu_PGER[2] * (1 - p_fT_DAYLEN))
    end

    aux_du_TRANI=zeros(NLAYER)

    for i = 1:NLAYER
        aux_du_TRANI[i] = (p_fu_ATRI[1, i] * p_fT_DAYLEN + p_fu_ATRI[2, i] * (1 - p_fT_DAYLEN))
    end

    return (p_fu_PTRAN, # average potential transpiration rate for day (mm/d)
            p_fu_GEVP,  # average ground evaporation for day (mm/d)
            p_fu_PINT,  # average potential interception for day (mm/d)
            p_fu_GIVP,  # average ground evaporation for day with interception (mm/d)
            p_fu_PSLVP, # average potential evaporation rate from soil for day (mm/d) TODO(bernhard) seems unused further on
            aux_du_TRANI) # average transpiration rate for day from layer (mm/d)
end


function MSBPREINT(#arguments:
                   p_fT_PREC, p_DTP, p_fT_SNOFRC, p_NPINT, p_fu_PINT, p_fu_TA,
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

    if (p_NPINT > 1.0)
        # more than one precip interval in day
        error("Case with multiple precipitation intervals (using PRECDAT and precip_interval != 1) is not implemented.")
        # # snow interception
        # if (p_fu_PINT < 0 && p_fu_TA > 0)
        #     # prevent frost when too warm, carry negative p_fu_PINT to rain
        #     aux_du_SINT, aux_du_ISVP = LWFBrook90.EVP.INTER(p_fT_SFAL, 0, p_fu_LAI, p_fu_SAI, p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS, p_DTP, u_INTS)
        # else
        #     aux_du_SINT, aux_du_ISVP = LWFBrook90.EVP.INTER(p_fT_SFAL, p_fu_PINT, p_fu_LAI, p_fu_SAI, p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS, p_DTP, u_INTS)
        # end
        # # rain interception,  note potential interception rate is PID-aux_du_ISVP (mm/day)
        # aux_du_RINT, aux_du_IRVP = LWFBrook90.EVP.INTER(p_fT_RFAL, p_fu_PINT - aux_du_ISVP, p_fu_LAI, p_fu_SAI, p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS, p_DTP, u_INTR)
    else
        # one precip interval in day, use storm p_DURATN and INTER24
        # snow interception
        if (p_fu_PINT < 0 && p_fu_TA > 0)
            # prevent frost when too warm, carry negative p_fu_PINT to rain
            aux_du_SINT, aux_du_ISVP = LWFBrook90.EVP.INTER24(p_fT_SFAL, 0, p_fu_LAI, p_fu_SAI, p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS, p_DURATN, u_INTS, MONTHN)
        else
            aux_du_SINT, aux_du_ISVP = LWFBrook90.EVP.INTER24(p_fT_SFAL, p_fu_PINT, p_fu_LAI, p_fu_SAI, p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS, p_DURATN, u_INTS, MONTHN)
        end
        # rain interception,  note potential interception rate is PID-aux_du_ISVP (mm/day)
        aux_du_RINT, aux_du_IRVP = LWFBrook90.EVP.INTER24(p_fT_RFAL, p_fu_PINT - aux_du_ISVP, p_fu_LAI, p_fu_SAI, p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS, p_DURATN, u_INTR, MONTHN)
    end

    # throughfall
    p_fu_RTHR = p_fT_RFAL - aux_du_RINT
    p_fu_STHR = p_fT_SFAL - aux_du_SINT

    # reduce transpiration for fraction of precip interval that canopy is wet
    p_fu_WETFR = min(1.0, (aux_du_IRVP + aux_du_ISVP) / p_fu_PINT)
    p_fu_PTRAN = (1.0 - p_fu_WETFR) * p_fu_PTRAN
    for i = 1:NLAYER
        aux_du_TRANI[i] = (1.0 - p_fu_WETFR) * aux_du_TRANI[i]
    end
    if (u_SNOW <= 0 && p_fu_STHR <= 0)
        # no snow, soil evaporation weighted for p_fu_WETFR
        aux_du_SLVP = p_fu_WETFR * p_fu_GIVP + (1.0 - p_fu_WETFR) * p_fu_GEVP
        p_fu_RNET = p_fu_RTHR
        aux_du_RSNO = 0.0
        aux_du_SNVP = 0.0
        aux_du_SMLT = 0.0
    else
        if (u_SNOW <= 0 && p_fu_STHR > 0)
            # new snow only, zero CC and SNOWLQ assumed
            u_CC = 0.0
            u_SNOWLQ = 0.0
        end
        # snow accumulation and melt
        u_CC, u_SNOW, u_SNOWLQ, aux_du_RSNO, aux_du_SNVP, aux_du_SMLT =
          LWFBrook90.SNO.SNOWPACK(p_fu_RTHR, p_fu_STHR, p_fu_PSNVP, p_fu_SNOEN,
                   # States that are overwritten:
                   u_CC, u_SNOW, u_SNOWLQ,
                   p_DTP, p_fu_TA, p_MAXLQF, p_GRDMLT,
                   p_CVICE, p_LF, p_CVLQ)

        p_fu_RNET = p_fu_RTHR - aux_du_RSNO
        aux_du_SLVP = 0.0
    end

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
                    p_DSWMAX, u_aux_θ,
                    # for GWATER:
                    u_GWAT, p_GSC, p_GSP)

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

    # water supply rate to soil surface:
    p_fu_SLFL = p_fu_RNET + aux_du_SMLT - p_fu_SRFL

    ## Within soil compute flows from layers:
    #  a) vertical flux between layers VRFLI,
    #  b) downslope flow from the layers DSFLI
    aux_du_VRFLI = fill(NaN, NLAYER) # 0.000001 seconds (1 allocation: 144 bytes)
    aux_du_DSFLI = fill(NaN, NLAYER) # 0.000001 seconds (1 allocation: 144 bytes)

    # downslope flow rates
    if (p_LENGTH_SLOPE == 0 || p_DSLOPE == 0)
        # added in Version 4
        aux_du_DSFLI .= 0
    else
        aux_du_DSFLI .= LWFBrook90.WAT.DSLOP(p_DSLOPE, p_LENGTH_SLOPE, p_RHOWG, p_soil, u_aux_PSIM, p_fu_KK)
    end

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
        # if (IDAY >= 6 && i==NLAYER) # TODO(bernhard): this seemed like a no effect snippet
        #     p_DRAIN=p_DRAIN         # TODO(bernhard): this seemed like a no effect snippet
        # end                         # TODO(bernhard): this seemed like a no effect snippet
    end

    # first approximation on aux_du_VRFLI
    aux_du_VRFLI_1st_approx = aux_du_VRFLI
    #@debug "u_aux_PSITI[1]: $(u_aux_PSITI[1]), sum(u_aux_PSITI): $(sum(u_aux_PSITI))"
    #@debug "a) aux_du_VRFLI_1st_approx[1]: $(aux_du_VRFLI_1st_approx[1]), sum(aux_du_VRFLI_1st_approx): $(sum(aux_du_VRFLI_1st_approx))"

    # first approximation for iteration time step,time remaining or DTIMAX
    p_DTP    # first approximation for iteration time step,time remaining or DTIMAX
    DTRI = p_DTP - (t % p_DTP) # Time remaining
    # DTRI = 1.0 - (t % 1.0)   # as p_DTP is 1.0 days in a default simulation
    DTI = min(DTRI, p_DTIMAX)

    # net inflow to each layer including E and T withdrawal adjusted for interception
    aux_du_VRFLI, aux_du_INFLI, aux_du_BYFLI, du_NTFLI =
        LWFBrook90.WAT.INFLOW(NLAYER, DTI, p_INFRAC, p_fu_BYFRAC, p_fu_SLFL, aux_du_DSFLI, aux_du_TRANI,
                                    aux_du_SLVP, p_soil.p_SWATMAX, u_SWATI,
                                    aux_du_VRFLI_1st_approx) # 0.000003 seconds (4 allocations: 576 bytes)
    #@debug "b) aux_du_VRFLI[1]: $(aux_du_VRFLI[1]), sum(aux_du_VRFLI): $(sum(aux_du_VRFLI))"

    DPSIDW = LWFBrook90.KPT.FDPSIDWF(u_aux_WETNES, p_soil) # 0.000004 seconds (3 allocations: 432 bytes)

    # limit step size
    #       TODO(bernhard): This could alternatively be achieved with a stepsize limiter callback
    #                       see: https://diffeq.sciml.ai/stable/features/callback_library/#Stepsize-Limiters
    #       TODO(bernhard): Or it could alternatively be set manually using set_proposed_dt!
    #                       see: https://diffeq.sciml.ai/stable/basics/integrator/#SciMLBase.set_proposed_dt!
    #   ITER computes DTI so that the potential difference (due to aux_du_VRFLI)
    #   between adjacent layers does not change sign during the iteration time step
    DTINEW=LWFBrook90.WAT.ITER(NLAYER, FLAG_MualVanGen, DTI, LWFBrook90.CONSTANTS.p_DTIMIN, DPSIDW,
                                    du_NTFLI, u_aux_PSITI, u_aux_θ, p_DSWMAX, p_DPSIMAX,
                    p_soil) # 0.000002 seconds (2 allocations: 288 bytes)
    # recompute step
    if (DTINEW < DTI)
        # recalculate flow rates with new DTI
        DTI = DTINEW
        aux_du_VRFLI, aux_du_INFLI, aux_du_BYFLI, du_NTFLI =
            LWFBrook90.WAT.INFLOW(NLAYER, DTI, p_INFRAC, p_fu_BYFRAC, p_fu_SLFL, aux_du_DSFLI, aux_du_TRANI,
                                        aux_du_SLVP, p_soil.p_SWATMAX, u_SWATI,
                                        aux_du_VRFLI_1st_approx)
    end
    #@debug "c) aux_du_VRFLI[1]: $(aux_du_VRFLI[1]), sum(aux_du_VRFLI): $(sum(aux_du_VRFLI))"

    # groundwater flow and seepage loss
    du_GWFL, du_SEEP = LWFBrook90.WAT.GWATER(u_GWAT, p_GSC, p_GSP, aux_du_VRFLI[NLAYER])


    return (p_fu_SRFL, p_fu_SLFL, aux_du_DSFLI, aux_du_VRFLI, DTI, aux_du_INFLI, aux_du_BYFLI, du_NTFLI, du_GWFL, du_SEEP)
end


"""
    compute_isotope_INTS_INTR_SNOW!(
        p_δ2H_PREC, p_δ18O_PREC, p_fu_TADTM, p_EA,
        # for INTS (in: SINT; out: ISVP):
        u_INTS, aux_du_SINT, aux_du_ISVP, p_DTP, u_δ2H_INTS, u_δ18O_INTS,
        # for INTR (in: RINT; out: IRVP):
        u_INTR, aux_du_RINT, aux_du_IRVP, u_δ2H_INTR, u_δ18O_INTR,
        # for SNOW (in: STHR, RSNO (both δ_PREC); out: SMLT, SNVP (δ_SNOW and fractionated)):
        p_fu_STHR, aux_du_RSNO, aux_du_SMLT, aux_du_SNVP, u_δ2H_SNOW, u_δ18O_SNOW,
        # to compute isotopic signature of soil infiltration: SLFL
        p_fu_RNET)

In-line update of amounts of states INTS, INTR, and SNOW as well as their isotopic composition.
Compute mixing and evaporative fractionation, and also compute the isotopic composotion of
the resulting flux that infiltrates into the soil: δ18O_SLFL, δ2H_SLFL

The function is called in the daily callback.

The function returns (δ18O_SLFL, δ2H_SLFL) but also modifies the following input arguments in-place:
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
function compute_isotope_INTS_INTR_SNOW!(
    p_δ2H_PREC, p_δ18O_PREC, p_fu_TADTM, p_EA,
    # for INTS (in: SINT; out: ISVP):
    u_INTS, aux_du_SINT, aux_du_ISVP, p_DTP, u_δ2H_INTS, u_δ18O_INTS,
    # for INTR (in: RINT; out: IRVP):
    u_INTR, aux_du_RINT, aux_du_IRVP, u_δ2H_INTR, u_δ18O_INTR,
    # for SNOW (in: STHR, RSNO (both δ_PREC); out: SMLT, SNVP (δ_SNOW and fractionated)):
    u_SNOW, u_SNOW_new, p_fu_STHR, aux_du_RSNO, aux_du_SMLT, aux_du_SNVP, u_δ2H_SNOW, u_δ18O_SNOW,
    # to compute isotopic signature of soil infiltration: SLFL
    p_fu_RNET)
    # check chart "../docs/src/assets/b90flow.gif"

    ##################
    ##################
    ##################
    # Define conditions for isotope calculations:
    Tc = p_fu_TADTM  # °C, average daytime air temperature
    h = min(1.0,p_EA) # -, relative humidity of the atmosphere (vappress_atm/1 atm)
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
    # 2a) INTS (in: SINT*δ_SINT; out: ISVP*δ_ISVP)
    #          with δ_SINT = δ_PREC; δ_ISVP = f(f, α, ...)

    # Alternative formulation:
    # ε_eq2H  = (LWFBrook90.ISO.α¹⁸O_eq(Tc)-1)*1000
    # ε_eq18O = (LWFBrook90.ISO.α¹⁸O_eq(Tc)-1)*1000
    # ε_kin2H  = (LWFBrook90.ISO.α²H_dif - 1)*1000
    # ε_kin18O = (LWFBrook90.ISO.α¹⁸O_dif - 1)*1000
    # δ2H_E  = ((δ2H  - ε_eq2H) /LWFBrook90.ISO.α²H_eq(Tc)  - h*δ²H_a  - ε_kin2H)  / (1. - h + ε_kin2H/1000) # [‰] Benettin 2018 HESS eq. 1 and Gibson 2016 QSR eq. 3
    # δ18O_E = ((δ18O - ε_eq18O)/LWFBrook90.ISO.α¹⁸O_eq(Tc) - h*δ¹⁸O_a - ε_kin18O) / (1. - h + ε_kin18O/1000)# [‰] Benettin 2018 HESS eq. 1 and Gibson 2016 QSR eq. 3
    δ18O_empty = NaN
    δ2H_empty  = NaN

    # Operator step 1
    @assert aux_du_SINT >= 0 "aux_du_SINT should not be negative"
    if ((u_INTS == 0) & (aux_du_SINT == 0)) # initially no intercepted snow and no new is added
        u_δ18O_INTS_final = δ18O_empty
        u_δ2H_INTS_final  = δ2H_empty
        u_INTS_final      = 0
    else
        u_INTS_first = u_INTS + aux_du_SINT*p_DTP
        if ((u_INTS == 0) & (aux_du_SINT > 0)) # initially no intercepted snow but some is added
            u_δ18O_INTS_first = p_δ18O_PREC
            u_δ2H_INTS_first  = p_δ2H_PREC
        else # u_INTS is not zero, and some/none is added (aux_du_SINT=0 or aux_du_SINT>0)
            # some bug in this formulation: u_δ18O_INTS_first = u_δ18O_INTS + aux_du_SINT*p_DTP/u_INTS * (p_δ18O_PREC-u_δ18O_INTS)
            # some bug in this formulation: u_δ2H_INTS_first  = u_δ2H_INTS  + aux_du_SINT*p_DTP/u_INTS * (p_δ2H_PREC -u_δ2H_INTS )
            # u_δ18O_INTS_first = 1/u_INTS_first * (u_δ18O_INTS * u_INTS + p_DTP*(p_δ18O_PREC*aux_du_SINT))
            # u_δ2H_INTS_first  = 1/u_INTS_first * (u_δ2H_INTS  * u_INTS + p_DTP*(p_δ2H_PREC*aux_du_SINT))
            contrib_fraction = min(1,max(0, p_DTP * (aux_du_SINT)/u_INTS_first))
            u_δ18O_INTS_first = (u_δ18O_INTS * (1 - contrib_fraction) + p_δ18O_PREC * contrib_fraction)
            u_δ2H_INTS_first  = (u_δ2H_INTS  * (1 - contrib_fraction) + p_δ2H_PREC  * contrib_fraction)
        end

        # Operator step 2: now taking care of aux_du_ISVP
        u_INTS_final = u_INTS + (aux_du_SINT - aux_du_ISVP) * p_DTP
        f_INTS = min(1,max(0, u_INTS_final/u_INTS_first)) # fraction remaining after evaporation
        # if (f_INTS <= 0.3) #(f_INTS == 0)
        #     u_δ18O_INTS_final = δ18O_empty
        #     u_δ2H_INTS_final  = δ2H_empty
        # else
        f_INTS = min(1.0, max(0.5, f_INTS)) # TODO(bernhard): workaround for stability
            # ε_δ18O = 1/(LWFBrook90.ISO.α¹⁸O_eq(Tc) * LWFBrook90.ISO.α¹⁸O_dif^X_INTS) - 1
            # ε_δ2H  = 1/(LWFBrook90.ISO.α²H_eq(Tc)  * LWFBrook90.ISO.α²H_dif^X_INTS)  - 1
            # u_δ18O_INTS_final = 1000 * ( (1 + u_δ18O_INTS_first/1000)(f_INTS)^ε_δ18O - 1 )
            # u_δ2H_INTS_final  = 1000 * ( (1 + u_δ2H_INTS_first /1000)(f_INTS)^ε_δ2H  - 1 )
            u_δ18O_INTS_final = u_δ18O_INTS_first
            u_δ2H_INTS_final  = u_δ2H_INTS_first
            # TODO(Bernhard): investigate why below leads to instabilities
            # u_δ18O_INTS_final = 1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ18O_INTS_first/1000, δ¹⁸O_a/1000, f_INTS, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_INTS)
            # u_δ2H_INTS_final  = 1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ2H_INTS_first /1000,  δ²H_a/1000, f_INTS, h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_INTS)
        # end
    end

    # 2b) INTR (in: RINT*δ_RINT; out: IRVP*δ_IRVP)
    #          with δ_RINT = δ_PREC; δ_IRVP = f(f, α, ...)
    # Operator step 1

    @assert aux_du_RINT >= 0 "aux_du_RINT should not be negative"
    if ((u_INTR == 0) & (aux_du_RINT == 0)) # initially no intercepted rain and no new is added
        u_δ18O_INTR_final = δ18O_empty
        u_δ2H_INTR_final  = δ2H_empty
        u_INTR_final      = 0
    else
        u_INTR_first = u_INTR + aux_du_RINT*p_DTP
        if ((u_INTR == 0) & (aux_du_RINT > 0)) # initially no intercepted rain but some is added
            u_δ18O_INTR_first = p_δ18O_PREC
            u_δ2H_INTR_first  = p_δ2H_PREC
        else # u_INTR is not zero, and some/none is added (aux_du_RINT=0 or aux_du_RINT>0)
            # some bug in this formulation: u_δ18O_INTR_first = u_δ18O_INTR + aux_du_RINT*p_DTP/u_INTR * (p_δ18O_PREC-u_δ18O_INTR)
            # some bug in this formulation: u_δ2H_INTR_first  = u_δ2H_INTR  + aux_du_RINT*p_DTP/u_INTR * (p_δ2H_PREC -u_δ2H_INTR )
            # u_δ18O_INTR_first = 1/u_INTR_first * (u_δ18O_INTR * u_INTR + p_DTP*(p_δ18O_PREC*aux_du_RINT))
            # u_δ2H_INTR_first  = 1/u_INTR_first * (u_δ2H_INTR  * u_INTR + p_DTP*(p_δ2H_PREC*aux_du_RINT))
            contrib_fraction = min(1,max(0, p_DTP * (aux_du_RINT)/u_INTR_first))
            u_δ18O_INTR_first = (u_δ18O_INTR * (1 - contrib_fraction) + p_δ18O_PREC * contrib_fraction)
            u_δ2H_INTR_first  = (u_δ2H_INTR  * (1 - contrib_fraction) + p_δ2H_PREC  * contrib_fraction)
        end
        # Operator step 2: now taking care of aux_du_IRVP
        u_INTR_final = u_INTR + (aux_du_RINT - aux_du_IRVP) * p_DTP
        f_INTR = max(0, u_INTR_final/u_INTR_first) # fraction remaining after evaporation
        # if (f_INTR <= 0.3) #(f_INTR == 0)
        #     u_δ18O_INTR_final = δ18O_empty
        #     u_δ2H_INTR_final  = δ2H_empty
        # else
        f_INTR = min(1.0, max(0.3, f_INTR)) # TODO(bernhard): workaround for stability
            # ε_δ18O = 1/(LWFBrook90.ISO.α¹⁸O_eq(Tc) * LWFBrook90.ISO.α¹⁸O_dif^X_INTR) - 1
            # ε_δ2H  = 1/(LWFBrook90.ISO.α²H_eq(Tc)  * LWFBrook90.ISO.α²H_dif^X_INTR)  - 1
            # u_δ18O_INTR_final = 1000 * ( (1 + u_δ18O_INTR_first/1000)(f_INTR)^ε_δ18O - 1 )
            # u_δ2H_INTR_final  = 1000 * ( (1 + u_δ2H_INTR_first /1000)(f_INTR)^ε_δ2H  - 1 )
            u_δ18O_INTR_final = u_δ18O_INTR_first
            u_δ2H_INTR_final  = u_δ2H_INTR_first
            # TODO(Bernhard): investigate why below leads to instabilities
            # u_δ18O_INTR_final = 1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ18O_INTR_first/1000, δ¹⁸O_a/1000, f_INTR, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_INTR)
            # u_δ2H_INTR_final  = 1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ2H_INTR_first /1000,  δ²H_a/1000, f_INTR, h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_INTR)
        # end
    end

    # 2c) SNOW (in: STHR*δ_STHR, RSNO*δ_RSNO; out: SMLT*δ_SMLT, SNVP*δ_SNVP)
    #          with δ_STHR = δ_PREC, δ_RSNO = δ_PREC; δ_SMLT = δ_SNOW, δ_SNVP = f(f, α, ...)

    # NOTE: for SNOW the isotope balance is greatly simplified. The most precise
    #       approach would be to define mixing concentrattion within `SNOWPACK()` in `module_SNO.jl`
    # Operator step 1

    @assert p_fu_STHR + aux_du_RSNO >= 0 "p_fu_STHR + aux_du_RSNO should not be negative"
    if ((u_SNOW == 0) & (p_fu_STHR + aux_du_RSNO == 0)) # initially no snowpack and no new is added
        u_δ18O_SNOW_final = δ18O_empty
        u_δ2H_SNOW_final  = δ2H_empty
        u_SNOW_final      = 0
    else
        u_SNOW_first = u_SNOW + p_DTP * (p_fu_STHR + aux_du_RSNO - aux_du_SMLT)
        if ((u_SNOW == 0) & (p_fu_STHR + aux_du_RSNO > 0)) # initially no snowpack but some is added
            u_δ18O_SNOW_first = p_δ18O_PREC
            u_δ2H_SNOW_first  = p_δ2H_PREC
        else # u_SNOW is not zero, and some/none is added ((p_fu_STHR + aux_du_RSNO)=0 or (p_fu_STHR + aux_du_RSNO)>0)
            # p_fu_STHR, aux_du_RSNO, aux_du_SMLT, aux_du_SNVP
            # some bug in this formulation: u_δ18O_SNOW_first = u_δ18O_SNOW + (p_fu_STHR + aux_du_RSNO)*p_DTP/u_SNOW * (p_δ18O_PREC - u_δ18O_SNOW)
            # some bug in this formulation:                             # NOTE: because the outflow term for aux_du_SMLT has
            # some bug in this formulation:                             #       an isotope concentration of u_δ18O_SNOW and it is thus not needed:
            # some bug in this formulation:                             # - (aux_du_SMLT)*p_DTP/u_SNOW * (u_δ18O_SNOW-u_δ18O_SNOW)
            # some bug in this formulation: u_δ2H_SNOW_first = u_δ2H_SNOW + (p_fu_STHR + aux_du_RSNO)*p_DTP/u_SNOW * (p_δ2H_PREC - u_δ2H_SNOW)
            # some bug in this formulation:                             # NOTE: because the outflow term for aux_du_SMLT has
            # some bug in this formulation:                             #       an isotope concentration of u_δ2H_SNOW and it is thus not needed:
            # some bug in this formulation:                             # - (aux_du_SMLT)*p_DTP/u_SNOW * (u_δ2H_SNOW-u_δ2H_SNOW)
            # u_δ18O_SNOW_first = (u_δ18O_SNOW * (u_SNOW - p_DTP*aux_du_SMLT) + p_δ18O_PREC * p_DTP*(p_fu_STHR + aux_du_RSNO)) / u_SNOW_first
            # u_δ2H_SNOW_first  = (u_δ2H_SNOW  * (u_SNOW - p_DTP*aux_du_SMLT) + p_δ2H_PREC  * p_DTP*(p_fu_STHR + aux_du_RSNO)) / u_SNOW_first
            contrib_fraction = min(1,max(0, p_DTP * (p_fu_STHR + aux_du_RSNO)/u_SNOW_first))
            u_δ18O_SNOW_first = (u_δ18O_SNOW * (1 - contrib_fraction) + p_δ18O_PREC * contrib_fraction)
            u_δ2H_SNOW_first  = (u_δ2H_SNOW  * (1 - contrib_fraction) + p_δ2H_PREC  * contrib_fraction)
        end
        # Operator step 2: now taking care of aux_du_SNVP
        u_SNOW_final = u_SNOW + p_DTP * (p_fu_STHR + aux_du_RSNO - aux_du_SMLT - aux_du_SNVP)
        f_SNOW = min(1, max(0, u_SNOW_final/u_SNOW_first)) # fraction remaining after evaporation
        # if (f_SNOW <= 0.3) #(f_SNOW == 0)
        #     u_δ18O_SNOW_final = δ18O_empty
        #     u_δ2H_SNOW_final  = δ2H_empty
        # else
        f_SNOW = max(0.5, f_SNOW) # TODO(bernhard): workaround for stability
            # ε_δ18O = 1/(LWFBrook90.ISO.α¹⁸O_eq(Tc) * LWFBrook90.ISO.α¹⁸O_dif^X_SNOW) - 1
            # ε_δ2H  = 1/(LWFBrook90.ISO.α²H_eq(Tc)  * LWFBrook90.ISO.α²H_dif^X_SNOW)  - 1
            # u_δ18O_SNOW_final = 1000 * ( (1 + u_δ18O_SNOW_first/1000)(f_SNOW)^ε_δ18O - 1 )
            # u_δ2H_SNOW_final  = 1000 * ( (1 + u_δ2H_SNOW_first /1000)(f_SNOW)^ε_δ2H  - 1 )
            u_δ18O_SNOW_final = u_δ18O_SNOW_first
            u_δ2H_SNOW_final  = u_δ2H_SNOW_first
            # TODO(Bernhard): investigate why below leads to instabilities
            # u_δ18O_SNOW_final = 1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ18O_SNOW_first/1000, δ¹⁸O_a/1000, f_SNOW, h, LWFBrook90.ISO.α¹⁸O_eq(Tc), LWFBrook90.ISO.α¹⁸O_dif, γ, X_SNOW)
            # u_δ2H_SNOW_final  = 1000 * LWFBrook90.ISO.δ_CraigGordon.(u_δ2H_SNOW_first /1000,  δ²H_a/1000, f_SNOW, h, LWFBrook90.ISO.α²H_eq(Tc),  LWFBrook90.ISO.α²H_dif,  γ, X_SNOW)
        # end
    end


    # 3) also compute δ_SLFL as mix of δ_SMLT with δ_RNET (i.e. water that infiltrates)
    # if (aux_du_SMLT + p_fu_RNET == 0)
    #     δ18O_SLFL = δ18O_empty
    #     δ2H_SLFL  = δ2H_empty
    # elseif (aux_du_SMLT == 0)
    #     δ18O_SLFL = p_δ18O_PREC
    #     δ2H_SLFL  = p_δ2H_PREC
    # elseif (p_fu_RNET == 0)
    #     δ18O_SLFL = u_δ18O_SNOW_final  # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
    #     δ2H_SLFL  = u_δ2H_SNOW_final  # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
    # else # both fluxes are non-null and we need to compute their mix
    #     δ18O_SLFL = (u_δ18O_SNOW_final * aux_du_SMLT + p_δ18O_PREC * p_fu_RNET) / (aux_du_SMLT + p_fu_RNET)  # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
    #     δ2H_SLFL  = (u_δ2H_SNOW_final * aux_du_SMLT  + p_δ2H_PREC * p_fu_RNET)  / (aux_du_SMLT + p_fu_RNET)  # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
    # end
    # TODO(bernhard): deactivate the following workaround and activate above code
    δ18O_SLFL = p_δ18O_PREC
    δ2H_SLFL  = p_δ2H_PREC

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

    return δ18O_SLFL, δ2H_SLFL, u_INTS_final, u_δ18O_INTS_final, u_δ2H_INTS_final, u_INTR_final, u_δ18O_INTR_final, u_δ2H_INTR_final, u_SNOW_final, u_δ18O_SNOW_final, u_δ2H_SNOW_final
end


function compute_isotope_du_GWAT_SWATI(
    # for GWAT:
    u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT,
    # for SWATI:
    du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI, # (non-fractionating)
    aux_du_SLVP, p_fu_TADTM, p_EA, p_δ2H_PREC, p_δ18O_PREC, u_aux_WETNES, # (fractionating)
    u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, EffectiveDiffusivity_18O, EffectiveDiffusivity_2H,
    )

    ## FROM: 4a Julia Guswa-2002-Water_Resour_Res.pdf
    #    # Isotopic compositions of fluxes.
    #    δ2H_L = δ2H # [‰] no fractionation occurring, i.e. flux composition equal to storage composition
    #    δ2H_T = δ2H # [‰] no fractionation occurring, i.e. flux composition equal to storage composition
    #    δ18O_L = δ18O # [‰] no fractionation occurring, i.e. flux composition equal to storage composition
    #    δ18O_T = δ18O # [‰] no fractionation occurring, i.e. flux composition equal to storage composition
    #    δ2H_E  = ((δ2H  - ε_eq2H(T_cel)) /α_eq2H(T_cel)  - h*δ_A2H  - ε_kin2H(h))  / (1. - h + ε_kin2H(h) /1000)# [‰] Benettin 2018 HESS eq. 1 and Gibson 2016
    #    δ18O_E = ((δ18O - ε_eq18O(T_cel))/α_eq18O(T_cel) - h*δ_A18O - ε_kin18O(h)) / (1. - h + ε_kin18O(h)/1000)# [‰] Benettin 2018 HESS eq. 1 and Gibson 2016 QSR eq. 3
    #    # Rate of change due to balances
    #    dVdt = (-L -T -E) # dV/dt [mm/day]
    #    du[1] = dVdt / (n*Zr) # dS/dt [-/day]
    #    du[2] = - L/V*(δ2H_L - δ2H) - T/V*(δ2H_T - δ2H) - E/V*(δ2H_E - δ2H) # [‰]
    #    du[3] = - L/V*(δ18O_L-δ18O) - T/V*(δ18O_T-δ18O) - E/V*(δ18O_E-δ18O) # [‰]
    #    # equivalently: du[2] = - δ2H /V*dVdt - L/V*δ2H_L - T/V*δ2H_T - E/V*δ2H_E
    #    # equivalently: du[3] = - δ18O/V*dVdt - L/V*δ18O_L - T/V*δ18O_T - E/V*δ18O_E
            # u is a state vector with:
            # u[1] = S relative saturation (-)
            # u[2] = δ2H [‰] TODO(bernhard): #modify this to work with concentrations instead of ‰
            # u[3] = δ18O [‰] TODO(bernhard): #modify this to work with concentrations instead of ‰
    ######

    NLAYER = length(aux_du_VRFLI)
    # SWAT
        # Actual water mass balance equation over computational cell between interfaces z_upper and z_lower:
        # ∂/∂t ∫θ dz = ∫∂/∂z[K_h(∂h/∂z + 1)] dz - ∫S dz  # (where S are source/sink terms)
        # ∂/∂t SWATI = [K_h(∂h/∂z + 1)]_(z_lower)^(z_upper) - ∫S dz
        #            = [K_h(∂h/∂z + 1)]_(z_lower)^(z_upper) + (INFLI-TRANI-DSFL-SLVP)
        #            = q(z_upper) - q(z_lower) + INFLI - TRANI - DSFL - SLVP

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
        # ∂/∂t Cᵢ = [- Cᵢ/SWATI * ∂/∂t SWATI] + 1/SWATI * [diff(z_upper) - diff(z_lower) - qCᵢ(z_upper) + qCᵢ(z_lower) + INFLI*Cᵢ_{INFLI} - TRANI*Cᵢ_{TRANI} - DSFL*Cᵢ_{DSFL} - SLVP*Cᵢ_{SLVP}]
        #
        # The composition of the fluxes are the comoposition of the originating storage for mixing processes.
        # For diffusion process they are computed based on the diffusivity.
        # For evaporation processes they are computed using the Craig-Gordon model:
        #   Cᵢ = [ (δᵢ - ε_eq)/α_eq - h*δₐ - ε_kin ] / [ 1 - h + ε_kin/1000] # [‰]
        # source: Benettin 2018 HESS eq. 1 and Gibson 2016 QSR eq. 3
        #         and Craig Gordon 1965 eq. 23 (and text above)
        #
        # Above equation is not used, instead one is derived from Gonfiantini 2018:
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


        ##################
        ##################
        ##################
        # Define conditions for isotope calculations:
        Tc = p_fu_TADTM  # °C, average daytime air temperature
        h = min(1.0, p_EA) # -, relative humidity of the atmosphere (vappress_atm/1 atm)
        γ = 1.0          # -, thermodynamic activity coefficient of evaporating water
        # X_INTS = 0.5  # -, turbulence incex of the atmosphere above the evaporating water
        # X_INTR = 0.5  # -, turbulence incex of the atmosphere above the evaporating water
        # X_SNOW = 1.0  # -, turbulence incex of the atmosphere above the evaporating water

        # 1c) Atmospheric vapor composition assumed to be in equilibrium with precipitation
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

        # Compute composition of evaporating soil
        ### a) Fix constants
        α¹⁸O_eq = LWFBrook90.ISO.α¹⁸O_eq(Tc)
        α²H_eq  = LWFBrook90.ISO.α²H_eq(Tc)
        # ε¹⁸O_dif = (α¹⁸O_dif-1)*1000  # diffusive (i.e. kinetic fractionation)
        # ε²H_dif  = (α²H_dif-1)*1000   # diffusive (i.e. kinetic fractionation)
        # ε¹⁸O_eq = (α¹⁸O_eq-1)*1000    # equilibrium fractionation
        # ε²H_eq  = (α²H_eq-1)*1000     # equilibrium fractionation

        ### b) Compute δ signature of evaporating flux
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

        ### c) Composition of all other fluxes
        Cᵢ¹⁸ = u_δ18O_SWATI
        Cᵢ²H = u_δ2H_SWATI

        diff_upp = 0 #TODO(bernhard): include diffusive flux of isotopes
        diff_low = 0 #TODO(bernhard): include diffusive flux of isotopes
        qCᵢ¹⁸O_upp  = [0; aux_du_VRFLI[1:(NLAYER-1)]] .* [0; u_δ18O_SWATI[1:(NLAYER-1)]]
        qCᵢ¹⁸O_low  = aux_du_VRFLI[1:(NLAYER)]        .* u_δ18O_SWATI[1:(NLAYER)]
        qCᵢ²H_upp  = [0; aux_du_VRFLI[1:(NLAYER-1)]] .* [0; u_δ2H_SWATI[1:(NLAYER-1)]]
        qCᵢ²H_low  = aux_du_VRFLI[1:(NLAYER)]        .* u_δ2H_SWATI[1:(NLAYER)]

        Cᵢ¹⁸O_INFLI = p_δ18O_PREC # TODO(bernhard): for debugging, remove this again and replace with δ18O_INFLI
        # Cᵢ¹⁸O_INFLI = δ18O_INFLI
        # Cᵢ¹⁸O_INFLI = ifelse(sum(aux_du_INFLI) == 0, 0, δ18O_INFLI) # in case there is no inflow δ18O_INFLI was set to NaN, set it to zero for below equation
        Cᵢ¹⁸O_TRANI = Cᵢ¹⁸ # no fractionation occurring, i.e. outflux composition equal to storage composition
        Cᵢ¹⁸O_DSFL  = Cᵢ¹⁸ # no fractionation occurring, i.e. outflux composition equal to storage composition
        Cᵢ¹⁸O_SLVP  = [δ¹⁸O_SLVP; fill(0, NLAYER-1)]
        Cᵢ²H_INFLI = p_δ2H_PREC # TODO(bernhard): for debugging, remove this again and replace with δ2H_INFLI
        # Cᵢ²H_INFLI = δ2H_INFLI
        # Cᵢ²H_INFLI = ifelse(sum(aux_du_INFLI) == 0, 0, δ2H_INFLI) # in case there is no inflow δ2H_INFLI was set to NaN, set it to zero for below equation
        Cᵢ²H_TRANI = Cᵢ²H # no fractionation occurring, i.e. outflux composition equal to storage composition
        Cᵢ²H_DSFL  = Cᵢ²H # no fractionation occurring, i.e. outflux composition equal to storage composition
        Cᵢ²H_SLVP  = [δ²H_SLVP;  fill(0, NLAYER-1)]

        ### c) Some other terms in the isotope balance equation
        dVdt = du_NTFLI # [mm/day]
        # θ = u_aux_θ
        # dθdt = aux_du_VRFLI[NLAYER] - du_GWFL - du_SEEP # dV/dt [-/day] # TODO(bernhard)

        ### e) Isotope balance equation:
        # du_δ18O_SWATI = -Cᵢ¹⁸./u_SWATI .* dVdt .+ 1 ./ u_SWATI .* (
        #                         diff_upp .- diff_low .- qCᵢ¹⁸O_upp .+ qCᵢ¹⁸O_low .+
        #                         aux_du_INFLI.*Cᵢ¹⁸O_INFLI .- aux_du_TRANI.*Cᵢ¹⁸O_TRANI .- aux_du_DSFLI.*Cᵢ¹⁸O_DSFL .- aux_du_SLVP.*Cᵢ¹⁸O_SLVP
        #                     )
        # du_δ2H_SWATI = -Cᵢ²H./u_SWATI .* dVdt .+ 1 ./ u_SWATI .* (
        #                         diff_upp .- diff_low .- qCᵢ²H_upp .+ qCᵢ²H_low .+
        #                         aux_du_INFLI.*Cᵢ²H_INFLI .- aux_du_TRANI.*Cᵢ²H_TRANI .- aux_du_DSFLI.*Cᵢ²H_DSFL .- aux_du_SLVP.*Cᵢ²H_SLVP
        #                     )
        # NOTE: below max(0.001,u_SWATI) makes the code more robust
        du_δ18O_SWATI = -Cᵢ¹⁸./max.(0.001,u_SWATI) .* dVdt .+ 1 ./ max.(0.001,u_SWATI) .* (
                                diff_upp .- diff_low .+ qCᵢ¹⁸O_upp .- qCᵢ¹⁸O_low .+ # TODO(bernhard): check sign for qC
                                aux_du_INFLI.*Cᵢ¹⁸O_INFLI .- aux_du_TRANI.*Cᵢ¹⁸O_TRANI .- aux_du_DSFLI.*Cᵢ¹⁸O_DSFL .- aux_du_SLVP.*Cᵢ¹⁸O_SLVP
                            )
        du_δ2H_SWATI = -Cᵢ²H./max.(0.001,u_SWATI) .* dVdt .+ 1 ./ max.(0.001,u_SWATI) .* (
                                diff_upp .- diff_low .+ qCᵢ²H_upp .- qCᵢ²H_low .+ # TODO(bernhard): check sign for qC
                                aux_du_INFLI.*Cᵢ²H_INFLI .- aux_du_TRANI.*Cᵢ²H_TRANI .- aux_du_DSFLI.*Cᵢ²H_DSFL .- aux_du_SLVP.*Cᵢ²H_SLVP
                            )
    # GWAT
        δ18O_empty = NaN
        δ2H_empty  = NaN
        @assert aux_du_VRFLI[NLAYER] >= 0 "aux_du_VRFLI[NLAYER] should not be negative"

        if ((u_GWAT == 0) & (aux_du_VRFLI[NLAYER] == 0)) # initially no groundwater and no new is added
            du_δ18O_GWAT = δ18O_empty
            du_δ2H_GWAT  = δ2H_empty
        elseif aux_du_VRFLI[NLAYER] > 0 # initially no groundwater but some is added
            du_δ18O_GWAT = 0
            du_δ2H_GWAT  = 0

            u_δ18O_GWAT = u_δ18O_SWATI[end]
            u_δ2H_GWAT  = u_δ2H_SWATI[end]  # TODO(benrhard): make sure this is actually modifying the state vector
        else
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
            du_δ18O_GWAT = in_GWAT/V_GWAT*(δ18O_in_GWAT - u_δ18O_GWAT) # - out_GWAT*(δ18O_out_GWAT - u_δ18O_GWAT) # <- last part == 0
            du_δ2H_GWAT  = in_GWAT/V_GWAT*(δ2H_in_GWAT  - u_δ2H_GWAT)  # - out_GWAT*(δ2H_out_GWAT  - u_δ2H_GWAT)  # <- last part == 0
        end


    return du_δ18O_GWAT, du_δ2H_GWAT, du_δ18O_SWATI, du_δ2H_SWATI
    # return 0, 0, convert.(Float64, 1:NLAYER), convert.(Float64, 1:NLAYER)
end