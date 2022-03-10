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
    p_fu_SHEAT = 0
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
        p_fu_RSS = 0
    else
        p_fu_TSNOW = 0
        p_fu_PSNVP = 0
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
                   #
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
                   u_CC, u_SNOWLQ, p_fu_PSNVP, p_fu_SNOEN, p_MAXLQF, p_GRDMLT)

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
                   p_DTP, p_fu_TA, p_MAXLQF, p_GRDMLT)

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
            aux_du_RSNO, aux_du_SNVP, aux_du_SMLT,
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
    aux_du_VRFLI = fill(NaN, NLAYER)
    aux_du_DSFLI = fill(NaN, NLAYER)

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
                                    aux_du_VRFLI_1st_approx)
    #@debug "b) aux_du_VRFLI[1]: $(aux_du_VRFLI[1]), sum(aux_du_VRFLI): $(sum(aux_du_VRFLI))"

    DPSIDW = LWFBrook90.KPT.FDPSIDWF(u_aux_WETNES, p_soil)

    # limit step size
    #       TODO(bernhard): This could alternatively be achieved with a stepsize limiter callback
    #                       see: https://diffeq.sciml.ai/stable/features/callback_library/#Stepsize-Limiters
    #       TODO(bernhard): Or it could alternatively be set manually using set_proposed_dt!
    #                       see: https://diffeq.sciml.ai/stable/basics/integrator/#SciMLBase.set_proposed_dt!
    #   ITER computes DTI so that the potential difference (due to aux_du_VRFLI)
    #   between adjacent layers does not change sign during the iteration time step
    DTINEW=LWFBrook90.WAT.ITER(NLAYER, FLAG_MualVanGen, DTI, LWFBrook90.CONSTANTS.p_DTIMIN, DPSIDW,
                                    du_NTFLI, u_aux_PSITI, u_aux_θ, p_DSWMAX, p_DPSIMAX,
                    p_soil)
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