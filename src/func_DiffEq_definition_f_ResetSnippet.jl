
function apply_reset(integrator, DTINEW)
    # new!! update soil limited boundary flows and source/sink terms during iteration loop
    # print("\nUpdate boundary flows and source/sink terms during iteration loop in case they are soil limited.")
    # ***********************************************************

    # ***********************************************************
    # 0) Parse parameters
    ## 0) for debugging
    IDAY = floor(integrator.t) # TODO(bernhard) is just for debug, remove again after

    ## A) constant parameters:
    (p_DT, NLAYER, IMODEL, _, _, # 2nd _ = Reset
    p_SWATMX, p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_PSIG, p_KF,
    p_THSAT, p_θr, p_MvGα, p_MvGn, p_MvGl, p_Ksat,

    _, _, _, _, _, _,
    _, _, _, _,
    p_KSAT, _, _, _, _,
    _, _, p_THICK, p_STONEF,

    _) = integrator.p[1][1]

    (p_LAT, p_ESLOPE, p_L1, p_L2,
    p_SNODEN, p_MXRTLN, p_MXKPL, p_CS,
    p_Z0S, p_Z0G,
    p_ZMINH, p_CZS, p_CZR, p_HS, p_HR, p_LPC,
    p_RTRAD, p_FXYLEM,
    p_WNDRAT, p_FETCH, p_Z0W, p_ZW,
    p_RSTEMP,
    p_CVICE,
    p_LWIDTH, p_RHOTP, p_NN, p_KSNVP,
    p_ALBSN, p_ALB,
    p_RSSA, p_RSSB,
    p_CCFAC, p_MELFAC, p_LAIMLT, p_SAIMLT,

    p_inirdep, p_inirlen, p_rgroper, p_tini, p_frelden,

    p_WTOMJ, p_C1, p_C2, p_C3, p_CR,
    p_GLMIN, p_GLMAX, p_R5, p_CVPD, p_RM, p_TL, p_T1, p_T2, p_TH,
    p_PSICR, NOOUTF, p_PsiCrit,

    # for MSBPREINT:
    p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS,
    p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS,
    p_DURATN, p_MAXLQF, p_GRDMLT) = integrator.p[1][2]

    ## B) time dependent parameters
    p_DOY, p_MONTHN, p_SOLRAD, p_TMAX, p_TMIN, p_EA, p_UW, p_PRECIN, p_DTP, p_NPINT,
        _, p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_AGE = integrator.p[2] # _ = p_MESFL

    p_fT_DENSEF = max(0.050, p_DENSEF(integrator.t))

    # Compute time dependent root density parameters
    # TODO(bernhard): Do this outside of integration loop in define_DiffEq_parameters()
    p_fT_RELDEN = LWFBrook90Julia.WAT.LWFRootGrowth(p_frelden, p_tini, p_AGE(integrator.t), p_rgroper, p_inirdep, p_inirlen, NLAYER)

    # Compute rate of rain (mm/day)
    # TODO(bernhard): a) Do this outside of integration loop in define_DiffEq_parameters()
    #                 b) And simplify it directly to rate p_fT_PREC in both cases
    #                    i.e in case PREINT (PRECDAT) or in case
    if (isequal(p_DTP, 1))
        # p[2][9] # p_DTP
        # p[2][10] # p_NPINT

        # NOTE: Curently parameters overdetermine. We only need two out of the following three:
        # p_DTP = p_DT/p_NPINT
        # p_NPINT = p_DT/p_DTP

        p_fT_PREINT = p_PRECIN(integrator.t) / p_DTP # (mm/day)
        # TODO(bernhard): rename p_NPINT to p_NPINT
    else
        error("Case where input file PRECDAT is used is not implemented.
                Reading PRECDAT should result in PREINT (precipitation amount per interval)")
    end

    ## C) state dependent parameters:
    # Calculate parameters:
    #  - solar parameters depending on DOY
    #  - canopy parameters depending on DOY and u_SNOW
    #  - roughness parameters depending on u_SNOW
    #  - plant resistance components depending on u_SNOW
    #  - weather data depending on DOY and u_SNOW
    #  - fraction of precipitation as snowfall depending on DOY
    #  - snowpack temperature, potential snow evaporation and soil evaporation resistance depending on u_SNOW
    u_INTS     = integrator.u[2]
    u_INTR     = integrator.u[3]
    u_SNOW     = integrator.u[4]
    u_CC       = integrator.u[5]
    u_SNOWLQ   = integrator.u[6]
    u_SWATI    = integrator.u[7:(7+NLAYER-1)]

    # Derive (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) from u_SWATI
    (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        LWFBrook90Julia.KPT.derive_auxiliary_SOILVAR(u_SWATI, p_SWATMX, p_THSAT,
            p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_KF,
            p_θr, p_MvGα, p_MvGn, p_MvGl, p_Ksat,
            p_PSIG, NLAYER, IMODEL)

    p_fT_DAYLEN, p_fT_I0HDAY, p_fT_SLFDAY, p_fu_HEIGHT, p_fu_LAI, p_fu_SAI, _, # _ = p_fu_RTLEN
        _, p_fu_Z0GS, p_fu_Z0C, p_fu_DISPC, p_fu_Z0, p_fu_DISP, p_fu_ZA, # _ = p_fu_RPLANT,
        p_fu_RXYLEM, p_fu_RROOTI, p_fu_ALPHA,
        p_fu_SHEAT,p_fu_SOLRADC, p_fu_TA, p_fu_TADTM, p_fu_TANTM, p_fu_UADTM, p_fu_UANTM,
        p_fT_SNOFRC,_,p_fu_PSNVP, p_fu_ALBEDO,_, p_fu_SNOEN = # _ = p_fu_TSNOW, p_fu_RSS
        MSBSETVARS(IDAY, #TODO(bernhard) just for debug... remove again!
                    IMODEL,
                    # for SUNDS:
                    p_LAT, p_ESLOPE, p_DOY(integrator.t), p_L1, p_L2,
                    # for LWFBrook90_CANOPY:
                    p_HEIGHT(integrator.t), p_LAI(integrator.t), p_SAI(integrator.t), u_SNOW, p_SNODEN, p_MXRTLN, p_MXKPL, p_fT_DENSEF,
                    #
                    p_Z0S, p_Z0G,
                    # for ROUGH:
                    p_ZMINH, p_CZS, p_CZR, p_HS, p_HR, p_LPC, p_CS,
                    # for PLNTRES:
                    NLAYER, p_THICK, p_STONEF, p_fT_RELDEN, p_RTRAD, p_FXYLEM,
                    # for WEATHER:
                    p_TMAX(integrator.t), p_TMIN(integrator.t), p_EA(integrator.t), p_UW(integrator.t), p_WNDRAT, p_FETCH, p_Z0W, p_ZW, p_SOLRAD(integrator.t),
                    # for SNOFRAC:
                    p_RSTEMP,
                    #
                    u_CC, p_CVICE,
                    # for SNOVAP:
                    p_LWIDTH, p_RHOTP, p_NN, p_KSNVP,
                    #
                    p_ALBSN, p_ALB,
                    # for FRSS:
                    p_RSSA, p_RSSB, p_PSIF, u_aux_PSIM, p_PsiCrit, #p_PSIF[1], u_aux_PSIM[1]
                    # for SNOENRGY:
                    p_CCFAC, p_MELFAC, p_LAIMLT, p_SAIMLT)
    # End parse parameters
    # ***********************************************************
    # Now start Reset code:


    # 1) Update p_fu_RSS # Shuttleworth-Wallace soil surface resistance, s/m
    #    using recomputed u_aux_PSIM
    # Part of MSBSETVARS:
    if u_SNOW > 0.0
        p_fu_RSS = 0.0
    else
        p_fu_RSS = LWFBrook90Julia.PET.FRSS(p_RSSA, p_RSSB, p_PSIF[1], u_aux_PSIM[1], p_PsiCrit[1])
    end
    #print("\nIDAY:$(@sprintf("% 3d", IDAY)), p_RSSA: $(@sprintf("% 8.4f",p_RSSA)), p_RSSB: $(@sprintf("% 8.4f",p_RSSB)), p_PSIF[1]: $(@sprintf("% 8.4f",p_PSIF[1])), u_aux_PSIM[1]: $(@sprintf("% 8.4f",u_aux_PSIM[1])), p_PsiCrit[1]: $(@sprintf("% 8.4f",p_PsiCrit[1]))")


    # 2) Update p_fu_PTRAN, p_fu_GEVP, p_fu_PINT, p_fu_GIVP, p_fu_PSLVP, aux_du_TRANI
    #    using recomputed:     p_fu_RSS, u_aux_PSITI, p_fu_KK (total soil water potential, kPa; hydraulic conductivity, mm/d)
    #    Note that it is
    #    using not recomputed: p_fu_SOLRADC, p_fu_TADTM, p_fu_UADTM, p_fu_TANTM, p_fu_UANTM
    #                          p_fu_ALBEDO, p_fu_SHEAT, p_fu_LAI, p_fu_SAI, p_fu_ZA, p_fu_HEIGHT
    #                          p_fu_Z0, p_fu_Z0C, p_fu_Z0GS, p_fu_DISP, p_fu_DISPC,
    #                          p_fu_TA, p_fu_ALPHA, p_fu_RROOTI, p_fu_RXYLEM
    #* * * * *  B E G I N   D A Y - N I G H T   E T   L O O P  * * * * * * * * *
    # Compute day and night rates
    (p_fu_PTR, p_fu_GER, p_fu_PIR, p_fu_GIR, p_fu_ATRI, p_fu_PGER) =
        MSBDAYNIGHT(IDAY,
                    IMODEL,
                    p_fT_SLFDAY, p_fu_SOLRADC, p_WTOMJ, p_fT_DAYLEN, p_fu_TADTM, p_fu_UADTM, p_fu_TANTM, p_fu_UANTM,
                    p_fT_I0HDAY,
                    # for AVAILEN:
                    p_fu_ALBEDO, p_C1, p_C2, p_C3, p_EA(integrator.t), p_fu_SHEAT, p_CR, p_fu_LAI, p_fu_SAI,
                    # for SWGRA:
                    p_fu_ZA, p_fu_HEIGHT, p_fu_Z0, p_fu_DISP, p_fu_Z0C, p_fu_DISPC, p_fu_Z0GS, p_LWIDTH, p_RHOTP, p_NN,
                    # for SRSC:
                    p_fu_TA, p_GLMIN, p_GLMAX, p_R5, p_CVPD, p_RM, p_TL, p_T1, p_T2, p_TH,
                    # for SWPE:
                    p_fu_RSS,
                    # for TBYLAYER:
                    p_fu_ALPHA, p_fu_KK, p_fu_RROOTI, p_fu_RXYLEM, u_aux_PSITI, NLAYER, p_PSICR, NOOUTF)

    # Combine day and night rates to average daily rate
    (p_fu_PTRAN, p_fu_GEVP, p_fu_PINT, p_fu_GIVP, p_fu_PSLVP, aux_du_TRANI) =
        MSBDAYNIGHT_postprocess(IMODEL, NLAYER, p_fu_PTR, p_fu_GER, p_fu_PIR, p_fu_GIR, p_fu_ATRI, p_fT_DAYLEN, p_fu_PGER, p_DT)

    # 3) Update: p_fu_PTRAN, aux_du_TRANI, aux_du_SLVP, aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP
    #    using recomputed:     u_aux_PSIM, p_fu_PSLVP, p_fu_PINT, p_fu_PTRAN, aux_du_TRANI, p_fu_GIVP, p_fu_GEVP
    #    Note that it is
    #    using not recomputed: p_fu_TA, p_fu_LAI, p_fu_SAI, p_fu_PSNVP, p_fu_SNOEN,

    (# compute some fluxes as intermediate results:
        #p_fT_SFAL, p_fT_RFAL, p_fu_RNET, p_fu_PTRAN,
        _, _, _, p_fu_PTRAN, # NOTE: Don't modify these variable as this is the reset run!
        # compute changes in soil water storage:
        aux_du_TRANI, aux_du_SLVP,
        # compute change in interception storage:
        aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP,
        # compute change in snow storage:
        #aux_du_RSNO, aux_du_SNVP, aux_du_SMLT,
        _,_,_,     # NOTE: Don't modify these variable as this is the reset run!
        # compute updated states:
        #u_SNOW, u_CC, u_SNOWLQ) =
        _, _, _) = # NOTE: Don't modify these variable as this is the reset run!
        MSBPREINT(p_fT_PREINT, p_DTP, p_fT_SNOFRC, p_NPINT, p_fu_PINT, p_fu_TA,
                    # for INTER (snow)
                    u_INTS, p_fu_LAI, p_fu_SAI, p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS,
                    # for INTER (rain)
                    u_INTR, p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS,
                    # for INTER24 (snow + rain)
                    p_DURATN, p_MONTHN(integrator.t),
                    #
                    u_SNOW, p_fu_PTRAN, NLAYER, aux_du_TRANI, p_fu_GIVP, p_fu_GEVP,
                    # for SNOWPACK
                    u_CC, u_SNOWLQ, p_fu_PSNVP, p_fu_SNOEN, p_MAXLQF, p_GRDMLT)


        ########################
        # Apply corrected quantities

        # NOTE(bernhard): Note that in the original LWFBrook90R (2021-02-20) aux_du_TRANI and
        #                 aux_du_SLVP were only applied to the next time step after the current one
        #                 In LWFBrook90Julia, this takes immediate effect for the current time step.
        integrator.p[3][3] = aux_du_TRANI
        integrator.p[3][4] = aux_du_SLVP

        # NOTE(bernhard): Note that in the original LWFBrook90R (2021-02-20) PACCUM.h
        #                 is called before the Reset and therefore Reset quantities
        #                 are not updated in PACCUM for the current DTINEW
        #                 (e.g. TRANI, PSLVP, SLVP, SRFL, SLFL, GWFL, SEEP)
        #                 but they are updated for the next DTINEW
        #                 ==> In LWFBrook90Julia these quantities are updated already for the current DTINEW
        #
        #                 Update in PSUM.h and DACCUM.h happens after exiting iteration loop.
        #                 Upon exit of iteration loop no Reset recalculation is performed

        # IF Reset == 0 set below quantities need to be set entirely in daily callback only.
        # If reset == 1 set below quantities to zero in daily callback and add them up *DTINEW here.
        # Quantities entirely recomputed in Reset
        integrator.u[7+NLAYER+ 3] += DTINEW*(aux_du_RINT)
        integrator.u[7+NLAYER+ 4] += DTINEW*(aux_du_SINT)
        integrator.u[7+NLAYER+ 9] += DTINEW*(sum(aux_du_TRANI))
        integrator.u[7+NLAYER+10] += DTINEW*(aux_du_IRVP)
        integrator.u[7+NLAYER+11] += DTINEW*(aux_du_ISVP)
        integrator.u[7+NLAYER+12] += DTINEW*(aux_du_SLVP)
        integrator.u[7+NLAYER+14] += DTINEW*(p_fu_PINT)
        integrator.u[7+NLAYER+15] += DTINEW*(p_fu_PTRAN)
        integrator.u[7+NLAYER+16] += DTINEW*(p_fu_PSLVP)
        # Quantities only partially recomputed in Reset
        integrator.u[7+NLAYER+ 6] += DTINEW*(-aux_du_RINT)
                                    # cum_d_rnet: p_fT_RFAL and aux_du_RSNO is not recomputed in Reset
        integrator.u[7+NLAYER+ 8] += DTINEW*(aux_du_IRVP + aux_du_ISVP + aux_du_SLVP + sum(aux_du_TRANI))
                                    # cum_d_evap: aux_du_SNVP is not recomputed in Reset
end