"""
    define_LWFB90_cb()

Generate callback function cb needed for ODE() problem in DiffEq.jl package.

LWFBrook90 updates states INTS, INTR, SNOW, CC, SNOWLQ not continuously but only
once per day. This operator splitting (daily vs continuous update of ODEs) is
implemented by using this callback function which is called once per day.
"""
# B) Define callback
function define_LWFB90_cb()
    # cb_LaiRichardsCorrectorStep = FunctionCallingCallback(
    #     (u,t,integrator) -> nothing;    # TODO(bernhard): this can be used to implement the corrector step of Lai-2015 (or Li-2021-J_Hydrol.pdf) after DiffEq.jl is used for the predictor step
    #     func_everystep=true, func_start=false)

    cb_check_balance_errors = PeriodicCallback(LWFBrook90R_check_balance_errors!,  1.0;
                                initial_affect = false, # we do not need to check balance errors at the initial conditions...
                                save_positions=(false,false));

    cb_INTS_INTR_SNOW_amounts = PeriodicCallback(LWFBrook90R_updateAmounts_INTS_INTR_SNOW_CC_SNOWLQ!,  1.0;
                                initial_affect = true,
                                save_positions=(false,false));

    cb_INTS_INTR_SNOW_deltas = PeriodicCallback(LWFBrook90R_updateIsotopes_INTS_INTR_SNOW!,  1.0;
                               initial_affect = true,
                               save_positions=(false,false));

    cb_SWAT_GWAT_deltas = FunctionCallingCallback(LWFBrook90R_updateIsotopes_GWAT_SWAT!;
                                func_everystep = true,
                                func_start = false);

    cb_set = CallbackSet(
        # cb_LaiRichardsCorrectorStep,
        cb_check_balance_errors,
        cb_INTS_INTR_SNOW_amounts,
        cb_INTS_INTR_SNOW_deltas,
        cb_SWAT_GWAT_deltas
        )
    return cb_set
end

# A) Define updating function
function LWFBrook90R_updateAmounts_INTS_INTR_SNOW_CC_SNOWLQ!(integrator)
    # NOTE: we can make use of those:
    # integrator.t
    # integrator.p
    # integrator.u

    ############
    ### Compute parameters
    ## A) constant parameters:
    p_soil = integrator.p[1][1]
    (NLAYER, FLAG_MualVanGen, compute_intermediate_quantities, Reset,
    p_DTP, p_NPINT,

    _, _, _, _, _, _,
    _, _, _, _,
    _, _, _, _,
    _, _,

    _) = integrator.p[1][2]

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

    p_WTOMJ, p_C1, p_C2, p_C3, p_CR,
    p_GLMIN, p_GLMAX, p_R5, p_CVPD, p_RM, p_TL, p_T1, p_T2, p_TH,
    p_PSICR, NOOUTF,

    # for MSBPREINT:
    p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS,
    p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS,
    p_DURATN, p_MAXLQF, p_GRDMLT) = integrator.p[1][3]

    ## B) time dependent parameters
    p_DOY, p_MONTHN, p_SOLRAD, p_TMAX, p_TMIN, p_EA, p_UW, p_PREC,
        p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_AGE, p_RELDEN,
        p_δ18O_PREC, p_δ2H_PREC, ref_date = integrator.p[2]

    p_fT_DENSEF = max(0.050, p_DENSEF(integrator.t))

    ## C) state dependent parameters:
    # Calculate parameters:
    #  - solar parameters depending on DOY
    #  - canopy parameters depending on DOY and u_SNOW
    #  - roughness parameters depending on u_SNOW
    #  - plant resistance components depending on u_SNOW
    #  - weather data depending on DOY and u_SNOW
    #  - fraction of precipitation as snowfall depending on DOY
    #  - snowpack temperature, potential snow evaporation and soil evaporation resistance depending on u_SNOW

    idx_u_vector_amounts       = integrator.p[1][4][4]
    idx_u_vector_accumulators  = integrator.p[1][4][5]

    u_INTS     = integrator.u[2]
    u_INTR     = integrator.u[3]
    u_SNOW     = integrator.u[4]
    u_CC       = integrator.u[5]
    u_SNOWLQ   = integrator.u[6]
    u_SWATI    = integrator.u[idx_u_vector_amounts]

    LWFBrook90.KPT.SWCHEK!(u_SWATI, p_soil.p_SWATMAX, integrator.t)

    # Derive (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) from u_SWATI
    (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        LWFBrook90.KPT.derive_auxiliary_SOILVAR(u_SWATI, p_soil)

    IDAY = floor(integrator.t) # TODO(bernhard) is just for debug, remove again after
    p_fT_DAYLEN, p_fT_I0HDAY, p_fT_SLFDAY, p_fu_HEIGHT, p_fu_LAI, p_fu_SAI, p_fT_RTLEN,
        p_fT_RPLANT,p_fu_Z0GS, p_fu_Z0C, p_fu_DISPC, p_fu_Z0, p_fu_DISP, p_fu_ZA,
      p_fT_RXYLEM, p_fT_RROOTI, p_fT_ALPHA,
      p_fu_SHEAT,p_fu_SOLRADC, p_fu_TA, p_fu_TADTM, p_fu_TANTM, p_fu_UADTM, p_fu_UANTM,
      p_fT_SNOFRC,p_fu_TSNOW,p_fu_PSNVP, p_fu_ALBEDO,p_fu_RSS, p_fu_SNOEN =
      MSBSETVARS(IDAY, #TODO(bernhard) just for debug... remove again!
                 FLAG_MualVanGen, NLAYER, p_soil,
                 # for SUNDS:
                 p_LAT, p_ESLOPE, p_DOY(integrator.t), p_L1, p_L2,
                 # for LWFBrook90_CANOPY:
                 p_HEIGHT(integrator.t), p_LAI(integrator.t), p_SAI(integrator.t), u_SNOW, p_SNODEN, p_MXRTLN, p_MXKPL, p_fT_DENSEF,
                 #
                 p_Z0S, p_Z0G,
                 # for ROUGH:
                 p_ZMINH, p_CZS, p_CZR, p_HS, p_HR, p_LPC, p_CS,
                 # for PLNTRES:
                 p_RELDEN.(integrator.t, 1:NLAYER), p_RTRAD, p_FXYLEM,
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
                 p_RSSA, p_RSSB, u_aux_PSIM, #u_aux_PSIM[1]
                 # for SNOENRGY:
                 p_CCFAC, p_MELFAC, p_LAIMLT, p_SAIMLT)


    # Calculate average daily rate of potential and actual interception,
    # evaporation, and transpiration by considering weighted average of rate
    # during day and rate during night:
    #* * * * *  B E G I N   D A Y - N I G H T   E T   L O O P  * * * * * * * * *
    # Compute day and night rates
    (p_fu_PTR, p_fu_GER, p_fu_PIR, p_fu_GIR, p_fu_ATRI, p_fu_PGER) =
        MSBDAYNIGHT(IDAY,
                    FLAG_MualVanGen,
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
                    p_fT_ALPHA, p_fu_KK, p_fT_RROOTI, p_fT_RXYLEM, u_aux_PSITI, NLAYER, p_PSICR, NOOUTF)
                    # 0.000012 seconds (28 allocations: 1.938 KiB)
    # Combine day and night rates to average daily rate
    (p_fu_PTRAN, p_fu_GEVP, p_fu_PINT, p_fu_GIVP, p_fu_PSLVP, aux_du_TRANI) = # TODO(bernhard): p_fu_PSLVP is unused
        MSBDAYNIGHT_postprocess(FLAG_MualVanGen, NLAYER, p_fu_PTR, p_fu_GER, p_fu_PIR, p_fu_GIR, p_fu_ATRI, p_fT_DAYLEN, p_fu_PGER)
    #* * * * * * * *  E N D   D A Y - N I G H T   L O O P  * * * * * * * * * *

    ####################################################################
    # 1) Update snow accumulation/melt: u_SNOW, u_CC, u_SNOWLQ
    #    and compute fluxes to/from interception storage
    u_SNOW_old = u_SNOW
    (# compute some fluxes as intermediate results:
    p_fT_SFAL, p_fT_RFAL, p_fu_RNET, p_fu_PTRAN,
    # compute changes in soil water storage:
    aux_du_TRANI, aux_du_SLVP,
    # compute change in interception storage:
    aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP,
    # compute change in snow storage:
    aux_du_RSNO, aux_du_SNVP, aux_du_SMLT, p_fu_STHR,
    # compute updated states:
    u_SNOW_MSBupdate, u_CC, u_SNOWLQ) =
        MSBPREINT(p_PREC(integrator.t), p_DTP, p_fT_SNOFRC, p_NPINT, p_fu_PINT, p_fu_TA,
               # for INTER (snow)
               u_INTS, p_fu_LAI, p_fu_SAI, p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS,
               # for INTER (rain)
               u_INTR, p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS,
               # for INTER24 (snow + rain)
               p_DURATN, p_MONTHN(integrator.t),
               #
               u_SNOW, p_fu_PTRAN, NLAYER, aux_du_TRANI, p_fu_GIVP, p_fu_GEVP,
               # for SNOWPACK
               u_CC, u_SNOWLQ, p_fu_PSNVP, p_fu_SNOEN, p_MAXLQF, p_GRDMLT,
               # Constants
               LWFBrook90.CONSTANTS.p_CVICE, LWFBrook90.CONSTANTS.p_LF, LWFBrook90.CONSTANTS.p_CVLQ)
               # 0.000016 seconds (28 allocations: 3.609 KiB)

    ####################################################################
    # 2) Update states of interception storage over entire precipitation
    #    interval: u_INTS, u_INTR
    u_INTS = u_INTS + (aux_du_SINT - aux_du_ISVP) * p_DTP
    u_INTR = u_INTR + (aux_du_RINT - aux_du_IRVP) * p_DTP
    u_SNOW = u_SNOW_MSBupdate

    ####################################################################
    # 3) Update soil water using substeps smaller than precipitation
    #    interval. After each substep update: u_GWAT, u_SWATI
    # ... NOTE: this is not done in the callback function but in the main ODE solver

    ####################################################################
    # Return results from callback
    # update INTS
    integrator.u[2] = u_INTS

    # update INTR
    integrator.u[3] = u_INTR # INTR

    # update SNOW, CC, SNOWLQ
    integrator.u[4] = u_SNOW # SNOW
    integrator.u[5] = u_CC # CC
    integrator.u[6] = u_SNOWLQ # SNOWLQ

    # save intermediate results for use in ODE (function f()) or other callbacks
    # integrator.p[3][1] .= [δ18O_SLFL, δ2H_SLFL,
    #                        p_fu_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP,
    #                        p_fu_STHR, aux_du_RSNO, aux_du_SNVP,
    #                        aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, u_SNOW_old]
    integrator.p[3][1] .= [NaN, NaN,
                            p_fu_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP,
                            p_fu_STHR, aux_du_RSNO, aux_du_SNVP,
                            aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, u_SNOW_old]
    integrator.p[3][2] .= aux_du_TRANI


    ##########################################
    # Accumulate flows to compute daily sums
    # Note that below state variables serve only as accumulator but do not affect
    # the evolution of the system.
    if compute_intermediate_quantities

        # 1) Either set daily sum if rate is constant throughout precipitation interval: p_DTP*(...)
        # 2) or then set daily sum to zero and use ODE to accumulate flow.
        integrator.u[idx_u_vector_accumulators[ 1]] = p_DTP * (p_fT_RFAL + p_fT_SFAL)                 # RFALD + SFALD        # cum_d_prec
        integrator.u[idx_u_vector_accumulators[ 2]] = p_DTP * (p_fT_RFAL)                                                    # cum_d_rfal
        integrator.u[idx_u_vector_accumulators[ 3]] = p_DTP * (p_fT_SFAL)                                                    # cum_d_sfal
        integrator.u[idx_u_vector_accumulators[ 4]] = p_DTP * (aux_du_RINT)                                                  # cum_d_rint
        integrator.u[idx_u_vector_accumulators[ 5]] = p_DTP * (aux_du_SINT)                                                  # cum_d_sint
        integrator.u[idx_u_vector_accumulators[ 6]] = p_DTP * (aux_du_RSNO)                                                  # cum_d_rsno
        integrator.u[idx_u_vector_accumulators[ 7]] = p_DTP * (p_fT_RFAL - aux_du_RINT - aux_du_RSNO) # cum_d_RTHR - RSNOD   # cum_d_rnet
        integrator.u[idx_u_vector_accumulators[ 8]] = p_DTP * (aux_du_SMLT)                                                  # cum_d_smlt

        integrator.u[idx_u_vector_accumulators[ 9]] = p_DTP * (aux_du_IRVP + aux_du_ISVP + aux_du_SNVP + aux_du_SLVP + sum(aux_du_TRANI))  # cum_d_evap
        integrator.u[idx_u_vector_accumulators[10]] = p_DTP * (sum(aux_du_TRANI))                                                          # cum_d_tran
        integrator.u[idx_u_vector_accumulators[11]] = p_DTP * (aux_du_IRVP)                                                                # cum_d_irvp
        integrator.u[idx_u_vector_accumulators[12]] = p_DTP * (aux_du_ISVP)                                                                # cum_d_isvp
        integrator.u[idx_u_vector_accumulators[13]] = p_DTP * (aux_du_SLVP)                                                                # cum_d_slvp
        integrator.u[idx_u_vector_accumulators[14]] = p_DTP * (aux_du_SNVP)                                                                # cum_d_snvp
        integrator.u[idx_u_vector_accumulators[15]] = p_DTP * (p_fu_PINT)                                                                  # cum_d_pint
        integrator.u[idx_u_vector_accumulators[16]] = p_DTP * (p_fu_PTRAN)                                                                 # cum_d_ptran
        integrator.u[idx_u_vector_accumulators[17]] = p_DTP * (p_fu_PSLVP)                                                                 # cum_d_pslvp

        integrator.u[idx_u_vector_accumulators[18]] = 0 # flow,  is computed in ODE
        integrator.u[idx_u_vector_accumulators[19]] = 0 # seep,  is computed in ODE
        integrator.u[idx_u_vector_accumulators[20]] = 0 # srfl,  is computed in ODE
        integrator.u[idx_u_vector_accumulators[21]] = 0 # slfl,  is computed in ODE
        integrator.u[idx_u_vector_accumulators[22]] = 0 # byfl,  is computed in ODE
        integrator.u[idx_u_vector_accumulators[23]] = 0 # dsfl,  is computed in ODE
        integrator.u[idx_u_vector_accumulators[24]] = 0 # gwfl,  is computed in ODE
        integrator.u[idx_u_vector_accumulators[25]] = 0 # vrfln, is computed in ODE

        # integrator.u[idx_u_vector_accumulators[26]] = p_DTP*(p_fT_RFAL - aux_du_RINT) # cum_d_rthr
        # integrator.u[idx_u_vector_accumulators[27]] = p_DTP*(p_fT_SFAL - aux_du_SINT) # cum_d_sthr

        # below are computed in separate callback:
        # integrator.u[idx_u_vector_accumulators[28]] = new_SWAT
        # integrator.u[idx_u_vector_accumulators[29]] = new_totalWATER
        # integrator.u[idx_u_vector_accumulators[30]] = BALERD_SWAT
        # integrator.u[idx_u_vector_accumulators[31]] = BALERD_total

        # TODO(bernhard): use SavingCallback() for all quantities that have u=... and du=0
        #                 only keep du=... for quantities for which we compute cumulative sums
    end

    return nothing
    ##########################################
end

function LWFBrook90R_updateIsotopes_INTS_INTR_SNOW!(integrator)
    simulate_isotopes = integrator.p[1][4][3]

    if simulate_isotopes
        ## C) state dependent parameters:
        # Calculate parameters:
        #  - solar parameters depending on DOY
        #  - canopy parameters depending on DOY and u_SNOW
        #  - roughness parameters depending on u_SNOW
        #  - plant resistance components depending on u_SNOW
        #  - weather data depending on DOY and u_SNOW
        #  - fraction of precipitation as snowfall depending on DOY
        #  - snowpack temperature, potential snow evaporation and soil evaporation resistance depending on u_SNOW

        # idx_u_vector_amounts       = integrator.p[1][4][4]
        # idx_u_vector_accumulators  = integrator.p[1][4][5]

        u_INTS     = integrator.u[2]
        u_INTR     = integrator.u[3]
        u_SNOW     = integrator.u[4]
        # u_CC       = integrator.u[5]
        # u_SNOWLQ   = integrator.u[6]
        # u_SWATI    = integrator.u[idx_u_vector_amounts]

        ####################################################################
        # 2) Update states of isotopic compositions of interception storages and snowpack:
        #       u_δ18O_INTS, u_δ18O_INTR, u_δ18O_SNOW
        #       u_δ2H_INTS, u_δ2H_INTR, u_δ2H_SNOW
        ############
        ### Compute parameters
        ## A) constant parameters:
        (_, _, _, _,
        p_DTP, _,

        _, _, _, _, _, _,
        _, _, _, _,
        _, _, _, _,
        _, _,

        _) = integrator.p[1][2]

        ## B) time dependent parameters
        p_DOY, p_MONTHN, p_SOLRAD, p_TMAX, p_TMIN, p_EA, p_UW, p_PREC,
            p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_AGE, p_RELDEN,
            p_δ18O_PREC, p_δ2H_PREC, ref_date = integrator.p[2]

        ## C) state dependent parameters or intermediate results:
        # These were computed in the callback and are kept constant in between two
        # callbacks.
        (_, _, p_fu_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP,
            p_fu_STHR, aux_du_RSNO, aux_du_SNVP,
            aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, u_SNOW_old) = integrator.p[3][1]
        aux_du_TRANI = integrator.p[3][2]

        idx_u_scalar_isotopes_d18O = integrator.p[1][4][8]
        # idx_u_vector_isotopes_d18O = integrator.p[1][4][9]
        idx_u_scalar_isotopes_d2H  = integrator.p[1][4][10]
        # idx_u_vector_isotopes_d2H  = integrator.p[1][4][11]

        # u_δ18O_GWAT = integrator.u[idx_u_scalar_isotopes_d18O[1]]
        # u_δ2H_GWAT  = integrator.u[idx_u_scalar_isotopes_d2H[1] ]
        u_δ18O_INTS = integrator.u[idx_u_scalar_isotopes_d18O[2]]
        u_δ2H_INTS  = integrator.u[idx_u_scalar_isotopes_d2H[2] ]
        u_δ18O_INTR = integrator.u[idx_u_scalar_isotopes_d18O[3]]
        u_δ2H_INTR  = integrator.u[idx_u_scalar_isotopes_d2H[3] ]
        u_δ18O_SNOW = integrator.u[idx_u_scalar_isotopes_d18O[4]]
        u_δ2H_SNOW  = integrator.u[idx_u_scalar_isotopes_d2H[4] ]
        # u_δ18O_SWATI = integrator.u[idx_u_vector_isotopes_d18O]
        # u_δ2H_SWATI  = integrator.u[idx_u_vector_isotopes_d2H]

        # Variant 1) works but really much slower
        # δ18O_SLFL, δ2H_SLFL,
        # _,            u_δ18O_INTS, u_δ2H_INTS,
        # _,            u_δ18O_INTR, u_δ2H_INTR,
        # u_SNOW_iso_update, u_δ18O_SNOW, u_δ2H_SNOW =
        #     compute_isotope_U_of_INTS_INTR_SNOW_and_SLFL(
        # (δ18O_SLFL, δ2H_SLFL,
        # u_δ18O_INTS, u_δ2H_INTS,
        # u_δ18O_INTR, u_δ2H_INTR,
        # u_δ18O_SNOW, u_δ2H_SNOW) =
        #     compute_isotope_U_of_INTS_INTR_SNOW_and_SLFL(
        #         u_INTS, u_δ18O_INTS, u_δ2H_INTS, u_INTR, u_δ18O_INTR, u_δ2H_INTR, u_SNOW, u_δ18O_SNOW, u_δ2H_SNOW,
        #         p_δ2H_PREC(integrator.t), p_δ18O_PREC(integrator.t), p_fu_TADTM, p_EA(integrator.t),
        #         # for INTS (in: SINT; out: ISVP):
        #         aux_du_SINT, aux_du_ISVP, p_DTP,
        #         # for INTR (in: RINT; out: IRVP):
        #         aux_du_RINT, aux_du_IRVP,
        #         # for SNOW (in: STHR, RSNO (both δ_PREC); out: SMLT, SNVP (δ_SNOW and fractionated)):
        #         NaN, p_fu_STHR, aux_du_RSNO, aux_du_SMLT, aux_du_SNVP,
        #         # to compute isotopic signature of soil infiltration: SLFL
        #         p_fu_RNET)
        # # END variant 1

        # Variant 2) resulting in no effect (but is fast, so not generating any allocation problems)
        δ2H_SLFL   = p_δ2H_PREC(integrator.t)
        u_δ2H_INTS = p_δ2H_PREC(integrator.t)
        u_δ2H_INTR = p_δ2H_PREC(integrator.t)
        u_δ2H_SNOW = p_δ2H_PREC(integrator.t)
        δ18O_SLFL   = p_δ18O_PREC(integrator.t)
        u_δ18O_INTS = p_δ18O_PREC(integrator.t)
        u_δ18O_INTR = p_δ18O_PREC(integrator.t)
        u_δ18O_SNOW = p_δ18O_PREC(integrator.t)
        # # END variant 2

        # # Variant 3) trying to unwrap variant 1) by reducing allocations. Currently not working and also not faster than 1)
        # R_std_2H  = LWFBrook90.ISO.R_VSMOW²H
        # R_std_18O = LWFBrook90.ISO.R_VSMOW¹⁸O
        # # For SNOW ( do this before INTS as it uses u_δ2H_INTS as input)
        # inflow = [p_fu_STHR; aux_du_RSNO]
        # outflow = [aux_du_SMLT]
        # u⁺ = u_SNOW
        # u₀ = u_INTR - p_DTP * (sum(inflow) - sum(outflow)) # [mm]
        # E  = aux_du_SNVP
        # ## δ2H
        # δ₀ = u_δ2H_SNOW
        # x₀  = (δ₀/1000 + 1)*R_std_2H / (1 + (δ₀/1000 + 1)*R_std_2H)
        # δin = [u_δ2H_INTS;  p_δ2H_PREC(integrator.t)]
        # xin = (δin ./ 1000. .+ fill(1.,(2)))
        # xin = xin  .* R_std_2H ./ ((δin ./ 1000. .+ 1.) .* R_std_2H .+ 1. )
        # if u₀ < 0.001 && u⁺ > 0.001 # instead of u₀ == 0
        #     # Edge cases:
        #     x⁺ = sum(xin .* inflow) / sum(inflow)
        #     u_δ2H_SNOW = 1000 * (x⁺/(1-x⁺)/R_std_2H - 1)
        # elseif u⁺ > 0.001 # mm water = 1 µm water
        #     x⁺ = (x₀ * u₀ + p_DTP * (sum(xin .* inflow) - sum(x₀ .* outflow)) ) / u⁺
        #     u_δ2H_SNOW = 1000 * (x⁺/(1-x⁺)/R_std_2H - 1)
        # else
        #     x⁺ = NaN
        #     u_δ2H_SNOW = NaN
        # end
        # #TODO(bernhard): include fractionation due to evaporation
        # ## δ18O
        # δ₀ = u_δ18O_SNOW
        # x₀  = (δ₀/1000 + 1)*R_std_18O / (1 + (δ₀/1000 + 1)*R_std_18O)
        # δin = [u_δ18O_INTS,  p_δ18O_PREC(integrator.t)]
        # xin = (δin./1000. .+ 1).*R_std_18O ./ (1 .+ (δin./1000. .+ 1).*R_std_18O)
        # if u₀ < 0.001 && u⁺ > 0.001 # instead of u₀ == 0
        #     # Edge cases:
        #     x⁺ = sum(xin .* inflow) / sum(inflow)
        # elseif u⁺ > 0.001 # mm water = 1 µm water
        #     x⁺ = (x₀ * u₀ + p_DTP * (sum(xin .* inflow) - sum(x₀ .* outflow)) ) / u⁺
        # else
        #     x⁺ = NaN
        # end
        # u_δ18O_SNOW = 1000 * (x⁺/(1-x⁺)/R_std_18O - 1)
        # #TODO(bernhard): include fractionation due to evaporation
        # # For INTS
        # inflow = [aux_du_SINT]
        # outflow = [0]
        # u⁺ = u_INTS
        # u₀ = u_INTS - p_DTP * (sum(inflow) - sum(outflow)) # [mm]
        # E  = aux_du_ISVP
        # ## δ2H
        # δ₀ = u_δ2H_INTS
        # x₀  = (δ₀/1000 + 1)*R_std_2H / (1 + (δ₀/1000 + 1)*R_std_2H)
        # δin = p_δ2H_PREC(integrator.t)
        # xin = (δin./1000. .+ 1).*R_std_2H ./ (1 .+ (δin./1000. .+ 1).*R_std_2H)
        # if u₀ < 0.001 && u⁺ > 0.001 # instead of u₀ == 0
        #     # Edge cases:
        #     x⁺ = sum(xin .* inflow) / sum(inflow)
        # elseif u⁺ > 0.001 # mm water = 1 µm water
        #     x⁺ = (x₀ * u₀ + p_DTP * (sum(xin .* inflow) - sum(x₀ .* outflow)) ) / u⁺
        # else
        #     x⁺ = NaN
        # end
        # u_δ2H_INTS = 1000 * (x⁺/(1-x⁺)/R_std_2H - 1)
        # #TODO(bernhard): include fractionation due to evaporation
        # ## δ18O
        # δ₀ = u_δ18O_INTS
        # x₀  = (δ₀/1000 + 1)*R_std_18O / (1 + (δ₀/1000 + 1)*R_std_18O)
        # δin = p_δ18O_PREC(integrator.t)
        # xin = (δin./1000. .+ 1).*R_std_18O ./ (1 .+ (δin./1000. .+ 1).*R_std_18O)
        # if u₀ < 0.001 && u⁺ > 0.001 # instead of u₀ == 0
        #     # Edge cases:
        #     x⁺ = sum(xin .* inflow) / sum(inflow)
        # elseif u⁺ > 0.001 # mm water = 1 µm water
        #     x⁺ = (x₀ * u₀ + p_DTP * (sum(xin .* inflow) - sum(x₀ .* outflow)) ) / u⁺
        # else
        #     x⁺ = NaN
        # end
        # u_δ18O_INTS = 1000 * (x⁺/(1-x⁺)/R_std_18O - 1)
        # #TODO(bernhard): include fractionation due to evaporation
        # # For INTR
        # inflow = [aux_du_RINT]
        # outflow = [0]
        # u⁺ = u_INTR
        # u₀ = u_INTR - p_DTP * (sum(inflow) - sum(outflow)) # [mm]
        # E  = aux_du_IRVP
        # ## δ2H
        # δ₀ = u_δ2H_INTR
        # x₀  = (δ₀/1000 + 1)*R_std_2H / (1 + (δ₀/1000 + 1)*R_std_2H)
        # δin = p_δ2H_PREC(integrator.t)
        # xin = (δin./1000. .+ 1).*R_std_2H ./ (1 .+ (δin./1000. .+ 1).*R_std_2H)
        # if u₀ < 0.001 && u⁺ > 0.001 # instead of u₀ == 0
        #     # Edge cases:
        #     x⁺ = sum(xin .* inflow) / sum(inflow)
        # elseif u⁺ > 0.001 # mm water = 1 µm water
        #     x⁺ = (x₀ * u₀ + p_DTP * (sum(xin .* inflow) - sum(x₀ .* outflow)) ) / u⁺
        # else
        #     x⁺ = NaN
        # end
        # u_δ2H_INTR = 1000 * (x⁺/(1-x⁺)/R_std_2H - 1)
        # #TODO(bernhard): include fractionation due to evaporation
        # ## δ18O
        # δ₀ = u_δ18O_INTR
        # x₀  = (δ₀/1000 + 1)*R_std_18O / (1 + (δ₀/1000 + 1)*R_std_18O)
        # δin = p_δ18O_PREC(integrator.t)
        # xin = (δin./1000. .+ 1).*R_std_18O ./ (1 .+ (δin./1000. .+ 1).*R_std_18O)
        # if u₀ < 0.001 && u⁺ > 0.001 # instead of u₀ == 0
        #     # Edge cases:
        #     x⁺ = sum(xin .* inflow) / sum(inflow)
        # elseif u⁺ > 0.001 # mm water = 1 µm water
        #     x⁺ = (x₀ * u₀ + p_DTP * (sum(xin .* inflow) - sum(x₀ .* outflow)) ) / u⁺
        # else
        #     x⁺ = NaN
        # end
        # u_δ18O_INTR = 1000 * (x⁺/(1-x⁺)/R_std_18O - 1)
        # #TODO(bernhard): include fractionation due to evaporation
        # δ2H_SLFL   = p_δ2H_PREC(integrator.t)
        # δ18O_SLFL   = p_δ18O_PREC(integrator.t)
        # # END variant 1

        ####################################################################
        # Return results from callback

        # save intermediate results for use in ODE (function f()) or other callbacks
        integrator.p[3][1][1:2] .= [δ18O_SLFL, δ2H_SLFL]

        # update δ values of INTS, INTR, SNOW
        # do not update δ values of GWAT and SWATI (is done in f())

        # integrator.u[idx_u_scalar_isotopes_d18O[1]] = u_δ18O_GWAT
        # integrator.u[idx_u_scalar_isotopes_d2H[1] ] = u_δ2H_GWAT
        integrator.u[idx_u_scalar_isotopes_d18O[2]] = u_δ18O_INTS
        integrator.u[idx_u_scalar_isotopes_d2H[2] ] = u_δ2H_INTS
        integrator.u[idx_u_scalar_isotopes_d18O[3]] = u_δ18O_INTR
        integrator.u[idx_u_scalar_isotopes_d2H[3] ] = u_δ2H_INTR
        integrator.u[idx_u_scalar_isotopes_d18O[4]] = u_δ18O_SNOW
        integrator.u[idx_u_scalar_isotopes_d2H[4] ] = u_δ2H_SNOW
        # integrator.u[idx_u_vector_isotopes_d18O]    = u_δ18O_SWATI
        # integrator.u[idx_u_vector_isotopes_d2H]     = u_δ2H_SWATI
    end

    return nothing
    ##########################################
end



function LWFBrook90R_updateIsotopes_GWAT_SWAT!(u, t, integrator)

    use_numerical_solution = false
    simulate_isotopes = integrator.p[1][4][3]

    if simulate_isotopes
        idx_u_vector_amounts       = integrator.p[1][4][4]
        idx_u_vector_accumulators  = integrator.p[1][4][5]
        # idx_u_scalar_amounts     = integrator.p[1][4][6]

        u_GWAT     = integrator.u[1]
        u_INTS     = integrator.u[2]
        u_INTR     = integrator.u[3]
        u_SNOW     = integrator.u[4]
        u_CC       = integrator.u[5]
        u_SNOWLQ   = integrator.u[6]
        u_SWATI    = integrator.u[idx_u_vector_amounts]

        ####################################################################
        # 2) Update states of isotopic compositions of interception storages and snowpack:
        #       u_δ18O_INTS, u_δ18O_INTR, u_δ18O_SNOW
        #       u_δ2H_INTS, u_δ2H_INTR, u_δ2H_SNOW
        ############
        ### Compute parameters
        ## B) time dependent parameters
        p_DOY, p_MONTHN, p_SOLRAD, p_TMAX, p_TMIN, p_EA, p_UW, p_PREC,
            p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_AGE, p_RELDEN,
            p_δ18O_PREC, p_δ2H_PREC, ref_date = integrator.p[2]

        ## C) state dependent parameters or intermediate results:
        # These were computed in the callback and are kept constant in between two
        # callbacks.
        (δ18O_SLFL, δ2H_SLFL, p_fu_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP,
            p_fu_STHR, aux_du_RSNO, aux_du_SNVP,
            aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, u_SNOW_old) = integrator.p[3][1]
        aux_du_TRANI = integrator.p[3][2]
        du_NTFLI     = integrator.p[3][3][:,1]
        aux_du_VRFLI = integrator.p[3][3][:,2]
        aux_du_DSFLI = integrator.p[3][3][:,3]
        aux_du_INFLI = integrator.p[3][3][:,4]
        u_aux_WETNES = integrator.p[3][3][:,5]
        du_GWFL      = integrator.p[3][4][1]
        du_SEEP      = integrator.p[3][4][2]

        idx_u_scalar_isotopes_d18O = integrator.p[1][4][8]
        idx_u_vector_isotopes_d18O = integrator.p[1][4][9]
        idx_u_scalar_isotopes_d2H  = integrator.p[1][4][10]
        idx_u_vector_isotopes_d2H  = integrator.p[1][4][11]

        u_δ18O_GWAT = integrator.u[idx_u_scalar_isotopes_d18O[1]]
        u_δ2H_GWAT  = integrator.u[idx_u_scalar_isotopes_d2H[1] ]
        # u_δ18O_INTS = integrator.u[idx_u_scalar_isotopes_d18O[2]]
        # u_δ2H_INTS  = integrator.u[idx_u_scalar_isotopes_d2H[2] ]
        # u_δ18O_INTR = integrator.u[idx_u_scalar_isotopes_d18O[3]]
        # u_δ2H_INTR  = integrator.u[idx_u_scalar_isotopes_d2H[3] ]
        # u_δ18O_SNOW = integrator.u[idx_u_scalar_isotopes_d18O[4]]
        # u_δ2H_SNOW  = integrator.u[idx_u_scalar_isotopes_d2H[4] ]
        u_δ18O_SWATI = integrator.u[idx_u_vector_isotopes_d18O]
        u_δ2H_SWATI  = integrator.u[idx_u_vector_isotopes_d2H]

        EffectiveDiffusivity_18O = 0  # TODO(bernhard): compute correct values using eq 33 Zhou-2021-Environ_Model_Softw.pdf
        EffectiveDiffusivity_2H  = 0  # TODO(bernhard): compute correct values using eq 33 Zhou-2021-Environ_Model_Softw.pdf
        δ18O_INFLI = δ18O_SLFL
        δ2H_INFLI  = δ2H_SLFL

        if use_numerical_solution
            du_δ18O_GWAT, du_δ2H_GWAT, du_δ18O_SWATI, du_δ2H_SWATI =
                compute_isotope_du_GWAT_SWATI(
                    # for GWAT:
                    u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT,
                    # for SWATI:
                    du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI,  # (non-fractionating)
                    aux_du_SLVP, p_fu_TADTM, p_EA(t), p_δ2H_PREC(t), p_δ18O_PREC(t), u_aux_WETNES, # (fractionating)
                    u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, 0, 0) #EffectiveDiffusivity_18O, EffectiveDiffusivity_2H)

            # update δ values of GWAT and SWATI
            # du[idx_u_scalar_isotopes_d18O[1]]   = du_δ18O_GWAT     #TODO(bernhard)
            # du[idx_u_scalar_isotopes_d2H[1] ]   = du_δ2H_GWAT      #TODO(bernhard)
            # du[idx_u_vector_isotopes_d18O]      = du_δ18O_SWATI    #TODO(bernhard)
            # du[idx_u_vector_isotopes_d2H]       = du_δ2H_SWATI     #TODO(bernhard)

            # Since we update this in the callback we can't overwrite `du` and let DiffEq.jl do
            # the work, so we need to update `u` instead of `du`
            # Using a simple forward euler update:
            u[idx_u_scalar_isotopes_d18O[1]]   += integrator.dt * du_δ18O_GWAT     #TODO(bernhard)
            u[idx_u_scalar_isotopes_d2H[1] ]   += integrator.dt * du_δ2H_GWAT      #TODO(bernhard)
            u[idx_u_vector_isotopes_d18O]      .+= integrator.dt * du_δ18O_SWATI    #TODO(bernhard)
            u[idx_u_vector_isotopes_d2H]       .+= integrator.dt * du_δ2H_SWATI     #TODO(bernhard)
        else # using analytical solution
            u_δ18O_GWAT, u_δ2H_GWAT, u_δ18O_SWATI, u_δ2H_SWATI =
                compute_isotope_u_GWAT_SWATI(integrator,
                    # for GWAT:
                    u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT, du_GWFL, du_SEEP,
                    # for SWATI:
                    du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI,  # (non-fractionating)
                    aux_du_SLVP, p_fu_TADTM, p_EA(t), p_δ2H_PREC(t), p_δ18O_PREC(t), u_aux_WETNES, # (fractionating)
                    u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, 0, 0) #EffectiveDiffusivity_18O, EffectiveDiffusivity_2H)

            u[idx_u_scalar_isotopes_d18O[1]]   = u_δ18O_GWAT
            u[idx_u_scalar_isotopes_d2H[1] ]   = u_δ2H_GWAT
            u[idx_u_vector_isotopes_d18O]      = u_δ18O_SWATI
            u[idx_u_vector_isotopes_d2H]       = u_δ2H_SWATI
            #TODO(bernhard): this update generates an error message with adaptive solvers:
            # Warning: dt <= dtmin. Aborting. There is either an error in your model specification or the true solution is unstable.
            # └ @ SciMLBase ~/.julia/packages/SciMLBase/BoNUy/src/integrator_interface.jl:345
            # TODO(bernhard): we can force the solving by using force_dtmin = true
        end
    end

    return nothing
end




function LWFBrook90R_check_balance_errors!(integrator)

    # Compute daily water balance errors

    (NLAYER, FLAG_MualVanGen, compute_intermediate_quantities, Reset,
    p_DTP, p_NPINT,

    _, _, _, _, _, _,
    _, _, _, _,
    _, _, _, _,
    _, _,

    _) = integrator.p[1][2]

    idx_u_vector_amounts       = integrator.p[1][4][4]
    idx_u_vector_accumulators  = integrator.p[1][4][5]

    u_GWAT     = integrator.u[2]
    u_INTS     = integrator.u[2]
    u_INTR     = integrator.u[3]
    u_SNOW     = integrator.u[4]
    # u_CC       = integrator.u[5]
    # u_SNOWLQ   = integrator.u[6]
    u_SWATI    = integrator.u[idx_u_vector_amounts]

    if compute_intermediate_quantities
        # a) Get change in total storages
        # Get old total water volumes and compute new total water volumes
        old_SWAT       = integrator.uprev[idx_u_vector_accumulators[28]]
        old_totalWATER = integrator.uprev[idx_u_vector_accumulators[29]]
        new_SWAT       = sum(u_SWATI) # total soil water in all layers, mm
        new_totalWATER = u_INTR + u_INTS + u_SNOW + new_SWAT + u_GWAT # total water in all compartments, mm

        # Save current totals into caches to check error for next period
        integrator.u[idx_u_vector_accumulators[28]] = new_SWAT
        integrator.u[idx_u_vector_accumulators[29]] = new_totalWATER

        # b) Get fluxes that caused change in total storages
        # Across outer boundaries of above and belowground system:
        idx_cum_d_prec = idx_u_vector_accumulators[ 1]
        idx_cum_d_evap = idx_u_vector_accumulators[ 9]
        idx_cum_flow   = idx_u_vector_accumulators[18]
        idx_cum_seep   = idx_u_vector_accumulators[19]

        # Across outer boundaries of belowground system (SWATI):
        # SLFL  # idx_u_vector_accumulators[21] # cum_d_slfl
        # BYFL  # idx_u_vector_accumulators[22] # cum_d_byfl
        # INFLI # as sum(INFLI) =  SLFL - sum(BYFL)
        idx_cum_slfl = idx_u_vector_accumulators[21]
        idx_cum_byfl = idx_u_vector_accumulators[22]
        idx_cum_dsfl  = idx_u_vector_accumulators[23]# cum_d_dsfl
        idx_cum_vrfln = idx_u_vector_accumulators[25] # cum_d_vrfln
        idx_cum_tran  = idx_u_vector_accumulators[10] # cum_d_tran
        idx_cum_slvp  = idx_u_vector_accumulators[13] # cum_d_slvp


        # c) Compute balance error
        # BALERD_SWAT  = old_SWAT - new_SWAT + INFLI - DSFLI - VRFL(NLAYER) - TRANI - SLVP
        BALERD_SWAT  = old_SWAT - new_SWAT +
            (integrator.u[idx_cum_slfl] - integrator.u[idx_cum_byfl]) - integrator.u[idx_cum_dsfl] - integrator.u[idx_cum_vrfln] - integrator.u[idx_cum_tran] - integrator.u[idx_cum_slvp]
        # BALERD_total = old_totalWATER - new_totalWATER + PRECD - EVAPD - FLOWD - SEEPD
        BALERD_total = old_totalWATER - new_totalWATER +
            integrator.u[idx_cum_d_prec] - integrator.u[idx_cum_d_evap] - integrator.u[idx_cum_flow] - integrator.u[idx_cum_seep]

        # d) Store balance errors into state vector
        integrator.u[idx_u_vector_accumulators[30]] = BALERD_SWAT
        integrator.u[idx_u_vector_accumulators[31]] = BALERD_total
    end

    return nothing
end
