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

    cb_SWAT_GWAT_deltas = FunctionCallingCallback(
                                LWFBrook90R_updateIsotopes_GWAT_SWAT_AdvecDiff!;
                                # LWFBrook90R_updateIsotopes_GWAT_SWAT!;
                                func_everystep = true,
                                func_start = false);

    #TODO(bernhard) Implement swchek from LWFBrook90 as ContinuousCallback
    # swcheck_cb = ContinuousCallback()

    cb_set = CallbackSet(
        # 1) Belowground:
        # Check and correct flux (Richards equation)
        # # not implemented swcheck_cb,
        # # not implemented: cb_LaiRichardsCorrectorStep,
        cb_check_balance_errors,
        # Continuous update of belowground temperature
        # # cb_Temp not implemented
        # Continuous update of belowground concentrations (Adv.-Disp. equation)
        cb_SWAT_GWAT_deltas,

        # 2) Aboveground
        # Daily update of aboveground storages and concentrations
        cb_INTS_INTR_SNOW_amounts,
        cb_INTS_INTR_SNOW_deltas
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
    @unpack p_soil = integrator.p;
    @unpack NLAYER, FLAG_MualVanGen, compute_intermediate_quantities, Reset, p_DTP, p_NPINT = integrator.p;
    @unpack p_LAT, p_ESLOPE, p_L1, p_L2,
        p_SNODEN, p_MXRTLN, p_MXKPL,
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
        p_DURATN, p_MAXLQF, p_GRDMLT,

        p_VXYLEM, p_DISPERSIVITY = integrator.p;

    ## B) time dependent parameters
    @unpack p_DOY, p_MONTHN, p_GLOBRAD, p_TMAX, p_TMIN, p_VAPPRES, p_WIND, p_PREC,
        p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_fT_RELDEN,
        p_δ18O_PREC, p_δ2H_PREC, REFERENCE_DATE = integrator.p;

    ## C) state dependent parameters:
    # Calculate parameters:
    #  - solar parameters depending on DOY
    #  - canopy parameters depending on DOY and u_SNOW
    #  - roughness parameters depending on u_SNOW
    #  - plant resistance components depending on u_SNOW
    #  - weather data depending on DOY and u_SNOW
    #  - fraction of precipitation as snowfall depending on DOY
    #  - snowpack temperature, potential snow evaporation and soil evaporation resistance depending on u_SNOW
    @unpack p_fT_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP,
        p_fu_STHR, aux_du_RSNO, aux_du_SNVP,
        aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, u_SNOW_old,
        aux_du_TRANI = integrator.p;
    @unpack p_fu_δ18O_SLFL, p_fu_δ2H_SLFL = integrator.p;

    # integrator.p.aux_du_TRANI
    # Parse states
    u_INTS     = integrator.u.INTS.mm
    u_INTR     = integrator.u.INTR.mm
    u_SNOW     = integrator.u.SNOW.mm
    u_CC       = integrator.u.CC.MJm2
    u_SNOWLQ   = integrator.u.SNOWLQ.mm
    u_SWATI    = integrator.u.SWATI.mm

    LWFBrook90.KPT.SWCHEK!(u_SWATI, p_soil.p_SWATMAX, integrator.t)

    # Derive (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) from u_SWATI
    (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        LWFBrook90.KPT.derive_auxiliary_SOILVAR(u_SWATI, p_soil)

    # MSBSETVARS: A) get: 1) sunshine durations, 2) SFAL, 3) plant resistance (TODO: could be done before simulation)
    # MSBSETVARS: B) get_windspeed_from_canopy_and_snowpack: (p_fu_UADTM, p_fu_UANTM) = f(...)
    # MSBSETVARS: C) get_snowpack_evaporation_from_canopy,snowpack,atmosphere,: (p_fu_UADTM, p_fu_UANTM) = f(...)
    p_fT_DAYLEN, p_fT_I0HDAY, p_fT_SLFDAY, p_fu_HEIGHT, p_fu_LAI, p_fu_SAI,
        p_fu_Z0GS, p_fu_Z0C, p_fu_DISPC, p_fu_Z0, p_fu_DISP, p_fu_ZA,
        p_fT_RXYLEM, p_fT_RROOTI, p_fT_ALPHA,
        p_fu_SHEAT, p_fT_SOLRADC, p_fT_TA, p_fT_TADTM[1], p_fT_TANTM, p_fu_UADTM, p_fu_UANTM,
        p_fT_SNOFRC, p_fu_TSNOW, p_fu_PSNVP, p_fu_ALBEDO,p_fu_RSS, p_fu_SNOEN =
            MSBSETVARS(FLAG_MualVanGen, NLAYER, p_soil,
                       # for SUNDS:
                       p_LAT, p_ESLOPE, p_DOY(integrator.t), p_L1, p_L2,
                       # for CANOPY:
                       p_HEIGHT(integrator.t), p_LAI(integrator.t), p_SAI(integrator.t), u_SNOW, p_SNODEN, p_MXRTLN, p_MXKPL, p_DENSEF(integrator.t),
                       #
                       p_Z0S, p_Z0G,
                       # for ROUGH:
                       p_ZMINH, p_CZS, p_CZR, p_HS, p_HR, p_LPC,
                       # for PLNTRES:
                       p_fT_RELDEN.(integrator.t, 1:NLAYER), p_RTRAD, p_FXYLEM,
                       # for WEATHER:
                       p_TMAX(integrator.t), p_TMIN(integrator.t), p_VAPPRES(integrator.t), p_WIND(integrator.t), p_WNDRAT, p_FETCH, p_Z0W, p_ZW, p_GLOBRAD(integrator.t),
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
    (p_fu_PTR, p_fu_GER, p_fu_PIR, p_fu_GIR, p_fu_ATRI) =
        MSBDAYNIGHT(p_fT_SLFDAY, p_fT_SOLRADC, p_WTOMJ, p_fT_DAYLEN, p_fT_TADTM[1], p_fu_UADTM, p_fT_TANTM, p_fu_UANTM,
                    p_fT_I0HDAY,
                    # for AVAILEN:
                    p_fu_ALBEDO, p_C1, p_C2, p_C3, p_VAPPRES(integrator.t), p_fu_SHEAT, p_CR, p_fu_LAI, p_fu_SAI,
                    # for SWGRA:
                    p_fu_ZA, p_fu_HEIGHT, p_fu_Z0, p_fu_DISP, p_fu_Z0C, p_fu_DISPC, p_fu_Z0GS, p_LWIDTH, p_RHOTP, p_NN,
                    # for SRSC:
                    p_fT_TA, p_GLMIN, p_GLMAX, p_R5, p_CVPD, p_RM, p_TL, p_T1, p_T2, p_TH,
                    # for SWPE:
                    p_fu_RSS,
                    # for TBYLAYER:
                    p_fT_ALPHA, p_fu_KK, p_fT_RROOTI, p_fT_RXYLEM, u_aux_PSITI, NLAYER, p_PSICR, NOOUTF)
                    # 0.000012 seconds (28 allocations: 1.938 KiB)
    # Combine day and night rates to average daily rate
    (p_fu_PTRAN, p_fu_GEVP, p_fu_PINT, p_fu_GIVP, aux_du_TRANI[:]) =
        MSBDAYNIGHT_postprocess(NLAYER, p_fu_PTR, p_fu_GER, p_fu_PIR, p_fu_GIR, p_fu_ATRI, p_fT_DAYLEN)
    #* * * * * * * *  E N D   D A Y - N I G H T   L O O P  * * * * * * * * * *
    ####################################################################
    # 1) Update snow accumulation/melt: u_SNOW, u_CC, u_SNOWLQ
    #    and compute fluxes to/from interception storage
    u_SNOW_old[1] = u_SNOW
    (# compute some fluxes as intermediate results to be used in RHS function f:
    p_fT_SFAL, p_fT_RFAL, p_fu_RNET[1], p_fu_PTRAN,
    # compute changes in soil water storage:
    aux_du_TRANI[:], aux_du_SLVP[1],
    # compute change in interception storage:
    aux_du_SINT[1], aux_du_ISVP[1], aux_du_RINT[1], aux_du_IRVP[1],
    # compute change in snow storage:
    aux_du_RSNO[1], aux_du_SNVP[1], aux_du_SMLT[1], p_fu_STHR[1],
    # compute updated states:
    u_SNOW_MSBupdate, u_CC, u_SNOWLQ) =
        MSBPREINT(p_PREC(integrator.t), p_DTP, p_fT_SNOFRC, p_NPINT, p_fu_PINT, p_fT_TA,
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
    u_INTS = u_INTS + (aux_du_SINT[1] - aux_du_ISVP[1]) * p_DTP
    u_INTR = u_INTR + (aux_du_RINT[1] - aux_du_IRVP[1]) * p_DTP
    u_SNOW = u_SNOW_MSBupdate

    ####################################################################
    # 3) Update soil water using substeps smaller than precipitation
    #    interval. After each substep update: u_GWAT, u_SWATI
    # ... NOTE: this is not done in the callback function but in the main ODE solver

    ####################################################################
    # Return results from callback
    # update INTS
    integrator.u.INTS.mm = u_INTS

    # update INTR
    integrator.u.INTR.mm   = u_INTR

    # update SNOW, CC, SNOWLQ
    integrator.u.SNOW.mm   = u_SNOW
    integrator.u.CC.MJm2   = u_CC
    integrator.u.SNOWLQ.mm = u_SNOWLQ

    # save intermediate results for use in ODE (function f()) or other callbacks
    # integrator.p[3][1] .= [p_fu_δ18O_SLFL, p_fu_δ2H_SLFL,
    #                        p_fT_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP,
    #                        p_fu_STHR, aux_du_RSNO, aux_du_SNVP,
    #                        aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, u_SNOW_old]
    # integrator.p[3][1] .= [NaN, NaN,        # NaNs are directly after overwritten by callback for isotopes
    #                         p_fT_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP,
    #                         p_fu_STHR, aux_du_RSNO, aux_du_SNVP,
    #                         aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, u_SNOW_old]
    # integrator.p[3][2] .= aux_du_TRANI
    integrator.p.p_fu_δ18O_SLFL .= NaN # NaNs are directly after overwritten by callback for isotopes
    integrator.p.p_fu_δ2H_SLFL  .= NaN # NaNs are directly after overwritten by callback for isotopes
    ##########################################
    # Accumulate flows to compute daily sums
    # Note that below state variables serve only as accumulator but do not affect
    # the evolution of the system.

    integrator.u.TRANI.mmday .= p_DTP .* aux_du_TRANI
    integrator.u.RWU.mmday    = p_DTP  * sum(aux_du_TRANI)
    if compute_intermediate_quantities

        # 1) Either set daily sum if rate is constant throughout precipitation interval: p_DTP*(...)
        # 2) or then set daily sum to zero and use ODE to accumulate flow.
        integrator.u.accum.cum_d_prec       = p_DTP * (p_fT_RFAL + p_fT_SFAL)                 # RFALD + SFALD        # cum_d_prec
        integrator.u.accum.cum_d_rfal       = p_DTP * (p_fT_RFAL)                                                    # cum_d_rfal
        integrator.u.accum.cum_d_sfal       = p_DTP * (p_fT_SFAL)                                                    # cum_d_sfal
        integrator.u.accum.cum_d_rint       = p_DTP * (aux_du_RINT[1])                                                  # cum_d_rint
        integrator.u.accum.cum_d_sint       = p_DTP * (aux_du_SINT[1])                                                  # cum_d_sint
        integrator.u.accum.cum_d_rsno       = p_DTP * (aux_du_RSNO[1])                                                  # cum_d_rsno
        integrator.u.accum.cum_d_rnet       = p_DTP * (p_fT_RFAL - aux_du_RINT[1] - aux_du_RSNO[1]) # cum_d_RTHR - RSNOD   # cum_d_rnet
        integrator.u.accum.cum_d_smlt       = p_DTP * (aux_du_SMLT[1])                                                  # cum_d_smlt
        integrator.u.accum.cum_d_evap       = p_DTP * (aux_du_IRVP[1] + aux_du_ISVP[1] + aux_du_SNVP[1] + aux_du_SLVP[1] + sum(aux_du_TRANI))  # cum_d_evap
        integrator.u.accum.cum_d_tran       = p_DTP * (sum(aux_du_TRANI))                                                          # cum_d_tran
        integrator.u.accum.cum_d_irvp       = p_DTP * (aux_du_IRVP[1])                                                                # cum_d_irvp
        integrator.u.accum.cum_d_isvp       = p_DTP * (aux_du_ISVP[1])                                                                # cum_d_isvp
        integrator.u.accum.cum_d_slvp       = p_DTP * (aux_du_SLVP[1])                                                                # cum_d_slvp
        integrator.u.accum.cum_d_snvp       = p_DTP * (aux_du_SNVP[1])                                                                # cum_d_snvp
        integrator.u.accum.cum_d_pint       = p_DTP * (p_fu_PINT)                                                                  # cum_d_pint
        integrator.u.accum.cum_d_ptran      = p_DTP * (p_fu_PTRAN)                                                                 # cum_d_ptran
        integrator.u.accum.cum_d_pslvp      = 0 # p_DTP * (p_fu_PSLVP)                                                             # cum_d_pslvp # Deactivated as p_fu_PSLVP is never used
        integrator.u.accum.flow             = 0 # flow,  is computed in ODE
        integrator.u.accum.seep             = 0 # seep,  is computed in ODE
        integrator.u.accum.srfl             = 0 # srfl,  is computed in ODE
        integrator.u.accum.slfl             = 0 # slfl,  is computed in ODE
        integrator.u.accum.byfl             = 0 # byfl,  is computed in ODE
        integrator.u.accum.dsfl             = 0 # dsfl,  is computed in ODE
        integrator.u.accum.gwfl             = 0 # gwfl,  is computed in ODE
        integrator.u.accum.vrfln            = 0 # vrfln, is computed in ODE
        integrator.u.accum.cum_d_rthr     = p_DTP*(p_fT_RFAL - aux_du_RINT[1]) # cum_d_rthr
        integrator.u.accum.cum_d_sthr     = p_DTP*(p_fT_SFAL - aux_du_SINT[1]) # cum_d_sthr
        # integrator.u.accum.ε_prev_t          = t              # Update in separate daily callback
        # integrator.u.accum.ε_prev_totalSWAT  = new_totalSWAT  # Update in separate daily callback
        # integrator.u.accum.ε_prev_totalWATER = new_totalWATER # Update in separate daily callback
        # integrator.u.accum.BALERD_SWAT       = BALERD_SWAT    # Update in separate daily callback
        # integrator.u.accum.BALERD_total      = BALERD_total   # Update in separate daily callback

        # TODO(bernhard): use SavingCallback() for all quantities that have u=... and du=0
        #                 only keep du=... for quantities for which we compute cumulative sums
    end

    return nothing
    ##########################################
end

function LWFBrook90R_updateIsotopes_INTS_INTR_SNOW!(integrator)


    @unpack simulate_isotopes = integrator.p

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

        # Parse states
        u_INTS     = integrator.u.INTS.mm
        u_INTR     = integrator.u.INTR.mm
        u_SNOW     = integrator.u.SNOW.mm
        # u_CC       = integrator.u.CC.MJm2
        # u_SNOWLQ   = integrator.u.SNOWLQ.mm
        # u_SWATI    = integrator.u.SWATI.mm

        ####################################################################
        # 2) Update states of isotopic compositions of interception storages and snowpack:
        #       u_δ18O_INTS, u_δ18O_INTR, u_δ18O_SNOW
        #       u_δ2H_INTS, u_δ2H_INTR, u_δ2H_SNOW
        ############
        ### Compute parameters
        ## A) constant parameters:
        @unpack p_DTP = integrator.p # TODO: is p_DTP still used??

        ## B) time dependent parameters
        # @unpack p_DOY, p_MONTHN, p_GLOBRAD, p_TMAX, p_TMIN, p_VAPPRES, p_WIND, p_PREC,
        #     p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_fT_RELDEN,
        #     p_δ18O_PREC, p_δ2H_PREC, REFERENCE_DATE = integrator.p
        @unpack p_δ18O_PREC, p_δ2H_PREC = integrator.p

        ## C) state dependent parameters or intermediate results:
        # These were computed in the callback and are kept constant in between two
        # callbacks.
        @unpack p_fT_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP,
            p_fu_STHR, aux_du_RSNO, aux_du_SNVP,
            aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, u_SNOW_old = integrator.p
        @unpack aux_du_TRANI = integrator.p
        @unpack p_fu_δ18O_SLFL, p_fu_δ2H_SLFL = integrator.p;

        # u_δ18O_GWAT = integrator.u.GWAT.d18O
        # u_δ2H_GWAT  = integrator.u.GWAT.d2H
        u_δ18O_INTS = integrator.u.INTS.d18O
        u_δ2H_INTS  = integrator.u.INTS.d2H
        u_δ18O_INTR = integrator.u.INTR.d18O
        u_δ2H_INTR  = integrator.u.INTR.d2H
        u_δ18O_SNOW = integrator.u.SNOW.d18O
        u_δ2H_SNOW  = integrator.u.SNOW.d2H
        # u_δ18O_SWATI = integrator.u.SWATI.d18O
        # u_δ2H_SWATI  = integrator.u.SWATI.d18O
        # u_δ18O_RWU = integrator.u.RWU.d18O
        # u_δ2H_RWU  = integrator.u.RWU.d2H
        # u_δ¹⁸O_Xylem = integrator.u.XYLEM.d18O
        # u_δ²H_Xylem  = integrator.u.XYLEM.d2H

        # Variant 1) works but really much slower
        # p_fu_δ18O_SLFL, p_fu_δ2H_SLFL,
        # _,                u_δ18O_INTS, u_δ2H_INTS,
        # _,            u_δ18O_INTR, u_δ2H_INTR,
        # u_SNOW_iso_update, u_δ18O_SNOW, u_δ2H_SNOW =
        #     compute_isotope_U_of_INTS_INTR_SNOW_and_SLFL(
        # (p_fu_δ18O_SLFL, p_fu_δ2H_SLFL,
        # u_δ18O_INTS, u_δ2H_INTS,
        # u_δ18O_INTR, u_δ2H_INTR,
        # u_δ18O_SNOW, u_δ2H_SNOW) =
        #     compute_isotope_U_of_INTS_INTR_SNOW_and_SLFL(
        #         u_INTS, u_δ18O_INTS, u_δ2H_INTS, u_INTR, u_δ18O_INTR, u_δ2H_INTR, u_SNOW, u_δ18O_SNOW, u_δ2H_SNOW,
        #         p_δ2H_PREC(integrator.t), p_δ18O_PREC(integrator.t), p_fT_TADTM[1], p_VAPPRES(integrator.t),
        #         # for INTS (in: SINT; out: ISVP):
        #         aux_du_SINT[1], aux_du_ISVP[1], p_DTP,
        #         # for INTR (in: RINT; out: IRVP):
        #         aux_du_RINT[1], aux_du_IRVP[1],
        #         # for SNOW (in: STHR, RSNO (both δ_PREC); out: SMLT, SNVP (δ_SNOW and fractionated)):
        #         NaN, p_fu_STHR[1], aux_du_RSNO[1], aux_du_SMLT[1], aux_du_SNVP[1],
        #         # to compute isotopic signature of soil infiltration: SLFL
        #         p_fu_RNET[1])
        # # END variant 1

        # Variant 2) resulting in no effect (but is fast, so not generating any allocation problems)
        p_fu_δ2H_SLFL[1]   = p_δ2H_PREC(integrator.t)
        u_δ2H_INTS = p_δ2H_PREC(integrator.t)
        u_δ2H_INTR = p_δ2H_PREC(integrator.t)
        u_δ2H_SNOW = p_δ2H_PREC(integrator.t)
        p_fu_δ18O_SLFL[1]   = p_δ18O_PREC(integrator.t)
        u_δ18O_INTS = p_δ18O_PREC(integrator.t)
        u_δ18O_INTR = p_δ18O_PREC(integrator.t)
        u_δ18O_SNOW = p_δ18O_PREC(integrator.t)
        # # END variant 2

        # # Variant 3) trying to unwrap variant 1) by reducing allocations. Currently not working and also not faster than 1)
        # R_std_2H  = LWFBrook90.ISO.R_VSMOW²H
        # R_std_18O = LWFBrook90.ISO.R_VSMOW¹⁸O
        # # For SNOW ( do this before INTS as it uses u_δ2H_INTS as input)
        # inflow = [p_fu_STHR[1]; aux_du_RSNO[1]]
        # outflow = [aux_du_SMLT[1]]
        # u⁺ = u_SNOW
        # u₀ = u_INTR - p_DTP * (sum(inflow) - sum(outflow)) # [mm]
        # E  = aux_du_SNVP[1]
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
        # inflow = [aux_du_SINT[1]]
        # outflow = [0]
        # u⁺ = u_INTS
        # u₀ = u_INTS - p_DTP * (sum(inflow) - sum(outflow)) # [mm]
        # E  = aux_du_ISVP[1]
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
        # inflow = [aux_du_RINT[1]]
        # outflow = [0]
        # u⁺ = u_INTR
        # u₀ = u_INTR - p_DTP * (sum(inflow) - sum(outflow)) # [mm]
        # E  = aux_du_IRVP[1]
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
        # p_fu_δ2H_SLFL   = p_δ2H_PREC(integrator.t)
        # p_fu_δ18O_SLFL   = p_δ18O_PREC(integrator.t)
        # # END variant 1

        ####################################################################
        # Return results from callback

        # save intermediate results for use in ODE (function f()) or other callbacks
        # integrator.p[3][1][1:2] .= [p_fu_δ18O_SLFL, p_fu_δ2H_SLFL]
        # not needed with @unpack

        # update δ values of INTS, INTR, SNOW
        # do not update δ values of GWAT and SWATI (is done in f())

        # integrator.u.GWAT.d18O = u_δ18O_GWAT
        # integrator.u.GWAT.d2H  = u_δ2H_GWAT
        integrator.u.INTS.d18O = u_δ18O_INTS
        integrator.u.INTS.d2H  = u_δ2H_INTS
        integrator.u.INTR.d18O = u_δ18O_INTR
        integrator.u.INTR.d2H  = u_δ2H_INTR
        integrator.u.SNOW.d18O = u_δ18O_SNOW
        integrator.u.SNOW.d2H  = u_δ2H_SNOW
        # integrator.u[idx_u_vector_isotopes_d18O]    = u_δ18O_SWATI
        # integrator.u[idx_u_vector_isotopes_d2H]     = u_δ2H_SWATI
    end

    return nothing
    ##########################################
end



function LWFBrook90R_updateIsotopes_GWAT_SWAT!(u, t, integrator)

    # Define flag to switch between different methods
    use_method = ["numerical_ForwardEuler",
                  "analytical"][1]
    # Note that the method "numerical_ForwardEuler" has been implemented more efficiently and
    # with diffusive fluxes in LWFBrook90R_updateIsotopes_GWAT_SWAT_AdvecDiff!().
    # The code below relating to "numerical_ForwardEuler" is redundant... the code relating
    # to analytical however might still be used.

    @unpack simulate_isotopes = integrator.p
    if simulate_isotopes
        u_GWAT     = integrator.u.GWAT.mm
        # u_INTS     = integrator.u.INTS.mm
        # u_INTR     = integrator.u.INTR.mm
        # u_SNOW     = integrator.u.SNOW.mm
        # u_CC       = integrator.u.CC.MJm2
        # u_SNOWLQ   = integrator.u.SNOWLQ.mm
        u_SWATI    = integrator.u.SWATI.mm

        ####################################################################
        # 2) Update states of isotopic compositions of interception storages and snowpack:
        #       u_δ18O_INTS, u_δ18O_INTR, u_δ18O_SNOW
        #       u_δ2H_INTS, u_δ2H_INTR, u_δ2H_SNOW
        ############
        ### Compute parameters
        ## B) time dependent parameters
        # p_DOY, p_MONTHN, p_GLOBRAD, p_TMAX, p_TMIN, p_VAPPRES, p_WIND, p_PREC,
        #     p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_fT_RELDEN,
        #     p_δ18O_PREC, p_δ2H_PREC, REFERENCE_DATE = integrator.p[2]
        @unpack p_DOY, p_MONTHN, p_GLOBRAD, p_TMAX, p_TMIN, p_VAPPRES, p_WIND, p_PREC,
            p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_fT_RELDEN,
            p_δ18O_PREC, p_δ2H_PREC, REFERENCE_DATE = integrator.p

        ## C) state dependent parameters or intermediate results:
        # These were computed in the callback and are kept constant in between two
        # callbacks.
        # TODO(bernharf): unpack these....
        # (p_fu_δ18O_SLFL, p_fu_δ2H_SLFL, p_fT_TADTM[1], p_fu_RNET[1], aux_du_SMLT[1], aux_du_SLVP[1],
        #     p_fu_STHR[1], aux_du_RSNO[1], aux_du_SNVP[1],
        #     aux_du_SINT[1], aux_du_ISVP[1], aux_du_RINT[1], aux_du_IRVP[1], u_SNOW_old[1]) = integrator.p[3][1]
        # aux_du_TRANI = integrator.p[3][2]
        # du_NTFLI     = integrator.p[3][3][:,1]
        # aux_du_VRFLI = integrator.p[3][3][:,2]
        # aux_du_DSFLI = integrator.p[3][3][:,3]
        # aux_du_INFLI = integrator.p[3][3][:,4]
        # u_aux_WETNES = integrator.p[3][3][:,5]
        # du_GWFL      = integrator.p[3][4][1]
        # du_SEEP      = integrator.p[3][4][2]
        @unpack p_fu_δ18O_SLFL, p_fu_δ2H_SLFL, p_fT_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP,
            p_fu_STHR, aux_du_RSNO, aux_du_SNVP,
            aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, u_SNOW_old = integrator.p
        @unpack aux_du_TRANI = integrator.p         # mm/day
        @unpack aux_du_DSFLI, aux_du_VRFLI, aux_du_INFLI, du_NTFLI = integrator.p # all in mm/day
        @unpack u_aux_WETNES = integrator.p         # all in mm/day
        @unpack du_GWFL, du_SEEP = integrator.p     # all in mm/day

        u_δ18O_GWAT  = integrator.u.GWAT.d18O
        u_δ2H_GWAT   = integrator.u.GWAT.d2H
        u_δ18O_SWATI = integrator.u.SWATI.d18O
        u_δ2H_SWATI  = integrator.u.SWATI.d2H

        EffectiveDiffusivity_18O = 0  # TODO(bernhard): compute correct values using eq 33 Zhou-2021-Environ_Model_Softw.pdf
        EffectiveDiffusivity_2H  = 0  # TODO(bernhard): compute correct values using eq 33 Zhou-2021-Environ_Model_Softw.pdf
        δ18O_INFLI = p_fu_δ18O_SLFL
        δ2H_INFLI  = p_fu_δ2H_SLFL

        if use_method == "numerical_ForwardEuler"
            du_δ18O_GWAT, du_δ2H_GWAT, du_δ18O_SWATI, du_δ2H_SWATI =
                compute_isotope_du_GWAT_SWATI(
                    # for GWAT:
                    u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT,
                    # for SWATI:
                    du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI,  # (non-fractionating)
                    aux_du_SLVP[1], p_fT_TADTM[1], p_VAPPRES(integrator.t), p_δ2H_PREC(integrator.t), p_δ18O_PREC(integrator.t), u_aux_WETNES, # (fractionating)
                    u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, 0, 0) #EffectiveDiffusivity_18O, EffectiveDiffusivity_2H)

            # update δ values of GWAT and SWATI
            # du[integrator.p[1][4].row_idx_scalars.GWAT, integrator.p[1][4].col_idx_d18O]   = du_δ18O_GWAT     #TODO(bernhard)
            # du[integrator.p[1][4].row_idx_scalars.GWAT, integrator.p[1][4].col_idx_d2H ]   = du_δ2H_GWAT      #TODO(bernhard)
            # du[idx_u_vector_isotopes_d18O]      = du_δ18O_SWATI    #TODO(bernhard)
            # du[idx_u_vector_isotopes_d2H]       = du_δ2H_SWATI     #TODO(bernhard)

            # Since we update this in the callback we can't overwrite `du` and let DiffEq.jl do
            # the work, so we need to update `u` instead of `du`
            # Using a simple forward euler update:
            integrator.u.GWAT.d18O   += integrator.dt * du_δ18O_GWAT     #TODO(bernhard)
            integrator.u.GWAT.d2H    += integrator.dt * du_δ2H_GWAT      #TODO(bernhard)
            integrator.u.SWATI.d18O  .+= integrator.dt * du_δ18O_SWATI    #TODO(bernhard)
            integrator.u.SWATI.d2H   .+= integrator.dt * du_δ2H_SWATI     #TODO(bernhard)
        elseif use_method == "analytical"
            #TODO(bernhard): this update generates an error message with adaptive solvers:
            # Warning: dt <= dtmin. Aborting. There is either an error in your model specification or the true solution is unstable.
            # └ @ SciMLBase ~/.julia/packages/SciMLBase/BoNUy/src/integrator_interface.jl:345
            # TODO(bernhard): we can force the solving by using force_dtmin = true
            u_δ18O_GWAT, u_δ2H_GWAT, u_δ18O_SWATI, u_δ2H_SWATI =
                compute_isotope_u_GWAT_SWATI(integrator,
                    # for GWAT:
                    u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT, du_GWFL[1], du_SEEP[1],
                    # for SWATI:
                    du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI,  # (non-fractionating)
                    aux_du_SLVP[1], p_fT_TADTM[1], p_VAPPRES(integrator.t), p_δ2H_PREC(integrator.t), p_δ18O_PREC(integrator.t), u_aux_WETNES, # (fractionating)
                    u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, 0, 0) #EffectiveDiffusivity_18O, EffectiveDiffusivity_2H)

            integrator.u.GWAT.d18O   = u_δ18O_GWAT
            integrator.u.GWAT.d2H    = u_δ2H_GWAT
            integrator.u.SWATI.d18O   = u_δ18O_SWATI
            integrator.u.SWATI.d2H    = u_δ2H_SWATI
        else
            error("Unknown method for updating Isotopes in GWAT and SWAT specified.")
        end
    end

    return nothing
end


function LWFBrook90R_updateIsotopes_GWAT_SWAT_AdvecDiff!(u, t, integrator)
    # This callback function gets called after an integration time step:
    #   +______+
    #   |      |     the solve() function uses the f function to integrate uᵏ to uᵏ⁺¹
    #   tᵏ     tᵏ⁺¹
    #   uᵏ     uᵏ⁺¹
    #               i.e. states where f returns a du ≂̸ 0 are already updated uᵏ⁺¹ ≂̸ uᵏ
    #               for  states where f returns a du = 0 are still uᵏ⁺¹ = uᵏ
    #               in the callback we compute uᵏ⁺¹
    #
    #               The order of integration is:
    #                   a) in each time step f and then cb_SWAT_GWAT_deltas
    #                   b) whenever a full day is over: finish f, finish cb_SWAT_GWAT_deltas, then the daily callbacks

    @unpack simulate_isotopes = integrator.p
    if simulate_isotopes
        ##### This update could be done in a separate FunctionCallingCallback `update_auxiliaries`
        # Unpack pre-allocated caches and update them with the current u
        # Bind memory to variable names to avoid re-allocating.
        # Note that these values will be overwritten, this line is just about memory allocation.
        # (u_aux_WETNES,u_aux_PSIM,u_aux_PSITI,u_aux_θ,u_aux_θ_tminus1,p_fu_KK,
        #     aux_du_DSFLI,aux_du_VRFLI,aux_du_VRFLI_1st_approx,aux_du_INFLI,aux_du_BYFLI, du_NTFLI,
        #     p_fu_BYFRAC) = integrator.p[4][1]
        # @unpack(u_aux_WETNES,u_aux_PSIM,u_aux_PSITI,u_aux_θ,u_aux_θ_tminus1,p_fu_KK,
        #     aux_du_DSFLI,aux_du_VRFLI,aux_du_VRFLI_1st_approx,aux_du_INFLI,aux_du_BYFLI, du_NTFLI,
        #     p_fu_BYFRAC) = integrator.p
        @unpack aux_du_DSFLI,aux_du_VRFLI, aux_du_INFLI, du_NTFLI = integrator.p
        #TODO(bernhard): check that u_aux_θ_tminus1 and u_aux_θ are indeed different
        ##### END This update could be done in a separate FunctionCallingCallback


        # # load pre-allocated memory locations (values are unimportant and will be overwritten further down)
        # (θᵏ⁺¹, θᵏ, C_¹⁸Oᵏ⁺¹, C_¹⁸Oᵏ, C_²Hᵏ⁺¹, C_²Hᵏ, q, Tsoil_K, τw, Λ,
        # D⁰_¹⁸O, D⁰_²H, D_¹⁸O_ᵏ⁺¹, D_²H_ᵏ⁺¹,
        # C_¹⁸O_SLVP, C_²H_SLVP,
        # # diff¹⁸O_interfaces, diff²H_interfaces, qCᵢ¹⁸O_interfaces, qCᵢ²H_interfaces,
        # diff¹⁸O_upp, diff²H_upp, qCᵢ¹⁸O_upp, qCᵢ²H_upp,
        # diff¹⁸O_low, diff²H_low, qCᵢ¹⁸O_low, qCᵢ²H_low,
        # , ) =  integrator.p[4][2]
        # TODO(bernhard): above generates somehow still many allocations
        @unpack cache_for_ADE_28 = integrator.p
        θᵏ⁺¹,θᵏ,C_¹⁸Oᵏ⁺¹,C_¹⁸Oᵏ,C_²Hᵏ⁺¹,C_²Hᵏ,q,
        D⁰_¹⁸O,D⁰_²H,D_¹⁸O_ᵏ⁺¹,D_²H_ᵏ⁺¹,
        C_¹⁸O_SLVP,C_²H_SLVP,
        diff¹⁸O_upp,diff²H_upp,qCᵢ¹⁸O_upp,qCᵢ²H_upp,
        diff¹⁸O_low,diff²H_low,qCᵢ¹⁸O_low,qCᵢ²H_low,
        du_Cᵢ¹⁸_SWATI,du_Cᵢ²H_SWATI,du_δ18O_SWATI,du_δ2H_SWATI,
        Tsoil_K,τw,p_Λ = cache_for_ADE_28 #integrator.p[4][2]
        # diff¹⁸O_interfaces,diff²H_interfaces,qCᵢ¹⁸O_interfaces,qCᵢ²H_interfaces =
        #     integrator.p[4][3]

        @unpack p_soil = integrator.p
        # @unpack p_STONEF, p_THICK = integrator.p.p_soil;
        p_STONEF = p_soil.p_STONEF
        p_THICK  = p_soil.p_THICK

        # FLOW state variables (i.e. amounts, already updated from tᵏ to tᵏ⁺¹)
        u_SWATIᵏ⁺¹ = integrator.u.SWATI.mm # @view integrator.u[    integrator.p[1][4].row_idx_SWATI, 1]
        u_SWATIᵏ   = integrator.uprev.SWATI.mm # @view integrator.uprev[integrator.p[1][4].row_idx_SWATI, 1]                 # of time step before
        u_GWATᵏ⁺¹  = integrator.u.GWAT.mm #integrator.u[1]

        # θᵏ       = u_aux_θ_tminus1 # # TODO(bernhard): u_aux_θ_tminus1 is not saved, workaraound below:
        _,_,_,θᵏ⁺¹[:],_ = LWFBrook90.KPT.derive_auxiliary_SOILVAR(u_SWATIᵏ⁺¹, p_soil)
        _,_,_,θᵏ[:],_   = LWFBrook90.KPT.derive_auxiliary_SOILVAR(u_SWATIᵏ, p_soil)  # of time step before

        # TRANSPORT state variables (i.e. concentrations, not yet updated from tᵏ to tᵏ⁺¹)
        u_δ18O_GWAT   = integrator.u.GWAT.d18O
        u_δ2H_GWAT    = integrator.u.GWAT.d2H
        # u_δ18O_INTS = integrator.u.INTS.d18O
        # u_δ2H_INTS  = integrator.u.INTS.d2H
        # u_δ18O_INTR = integrator.u.INTR.d18O
        # u_δ2H_INTR  = integrator.u.INTR.d2H
        # u_δ18O_SNOW = integrator.u.SNOW.d18O
        # u_δ2H_SNOW  = integrator.u.SNOW.d2H
        u_δ18O_SWATI  = integrator.u.SWATI.d18O # @view
        u_δ2H_SWATI   = integrator.u.SWATI.d2H  # @view
        u_δ18O_XYLEM  = integrator.u.XYLEM.d18O # @view
        u_δ2H_XYLEM   = integrator.u.XYLEM.d2H  # @view

        # unpack other needed quantities
        ## A) constant parameters:
        @unpack NLAYER, p_VXYLEM, p_DISPERSIVITY = integrator.p
        # TODO(bernhard): Note to below:
        #                 Fluxes in integrator.p[3][:] might be just the values overwritten by the last call to f.
        #                 For a multi-stage method (Runge-Kutta) this might not be representative for the entire time step Δt
        #                 Nevertheless it is a good approximation. Gold standard would be to cumulate the fluxes over time step Δt
        #                 and save into auxiliary states, but this would grow the state vector by at least 5*NLAYER
        @unpack aux_du_TRANI = integrator.p         # mm/day
        @unpack aux_du_DSFLI, aux_du_VRFLI, aux_du_INFLI, du_NTFLI = integrator.p # all in mm/day
        # @unpack u_aux_WETNES = integrator.p         # all in mm/day
        @unpack du_GWFL, du_SEEP = integrator.p     # all in mm/day

        ## B) time dependent parameters
        # @unpack p_DOY, p_MONTHN, p_GLOBRAD, p_TMAX, p_TMIN, p_VAPPRES, p_WIND, p_PREC,
        #     p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_fT_RELDEN,
        #     p_δ18O_PREC, p_δ2H_PREC, REFERENCE_DATE = integrator.p
        @unpack p_δ18O_PREC, p_δ2H_PREC = integrator.p

        ## C) state dependent parameters or intermediate results:
        # These were computed in the callback and are kept constant in between two
        # callbacks.
        @unpack p_fu_δ18O_SLFL, p_fu_δ2H_SLFL, p_fT_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP,
            p_fu_STHR, aux_du_RSNO, aux_du_SNVP,
            aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, u_SNOW_old = integrator.p

        # Define quantities needed for transport equation (advection dispersion/diffusion equation)
        δ¹⁸Oᵏ   = integrator.uprev.SWATI.d18O # @view integrator.uprev[integrator.p[1][4].row_idx_SWATI, integrator.p[1][4].col_idx_d18O]  # of time step before
        δ²Hᵏ    = integrator.uprev.SWATI.d2H # @view integrator.uprev[integrator.p[1][4].row_idx_SWATI, integrator.p[1][4].col_idx_d2H]   # of time step before
        # ## unused: Note: transform the deltas into isotope-amount fraction x, a.k.a. isotopic abundance, a.k.a. atom fraction
        # Convert to concentrations:
        C_¹⁸Oᵏ   .= LWFBrook90.ISO.δ_to_C(δ¹⁸Oᵏ, LWFBrook90.ISO.R_VSMOW¹⁸O, LWFBrook90.ISO.Mi_¹⁸O) # kg/m3
        C_²Hᵏ    .= LWFBrook90.ISO.δ_to_C(δ²Hᵏ,  LWFBrook90.ISO.R_VSMOW²H,  LWFBrook90.ISO.Mi_²H)  # kg/m3

        q .= aux_du_VRFLI ./ 1000 # VRFLI is in mm/day, q is now in m/day
        # Compute quantities needed for advection diffusion equation
        # # For each Finite Volume cell:
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
            #
            #     from Braud et al. 2005 use the following approximation at cell boundaries
            #       a) for conc gradient: dC/dz(z_upper) = dC/dz(j-1/2) = (C(j) - C(j-1)) / dz(j-1) # finite differences
            #       b) for concentration: C(z_upper)     = C(j-1/2)     = (C(j) + C(j-1)) /2        # no specific upwind treatment
            #       c) for dispersivity:  D(z_upper)     = D(j-1/2)     = (D(j) + D(j-1)) /2        # no specific upwind treatment

            # ##################
            # ##################
            # ##################
            # 0) Define conditions for isotope calculations:
            Tc = p_fT_TADTM[1]  # °C, average daytime air temperature
            # h = min(1.0, p_fT_VAPPRES / LWFBrook90.PET.ESAT(Tc)[1]) # -, relative humidity of the atmosphere (vappress_atm/1 atm)
            # γ = 1.0          # -, thermodynamic activity coefficient of evaporating water
            # # X_INTS = 0.5  # -, turbulence incex of the atmosphere above the evaporating water
            # # X_INTR = 0.5  # -, turbulence incex of the atmosphere above the evaporating water
            # # X_SNOW = 1.0  # -, turbulence incex of the atmosphere above the evaporating water

            # # Atmospheric vapor composition assumed to be in equilibrium with precipitation
            # #     Benettin-2018-Hydrol_Earh_Syst_Sci citing Gibson et al. 2008
            # δₐ(δ_p, α_eq) = (δ_p - 1000 * (α_eq - 1))/α_eq
            # δ²H_a  = δₐ(p_δ2H_PREC,  LWFBrook90.ISO.α²H_eq(Tc))
            # δ¹⁸O_a = δₐ(p_δ18O_PREC, LWFBrook90.ISO.α¹⁸O_eq(Tc))

            # # For soil:
            # Xa = 0.5 # between molecular and turbulent
            # Xs = 1.0 # molecular only
            # X_SOIL = u_aux_WETNES[1] * Xa + (1 - u_aux_WETNES[1]) * Xs
            # # when fully saturated, WETNES == 1, water vapor leaves quickly to atmosphere
            # # when hardly saturated, WETNES -> 0, water vapor crosses soil pores until it reaches atmosphere
            # # Equivalent:
            # # X_SOIL = ((u_aux_θ[1]-θ_res)*Xa + (θ_sat-u_aux_θ[1])*Xs)/(θ_sat-θ_res)
            # # Taken from Zhou-2021-Environ_Model_Softw (X is called n_k there)
            # ##################
            # ##################
            # ##################


        # Define (constant) soil transport properties
        N = NLAYER
        Tsoil_K .= (p_fT_TADTM[1] + 273.15) .* ones(N) # °C, TODO: use solution from heat equation instead of approximation of p_fT_TADTM
        # τw = θᵏ⁺¹ .^ (7/3) ./ (θsat .^ 2)         # TODO: express tortuosity as function of θ, (Millington and Quirk 1961 as shown in Radcliffe et al. 2018, eq 6.6)
        τw .= 1.0 .* ones(N)   # -, tortuosity in liquid phase (w = water), using 1.0 will overestimate diffusion
        # τg = 1.0 .* ones(N) # -, tortuosity in vapor phase (g = gas), unused as no vapor transport is considered
        # p_Λ = p_soil.p_DISPER #(instead of p_DISPERSIVITY .* ones(N))         # m, dispersivity length (if we want to have different dispersivities per soil layer)
        D⁰_¹⁸O .= 0.96691 .* 10^-9 .* exp.(-535400 ./ Tsoil_K.^2 .+ 1393.3 ./ Tsoil_K .+ 2.1876)  .* 3600 .* 24 # m²/day molecular diffusion constant of ¹⁸O in liquid water (eq. A3, Zhou et al. 2021)
        D⁰_²H  .= 0.98331 .* 10^-9 .* exp.(-535400 ./ Tsoil_K.^2 .+ 1393.3 ./ Tsoil_K .+ 2.1876)  .* 3600 .* 24 # m²/day molecular diffusion constant of ²H  in liquid water (eq. A3, Zhou et al. 2021)

        # Effective dispersion coefficient (Eq. 33 Zhou et al. 2021)
        D_¹⁸O_ᵏ⁺¹ .= D⁰_¹⁸O .* τw .* θᵏ⁺¹ .+ p_DISPERSIVITY .* abs.(q) # m²/day
        D_²H_ᵏ⁺¹  .= D⁰_²H  .* τw .* θᵏ⁺¹ .+ p_DISPERSIVITY .* abs.(q) # m²/day

        # Define concentrations of source/sink terms in transport equation (TRANI, DSFLI, INFLI, SLVP)
        ### Define δ signature of in- and outflows
        # TODO(bernhard) for debugging:
        δ18O_INFLI = p_δ18O_PREC(integrator.t)#TODO(bernhard): debug remove workaround and set again = p_fu_δ18O_SLFL
        δ2H_INFLI  = p_δ2H_PREC(integrator.t) #TODO(bernhard): debug remove workaround and set again = p_fu_δ2H_SLFL

        C_¹⁸O_INFLI = LWFBrook90.ISO.δ_to_C.(δ18O_INFLI, LWFBrook90.ISO.R_VSMOW¹⁸O, LWFBrook90.ISO.Mi_¹⁸O) # for debugging use: LWFBrook90.ISO.δ_to_C.(p_δ18O_PREC(integrator.tprev), LWFBrook90.ISO.R_VSMOW¹⁸O, LWFBrook90.ISO.Mi_¹⁸O)
        # C_¹⁸O_TRANI = C_¹⁸Oᵏ # no fractionation occurring, i.e. outflux composition equal to storage composition
        # C_¹⁸O_DSFL  = C_¹⁸Oᵏ # no fractionation occurring, i.e. outflux composition equal to storage composition
        C_²H_INFLI  = LWFBrook90.ISO.δ_to_C.(δ2H_INFLI,  LWFBrook90.ISO.R_VSMOW²H,  LWFBrook90.ISO.Mi_²H) # for debugging use: LWFBrook90.ISO.δ_to_C.(p_δ2H_PREC(integrator.tprev), LWFBrook90.ISO.R_VSMOW²H,  LWFBrook90.ISO.Mi_²H)
        # C_²H_TRANI  = C_²Hᵏ # no fractionation occurring, i.e. outflux composition equal to storage composition
        # C_²H_DSFL   = C_²Hᵏ # no fractionation occurring, i.e. outflux composition equal to storage composition
        # TODO(bernhard): for TRANI and DSFL we are instead using Cᵏ⁺¹ (i.e. use the concentrations in the LHS in the implicit scheme...)

        ### Compute δ signature of evaporating flux
        δ¹⁸O_SLVP = u_δ18O_SWATI[1] # TODO(bernhard): disabling evaporation fractionation (TODO: should yield solution similar to Stumpp 2012 model)
        δ²H_SLVP  = u_δ2H_SWATI[1]  # TODO(bernhard): disabling evaporation fractionation (TODO: should yield solution similar to Stumpp 2012 model)
        # Equation derived based on Gonfiantini (see 60 lines above in comment)
        # # δ_E = 1000*(1 + 1/((γ - h)*(α_dif)^X) * (γ/α*(1+δ_w/1000) - h*(1+δ_A/1000)))
        # δ¹⁸O_SLVP = 1000*( 1 + (γ/α¹⁸O_eq*(1 + u_δ18O_SWATI[1] / 1000) - h*(1 + δ¹⁸O_a/1000)) /
        #                         ((γ - h)*(LWFBrook90.ISO.α¹⁸O_dif)^X_SOIL) )
        # δ²H_SLVP  = 1000*( 1 + (γ/α²H_eq*(1 + u_δ2H_SWATI[1] / 1000) - h*(1 + δ²H_a/1000)) /
        #                         ((γ - h)*(LWFBrook90.ISO.α²H_dif)^X_SOIL) )
        # (Above is an alternative to formulation in Benettin 2018 HESS eq. 1 and Gibson 2016)
        # Cᵢ_SLVP = ( (Cᵢ - ε¹⁸O_eq)/α¹⁸O_eq - h*δ¹⁸O_a - ε¹⁸O_dif ) /
        #             (1 - h + ε¹⁸O_dif/1000) # [‰]

        C_¹⁸O_SLVP .= 0
        C_²H_SLVP  .= 0
        C_¹⁸O_SLVP[1]  = LWFBrook90.ISO.δ_to_C(δ¹⁸O_SLVP, LWFBrook90.ISO.R_VSMOW¹⁸O, LWFBrook90.ISO.Mi_¹⁸O)
        C_²H_SLVP[1]   = LWFBrook90.ISO.δ_to_C(δ²H_SLVP,  LWFBrook90.ISO.R_VSMOW²H,  LWFBrook90.ISO.Mi_²H)
        # E¹⁸O = C_¹⁸O_SLVP * aux_du_SLVP[1] * 0.001 # kg/m3 * mm/day * 0.001 m/mm # (kg/m²/day¹)
        # E²H  = C_²H_SLVP  * aux_du_SLVP[1] * 0.001 # kg/m3 * mm/day * 0.001 m/mm # (kg/m²/day¹)

        ### Prepare terms to evaluate linear system to be solved
        Δt = integrator.t - integrator.tprev # days
        # TODO(bernhard): below dz and Δz can be computed once only when setting up the simulation
        Δz = integrator.p.p_soil.p_THICK / 1000 # m
        # z_center = (cumsum(Δz) .+ cumsum([0; Δz[1:NLAYER-1]]))/2 # m
        dz = diff((cumsum(Δz) .+ cumsum([0; Δz[1:NLAYER-1]]))/2)   # m

        ################
        ################
        ################
        # SOLVING VARIANT A (Braud et al. 2005, implicit time stepping) -> has a bug as it seems unstable (<- might it be linked to not using an upwind scheme?)
            # # #TODO(bernhard): for debugging of my adaptation of Braud et al. 2005 try
            # #                  to make a mistake and see if the solution becomes more
            # #                  stable if we do (mistake on purpose!)
            # # Variant 1:
            # θᵏ = θᵏ⁺¹
            # Δt = integrator.t - integrator.tprev # days
            # # θᵏ⁺¹ = (u_SWATIᵏ⁺¹ + (du_NTFLI)*integrator.dt) ./ p_THICK ./ (1.0 .- p_STONEF)
            # θᵏ⁺¹ = (max.(0.001, u_SWATIᵏ⁺¹) + (du_NTFLI)*Δt) ./ p_THICK ./ (1.0 .- p_STONEF)
            # # Variant 2:
            # θᵏ⁺¹ = u_SWATIᵏ⁺¹ ./ p_THICK ./ (1.0 .- p_STONEF)
            # Δt = integrator.t - integrator.tprev # days
            # θᵏ = (u_SWATIᵏ⁺¹ -  (du_NTFLI)*Δt)./ p_THICK ./ (1.0 .- p_STONEF)
            # # #TODO(bernhard): end debugging

        # # Setup linear system (implicit solver) (From Braud et al. 2005)

        # Define boundary conditions
        ### TOP: surface isotopic flux E = 0 (as evaporation in LWFBrook is implemented as sink term (see above E¹⁸O))
        BCFlux_top¹⁸O    = 0 # kg/m²/day¹
        BCFlux_top²H     = 0 # kg/m²/day¹
        ### BOTTOM:
        BCFlux_bottom¹⁸O = q[N] * C_¹⁸Oᵏ[N]  # kg/m²/day¹
        BCFlux_bottom²H  = q[N] * C_²Hᵏ[N]   # kg/m²/day¹

        # # Make matrix for implicit time stepping (solving Ax = b) (see Braud et al. 2005)
        # # TODO(bernharf): preallocate this matrix A (and vectors...)
        # A1_¹⁸O = ones(N)
        # A2_¹⁸O = ones(N)
        # A3_¹⁸O = ones(N)
        # b_¹⁸O  = ones(N)
        # A1_²H = ones(N)
        # A2_²H = ones(N)
        # A3_²H = ones(N)
        # b_²H  = ones(N)
        # # to concatenate to a Tridiagonal matrix as follows: A = Array(Tridiagonal(A1[2:end],A2,A3[1:end-1]))

        #     # FOR DEBUGGING: other working solution basically does: u⁺ = u₀ + dt * (sum(inflow) - sum(outflow)) # [mm]
        #     #                                                       x⁺ = (x₀ * u₀ + dt * (sum(xin .* inflow) - sum(x₀ .* outflow)) ) / u⁺
        #     #      where u is u_SWATI
        #     #      from derive_auxiliary_SOILVAR() we have
        #     #           θ = θr + S * (θs - θr)
        #     #             = θr +  min.(1, S) * (θs - θr)
        #     #             = θr +  min.(1, (θs * SWATI / SWATMAX - θr)/(θs - θr) ) * (θs - θr)
        #     #             = θr +  min.((θs - θr), (θs * SWATI / SWATMAX - θr) )
        #     #             = min.(θs, θs * SWATI / SWATMAX )
        #     #             = min.(θs, SWATI / (p_THICK .* (1.0 .- p_STONEF)) )
        #     #   is maybe the fact that we have a min() in there that change in θ and q are not fully compatible...?

        # ### Inner nodes:
        # for j = 2:(N-1)
        #     # Linear system from Braud et al. 2005:
        #     # A1[j] = - Diⱼ₋₁₎₂ / dzⱼ₋₁                    -   qⱼ₋₁₎₂ / 2
        #     # A2[j] = - Δzⱼ * θᵏ⁺¹ / Δt    +  qⱼ₊₁₎₂ / 2   -   qⱼ₋₁₎₂ / 2   +    Diⱼ₊₁₎₂ / dzⱼ   +   Diⱼ₋₁₎₂ / dzⱼ₋₁
        #     # A3[j] = - Diⱼ₊₁₎₂ / dzⱼ     +   qⱼ₋₁₎₂ / 2
        #     # b[j] = Δzⱼ * θᵏ * Ciⱼᵏ / Δt # Note: θᵏ instead of θᵏ⁺¹ (Braud eq. D.15)

        #     # Library of available expressions:
        #     # D_¹⁸O_ᵏ⁺¹[j-1]
        #     # D_¹⁸O_ᵏ⁺¹[j]
        #     # D_¹⁸O_ᵏ⁺¹[j+1]
        #     # q[j-1]
        #     # q[j]
        #     # q[j+1]
        #     # Δt
        #     # Δz[j-1], Δz[j], Δz[j+1]
        #     # dz[j-1], dz[j], dz[j+1]

        #     # A1_¹⁸O[j] =                                                                   - (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1]    -   q[j-1] / 2
        #     # A3_¹⁸O[j] =                    - (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j]                                                                     +   q[j] / 2
        #     # A2_¹⁸O[j] = Δz[j]/Δt * θᵏ⁺¹[j] + (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j]  + (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1]    -   q[j-1] / 2  +   q[j] / 2 + (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000) # Factor 1000 to get from mm/day to m/day

        #     # # RHS (using C_¹⁸Oᵏ instead of C_¹⁸Oᵏ⁺¹):
        #     # b_¹⁸O[j]  = Δz[j]/Δt * θᵏ[j]   * C_¹⁸Oᵏ[j]  + aux_du_INFLI[j]/1000 * C_¹⁸O_INFLI # - aux_du_SLVP[j]/1000 * C_¹⁸O_SLVP[j] # NOTE: SLVP only applies to cell 1

        #     # # explicit solver:
        #     # Variant 1: interface extrapolation as by Braud 2005, but just using explicit time stepping C_¹⁸Oᵏ
        #     # A1_¹⁸O[j] =  0
        #     # A3_¹⁸O[j] =  0
        #     # A2_¹⁸O[j] = Δz[j]/Δt * θᵏ⁺¹[j]   # + (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000) # Factor 1000 to get from mm/day to m/day
        #     # b_¹⁸O[j]  = aux_du_INFLI[j]/1000 * C_¹⁸O_INFLI +
        #     #             C_¹⁸Oᵏ[j-1] * (                                                                                 (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1] + q[j-1] / 2 ) +
        #     #             C_¹⁸Oᵏ[j] *   ( Δz[j]/Δt * θᵏ[j] -   (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j] -   q[j] / 2 - (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1] + q[j-1] / 2 - (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000)) +
        #     #             C_¹⁸Oᵏ[j+1] * (                      (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j] -   q[j] / 2 )
        #     # Variant 2: interface extrapolation based on upper cell (kind of an upwind scheme assuming flow is always downward...)
        #     A1_¹⁸O[j] =  0
        #     A3_¹⁸O[j] =  0
        #     A2_¹⁸O[j] = Δz[j]/Δt * θᵏ⁺¹[j]   # + (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000) # Factor 1000 to get from mm/day to m/day
        #     b_¹⁸O[j]  = aux_du_INFLI[j]/1000 * C_¹⁸O_INFLI +
        #                 C_¹⁸Oᵏ[j-1] * (                                                                                 (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1] + q[j-1] ) +
        #                 C_¹⁸Oᵏ[j] *   ( Δz[j]/Δt * θᵏ[j] -   (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j] -   q[j]     - (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1] + q[j-1]*0 - (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000)) +
        #                 C_¹⁸Oᵏ[j+1] * (                      (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j] -   q[j]*0 )

        # end
        # ### Outer nodes (including boundary conditions):
        # for j = 1
        #     # # # using a Dirchlet boundary condition
        #     # # A1[j] = 0
        #     # # A2[j] = 1
        #     # # A3[j] = 0
        #     # # b[j] = Ciₜₒₚ
        #     # # using a Neumann boundary condition (flux = E)
        #     # # A1[j] = 0
        #     # # A2[j] = - Δzⱼ * θᵏ / Δt     +   qⱼ₊₁₎₂ / 2                    +    Diⱼ₊₁₎₂ / dzⱼ
        #     # # A3[j] = - Diⱼ₊₁₎₂ / dzⱼ     +   qⱼ₋₁₎₂ / 2
        #     # # b[j]  = Δzⱼ * θᵏ / Δt * Cₛ¹⁸O + BCFlux_top¹⁸O
        #     # A1_¹⁸O[j] = 0
        #     # A3_¹⁸O[j] =                    - (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j]                                                                     +   q[j] / 2
        #     # A2_¹⁸O[j] = Δz[j]/Δt * θᵏ⁺¹[j] + (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j]                                                                     +   q[j] / 2 + (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000) # Factor 1000 to get from mm/day to m/day
        #     # b_¹⁸O[j]  = Δz[j]/Δt * θᵏ[j]   * C_¹⁸Oᵏ[j]  + aux_du_INFLI[j]/1000 * C_¹⁸O_INFLI - aux_du_SLVP[j]/1000 * C_¹⁸O_SLVP[j] + BCFlux_top¹⁸O

        #     # # explicit solver:
        #     # Variant 1: interface extrapolation as by Braud 2005, but just using explicit time stepping C_¹⁸Oᵏ
        #     # A1_¹⁸O[j] =  0
        #     # A3_¹⁸O[j] =  0
        #     # A2_¹⁸O[j] = Δz[j]/Δt * θᵏ⁺¹[j]   # + (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000) # Factor 1000 to get from mm/day to m/day
        #     # b_¹⁸O[j]  = aux_du_INFLI[j]/1000 * C_¹⁸O_INFLI +
        #     #             #C_¹⁸Oᵏ[j-1] * ( 0                                                                             ) +
        #     #             C_¹⁸Oᵏ[j] *   ( Δz[j]/Δt * θᵏ[j] -   (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j] -   q[j] / 2 - (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000)) +
        #     #             C_¹⁸Oᵏ[j+1] * (                      (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j] -   q[j] / 2 )
        #     # Variant 2: interface extrapolation based on upper cell (kind of an upwind scheme assuming flow is always downward...)
        #     A1_¹⁸O[j] =  0
        #     A3_¹⁸O[j] =  0
        #     A2_¹⁸O[j] = Δz[j]/Δt * θᵏ⁺¹[j]   # + (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000) # Factor 1000 to get from mm/day to m/day
        #     b_¹⁸O[j]  = aux_du_INFLI[j]/1000 * C_¹⁸O_INFLI +
        #                 #C_¹⁸Oᵏ[j-1] * ( 0                                                                             ) +
        #                 C_¹⁸Oᵏ[j] *   ( Δz[j]/Δt * θᵏ[j] -   (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j] -   q[j]   - (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000)) +
        #                 C_¹⁸Oᵏ[j+1] * (                      (D_¹⁸O_ᵏ⁺¹[j+1] + D_¹⁸O_ᵏ⁺¹[j]) / 2 / dz[j] -   q[j]*0 )
        # end
        # for j = N
        #     # # # using a Neuman boundary condition
        #     # # A1[j] = - Diⱼ₋₁₎₂ / dzⱼ       -   qⱼ₋₁₎₂ / 2
        #     # # A2[j] = - 2 * Δzⱼ * θᵏ / Δt   -   qⱼ₋₁₎₂ / 2   +    Diⱼ₋₁₎₂ / dzⱼ₋₁ # TODO(bernhard): where is the factor 2 coming from?
        #     # # A3[j] = - 0
        #     # # b[j] = 2*Δzⱼ * θᵏ * Ciⱼᵏ / Δt                                       # TODO(bernhard): where is the factor 2 coming from?

        #     # A1_¹⁸O[j] =                                                                   - (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1]    -   q[j-1] / 2
        #     # A3_¹⁸O[j] = 0
        #     # # TODO(bernhard): in Braud2005 there are factors 2 for A2 and b next to Δz[j] -> 2*Δz[j], is that correct or not?
        #     # # A2_¹⁸O[j] = 2*Δz[j]/Δt * θᵏ⁺¹[j]                                              + (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1]    -   q[j-1] / 2              + (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000) # Factor 1000 to get from mm/day to m/day
        #     # # b_¹⁸O[j]  = 2*Δz[j]/Δt * θᵏ[j]   * C_¹⁸Oᵏ[j]  -  BCFlux_bottom¹⁸O
        #     # A2_¹⁸O[j] = Δz[j]/Δt * θᵏ⁺¹[j]              +  q[N]               + (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1]    -   q[j-1] / 2              + (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000) # Factor 1000 to get from mm/day to m/day
        #     # b_¹⁸O[j]  = Δz[j]/Δt * θᵏ[j]   * C_¹⁸Oᵏ[j]
        #     # # A2_¹⁸O[j] = Δz[j]/Δt * θᵏ[j]              +  q[N]               + (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1]    -   q[j-1] / 2              + (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000) # Factor 1000 to get from mm/day to m/day
        #     # # b_¹⁸O[j]  = Δz[j]/Δt * (θᵏ[j] + Δt * du_NTFLI[j] / (Δz[j]*(1-p_STONEF[j])))   * C_¹⁸Oᵏ[j]  + aux_du_INFLI[j]/1000 * C_¹⁸O_INFLI

        #     # # explicit solver:
        #     # Variant 1: interface extrapolation as by Braud 2005, but just using explicit time stepping C_¹⁸Oᵏ
        #     # A1_¹⁸O[j] =  0
        #     # A3_¹⁸O[j] =  0
        #     # A2_¹⁸O[j] = Δz[j]/Δt * θᵏ⁺¹[j]   # + (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000) # Factor 1000 to get from mm/day to m/day
        #     # b_¹⁸O[j]  = aux_du_INFLI[j]/1000 * C_¹⁸O_INFLI +
        #     #             C_¹⁸Oᵏ[j-1] * (                    (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1] + q[j-1] / 2 ) +
        #     #             C_¹⁸Oᵏ[j] *   ( Δz[j]/Δt * θᵏ[j] - (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1] + q[j-1] / 2   -  q[N] - (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000)) #+
        #     #             #C_¹⁸Oᵏ[j+1] * ( 0                )
        #     # Variant 2: interface extrapolation based on upper cell (kind of an upwind scheme assuming flow is always downward...)
        #     A1_¹⁸O[j] =  0
        #     A3_¹⁸O[j] =  0
        #     A2_¹⁸O[j] = Δz[j]/Δt * θᵏ⁺¹[j]   # + (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000) # Factor 1000 to get from mm/day to m/day
        #     b_¹⁸O[j]  = aux_du_INFLI[j]/1000 * C_¹⁸O_INFLI +
        #                 C_¹⁸Oᵏ[j-1] * (                    (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1] + q[j-1] ) +
        #                 C_¹⁸Oᵏ[j] *   ( Δz[j]/Δt * θᵏ[j] - (D_¹⁸O_ᵏ⁺¹[j] + D_¹⁸O_ᵏ⁺¹[j-1]) / 2 / dz[j-1] + q[j-1]*0 -  q[N] - (aux_du_TRANI[j]/1000 + aux_du_DSFLI[j]/1000)) #+
        #                 #C_¹⁸Oᵏ[j+1] * ( 0                )

        # end
        # A_¹⁸O = Array(Tridiagonal(A1_¹⁸O[2:end],A2_¹⁸O,A3_¹⁸O[1:end-1]))

        # # Compute the updated concentrations:
        # C_¹⁸Oᵏ⁺¹ .= A_¹⁸O \ b_¹⁸O
        # u_δ18O_SWATI = LWFBrook90.ISO.C_to_δ(C_¹⁸Oᵏ⁺¹, LWFBrook90.ISO.R_VSMOW¹⁸O, LWFBrook90.ISO.Mi_¹⁸O) # in permil, transform from C to δ
        # # C_²Hᵏ⁺¹ = A_²H \ b_²H
        # # u_δ2H_SWATI  = LWFBrook90.ISO.C_to_δ(C_²Hᵏ⁺¹ , LWFBrook90.ISO.R_VSMOW²H , LWFBrook90.ISO.Mi_²H ) # in permil, transform from C to δ

        # # Compute balance for GWAT (for d2H and for d18O)
        # u[idx_u_vector_isotopes_d18O]      = u_δ18O_SWATI
        # # u[idx_u_vector_isotopes_d2H]       = u_δ2H_SWATI

        # # Update GWAT:
        # # TODO(bernhard): implement this
        # integrator.u.GWATd18O   = u_δ18O_GWAT
        # integrator.u.GWATd2H   = u_δ2H_GWAT

        # END SOLVING VARIANT A (Braud et al. 2005, implicit time stepping) -> has a bug as it seems unstable
        ################
        ################
        ################


        ################
        ################
        ################
        # SOLVING VARIANT B explicit time stepping
        # 1) Compute isotopic composition (or solute concentrations) of evaporating soil water
        ### a) Fix constants
        # α¹⁸O_eq = LWFBrook90.ISO.α¹⁸O_eq(Tc)
        # α²H_eq  = LWFBrook90.ISO.α²H_eq(Tc)
        # ε¹⁸O_dif = (α¹⁸O_dif-1)*1000  # diffusive (i.e. kinetic fractionation)
        # ε²H_dif  = (α²H_dif-1)*1000   # diffusive (i.e. kinetic fractionation)
        # ε¹⁸O_eq = (α¹⁸O_eq-1)*1000    # equilibrium fractionation
        # ε²H_eq  = (α²H_eq-1)*1000     # equilibrium fractionation

        ### b1) Compute δ signature of evaporating flux (for numerical solution)
        # Equation derived based on Gonfiantini (see 60 lines above in comment)
        # δ_E = 1000*(1 + 1/((γ - h)*(α_dif)^X) * (γ/α*(1+δ_w/1000) - h*(1+δ_A/1000)))
        # δ¹⁸O_SLVP = 1000*( 1 + (γ/α¹⁸O_eq*(1 + u_δ18O_SWATI[1] / 1000) - h*(1 + δ¹⁸O_a/1000)) /
        #                         ((γ - h)*(LWFBrook90.ISO.α¹⁸O_dif)^X_SOIL)
        # )
        # δ²H_SLVP  = 1000*( 1 + (γ/α²H_eq*(1 + u_δ2H_SWATI[1] / 1000) - h*(1 + δ²H_a/1000)) /
        #                         ((γ - h)*(LWFBrook90.ISO.α²H_dif)^X_SOIL)
        # )
        # (Above is an alternative to formulation in Benettin 2018 HESS eq. 1 and Gibson 2016)
        # Cᵢ_SLVP = ( (Cᵢ - ε¹⁸O_eq)/α¹⁸O_eq - h*δ¹⁸O_a - ε¹⁸O_dif ) /
        #             (1 - h + ε¹⁸O_dif/1000) # [‰]

        # 2) Compute isotopic composition (or solute concentrations) of all other water fluxes
        # Note: transform the fluxes into isotope-amount fraction x, a.k.a. isotopic abundance, a.k.a. atom fraction
        C_¹⁸Oᵏ .= LWFBrook90.ISO.δ_to_x.(δ¹⁸Oᵏ, LWFBrook90.ISO.R_VSMOW¹⁸O) # -, isotope amounts
        C_²Hᵏ  .= LWFBrook90.ISO.δ_to_x.(δ²Hᵏ,  LWFBrook90.ISO.R_VSMOW²H)  # -, isotope amounts
        # C_¹⁸Oᵏ   = LWFBrook90.ISO.δ_to_C(δ¹⁸Oᵏ, LWFBrook90.ISO.R_VSMOW¹⁸O, LWFBrook90.ISO.Mi_¹⁸O) # kg/m3, # concentrations
        # C_²Hᵏ    = LWFBrook90.ISO.δ_to_C(δ²Hᵏ,  LWFBrook90.ISO.R_VSMOW²H,  LWFBrook90.ISO.Mi_²H)  # kg/m3, # concentrations
        C_¹⁸O_Xylem = LWFBrook90.ISO.δ_to_x(u_δ18O_XYLEM, LWFBrook90.ISO.R_VSMOW¹⁸O) # -, isotope amounts
        C_²H_Xylem  = LWFBrook90.ISO.δ_to_x(u_δ2H_XYLEM,  LWFBrook90.ISO.R_VSMOW²H)  # -, isotope amounts


        #                     D_at_interface                                    .*   dC/dz_at_interface
        diff¹⁸O_upp .= [0; (D_¹⁸O_ᵏ⁺¹[1:NLAYER-1] + D_¹⁸O_ᵏ⁺¹[2:NLAYER]) /2]     .* [0; [(C_¹⁸Oᵏ[i] - C_¹⁸Oᵏ[i-1]) / dz[i-1] for i ∈ 2:length(C_¹⁸Oᵏ)]]    # flux in m/d (D: m2/d, dx/dz: m-1)
        diff¹⁸O_low .= [   (D_¹⁸O_ᵏ⁺¹[1:NLAYER-1] + D_¹⁸O_ᵏ⁺¹[2:NLAYER]) /2; 0]  .* [   [(C_¹⁸Oᵏ[i] - C_¹⁸Oᵏ[i-1]) / dz[i-1] for i ∈ 2:length(C_¹⁸Oᵏ)]; 0] # flux in m/d (D: m2/d, dx/dz: m-1)
        diff²H_upp  .= [0; (D_²H_ᵏ⁺¹[1:NLAYER-1]  + D_²H_ᵏ⁺¹[2:NLAYER] ) /2]     .* [0; [(C_²Hᵏ[i]  - C_²Hᵏ[i-1] ) / dz[i-1] for i ∈ 2:length(C_²Hᵏ )]]    # flux in m/d (D: m2/d, dx/dz: m-1)
        diff²H_low  .= [   (D_²H_ᵏ⁺¹[1:NLAYER-1]  + D_²H_ᵏ⁺¹[2:NLAYER] ) /2; 0]  .* [   [(C_²Hᵏ[i]  - C_²Hᵏ[i-1] ) / dz[i-1] for i ∈ 2:length(C_²Hᵏ )]; 0] # flux in m/d (D: m2/d, dx/dz: m-1)

        # TODO(bernhard): below C's are not at interface but in the middle of the cell... best would be to use an upwind scheme and use the correct C based on the sign of aux_du_VRFLI
        #                     q_at_interface                                    .*   C_at_interface
        qCᵢ¹⁸O_upp  .= [0; aux_du_VRFLI[1:(NLAYER-1)]] .* [0; C_¹⁸Oᵏ[1:(NLAYER-1)]]
        qCᵢ¹⁸O_low  .=     aux_du_VRFLI[1:(NLAYER)]    .*     C_¹⁸Oᵏ[1:(NLAYER)]
        qCᵢ²H_upp   .= [0; aux_du_VRFLI[1:(NLAYER-1)]] .* [0; C_²Hᵏ[1:(NLAYER-1)]]
        qCᵢ²H_low   .=     aux_du_VRFLI[1:(NLAYER)]    .*     C_²Hᵏ[1:(NLAYER)]

        Cᵢ¹⁸O_INFLI = LWFBrook90.ISO.δ_to_x.(p_δ18O_PREC(integrator.t), LWFBrook90.ISO.R_VSMOW¹⁸O)  # TODO(bernhard): for debugging, remove this again and replace with δ18O_INFLI
        # Cᵢ¹⁸O_INFLI = δ18O_INFLI
        # Cᵢ¹⁸O_INFLI = ifelse(sum(aux_du_INFLI) == 0, 0, δ18O_INFLI) # in case there is no inflow δ18O_INFLI was set to NaN, set it to zero for below equation
        Cᵢ¹⁸O_TRANI = ifelse.(aux_du_TRANI .< 0,
                              C_¹⁸O_Xylem,
                              C_¹⁸Oᵏ)
                             # Under root water uptake assume no fractionation occurring, i.e. outflux composition equal to storage composition C_¹⁸Oᵏ
                             # however if aux_du_TRANI is negative for some layers it means water from the root xylem flows out into the soil domain
                             # This is hereby accounted for the the isotope balance in soil water, as well as for the Xylem and δ_RWU
        Cᵢ¹⁸O_DSFL  = C_¹⁸Oᵏ # no fractionation occurring, i.e. outflux composition equal to storage composition
        C_¹⁸O_SLVP .= 0
        C_¹⁸O_SLVP[1]  = LWFBrook90.ISO.δ_to_x.(δ¹⁸O_SLVP, LWFBrook90.ISO.R_VSMOW¹⁸O)
        Cᵢ²H_INFLI = LWFBrook90.ISO.δ_to_x.(p_δ2H_PREC(integrator.t),  LWFBrook90.ISO.R_VSMOW²H )   # TODO(bernhard): for debugging, remove this again and replace with δ2H_INFLI
        # Cᵢ²H_INFLI = δ2H_INFLI
        # Cᵢ²H_INFLI = ifelse(sum(aux_du_INFLI) == 0, 0, δ2H_INFLI) # in case there is no inflow δ2H_INFLI was set to NaN, set it to zero for below equation
        Cᵢ²H_TRANI = C_²Hᵏ # no fractionation occurring, i.e. outflux composition equal to storage composition
        Cᵢ²H_TRANI = ifelse.(aux_du_TRANI .< 0,
                              C_²H_Xylem,
                              C_²Hᵏ)
                             # Under root water uptake assume no fractionation occurring, i.e. outflux composition equal to storage composition C_²Hᵏ
                             # however if aux_du_TRANI is negative for some layers it means water from the root xylem flows out into the soil domain
                             # This is hereby accounted for the the isotope balance in soil water, as well as for the Xylem and δ_RWU
        Cᵢ²H_DSFL  = C_²Hᵏ # no fractionation occurring, i.e. outflux composition equal to storage composition
        C_²H_SLVP    .= 0
        C_²H_SLVP[1]  = LWFBrook90.ISO.δ_to_x.(δ²H_SLVP,  LWFBrook90.ISO.R_VSMOW²H )

        # 3) Some other terms in the isotope balance equation
        dVdt = du_NTFLI # [mm/day]
        # θ = u_aux_θ
        # dθdt = aux_du_VRFLI[NLAYER] - du_GWFL - du_SEEP # dV/dt [-/day] # TODO(bernhard)
                        # du_NTFL was computed in MSBITERATE as following:
                                    #                      = ∫∂/∂z[K_h(∂h/∂z + 1)] dz   - ∫S dz
                                    # for i=1:    NTFLI[i] = 0          - VRFLI[i]      + INFLI[i] - aux_du_TRANI[i] - aux_du_DSFLI[i] - aux_du_SLVP
                                    # for i=x:    NTFLI[i] = VRFLI[i-1] - VRFLI[i]      + INFLI[i] - aux_du_TRANI[i] - aux_du_DSFLI[i]
        # 4a) Isotope balance equation:
        # SOIL WATER
        ### THEORY:
        ### ∂/∂t Cᵢ     = [- Cᵢ/SWATI * ∂/∂t SWATI]           + 1/SWATI *
        ###                      [diff(z_upper) - diff(z_lower) - qCᵢ(z_upper) + qCᵢ(z_lower) +
        ###                       INFLI*Cᵢ_{INFLI} - TRANI*Cᵢ_{TRANI} - DSFL*Cᵢ_{DSFL} - SLVP*Cᵢ_{SLVP}]
        du_Cᵢ¹⁸_SWATI .= 0 # assert vector is zero
        du_Cᵢ²H_SWATI .= 0 # assert vector is zero
        du_Cᵢ¹⁸_SWATI .= -C_¹⁸Oᵏ./u_SWATIᵏ⁺¹ .* dVdt .+ 1 ./ u_SWATIᵏ⁺¹ .* (
                                -diff¹⁸O_upp*1000 .+ diff¹⁸O_low*1000 .+ qCᵢ¹⁸O_upp .- qCᵢ¹⁸O_low .+
                                aux_du_INFLI.*Cᵢ¹⁸O_INFLI .- aux_du_DSFLI.*Cᵢ¹⁸O_DSFL .- aux_du_SLVP[1].*C_¹⁸O_SLVP .-
                                aux_du_TRANI.*Cᵢ¹⁸O_TRANI
                            )
        du_Cᵢ²H_SWATI .= -C_²Hᵏ./u_SWATIᵏ⁺¹ .* dVdt .+ 1 ./ u_SWATIᵏ⁺¹ .* (
                                -diff²H_upp*1000 .+ diff²H_low*1000 .+ qCᵢ²H_upp .- qCᵢ²H_low .+
                                aux_du_INFLI.*Cᵢ²H_INFLI .- aux_du_DSFLI.*Cᵢ²H_DSFL .- aux_du_SLVP[1].*C_²H_SLVP .-
                                aux_du_TRANI.*Cᵢ²H_TRANI
                            )

    #TODO(bernhard): Diffusion slows the code down considerably and makes it unstable (for NLAYER=81 only, NLAYER=41 seems okay...)

        # # NOTE: below max(0.001,u_SWATIᵏ⁺¹) makes the code more robust
        # du_Cᵢ¹⁸_SWATI = -C_¹⁸Oᵏ./max.(0.001,u_SWATIᵏ⁺¹) .* dVdt .+ 1 ./ max.(0.001,u_SWATIᵏ⁺¹) .* (
        #                         -diff¹⁸O_upp*1000 .+ diff¹⁸O_low*1000 .+ qCᵢ¹⁸O_upp .- qCᵢ¹⁸O_low .+
        #                         aux_du_INFLI.*Cᵢ¹⁸O_INFLI .- aux_du_TRANI.*Cᵢ¹⁸O_TRANI .- aux_du_DSFLI.*Cᵢ¹⁸O_DSFL .- aux_du_SLVP[1].*C_¹⁸O_SLVP
        #                     )
        # du_Cᵢ²H_SWATI = -C_²Hᵏ./max.(0.001,u_SWATIᵏ⁺¹) .* dVdt .+ 1 ./ max.(0.001,u_SWATIᵏ⁺¹) .* (
        #                         -diff²H_upp*1000 .+ diff²H_low*1000 .+ qCᵢ²H_upp .- qCᵢ²H_low .+
        #                         aux_du_INFLI.*Cᵢ²H_INFLI .- aux_du_TRANI.*Cᵢ²H_TRANI .- aux_du_DSFLI.*Cᵢ²H_DSFL .- aux_du_SLVP[1].*C_²H_SLVP
        #                     )
        # GROUND WATER (GWAT)
        δ18O_empty = NaN
        δ2H_empty  = NaN
        @assert aux_du_VRFLI[NLAYER] >= 0 || isnan(aux_du_VRFLI[NLAYER]) "aux_du_VRFLI[NLAYER] should be positive" # First time this callback is called it can be NaN.
        # @assert aux_du_VRFLI[NLAYER] >= 0 "aux_du_VRFLI[NLAYER] should be positive" # First time this callback is called it can be NaN.

        if ((u_GWATᵏ⁺¹ == 0) & (aux_du_VRFLI[NLAYER] == 0)) # initially no groundwater and no new is added
            du_δ18O_GWAT = δ18O_empty
            du_δ2H_GWAT  = δ2H_empty
        elseif ((u_GWATᵏ⁺¹ == 0) & (aux_du_VRFLI[NLAYER] > 0)) # initially no groundwater but some is added
            # If no groundwater initially (u_δ_GWAT were NaN)
            # In that case override u to fixed value and set du to zero
            #TODO(bernhard): this is not working in f, we'd need to do it in the callback to set u: idx_u_scalar_isotopes_d18O = integrator.p[1][4][8   ]
            #TODO(bernhard): this is not working in f, we'd need to do it in the callback to set u: idx_u_scalar_isotopes_d2H  = integrator.p[1][4][10   ]
            #TODO(bernhard): this is not working in f, we'd need to do it in the callback to set u: integrator.u.GWATd18O = u_δ18O_SWATI[end] # u_δ18O_GWAT
            #TODO(bernhard): this is not working in f, we'd need to do it in the callback to set u: integrator.u.GWATd2H  = u_δ2H_SWATI[end]  # u_δ2H_GWAT

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
            V_GWAT = u_GWATᵏ⁺¹
            in_GWAT = aux_du_VRFLI[NLAYER]
            # out_GWAT = (du_GWFL + du_SEEP)

            # isotope balance
            Cᵢ¹⁸O_GWAT = LWFBrook90.ISO.δ_to_x.(u_δ18O_GWAT,LWFBrook90.ISO.R_VSMOW¹⁸O)
            Cᵢ²H_GWAT  = LWFBrook90.ISO.δ_to_x.(u_δ2H_GWAT,LWFBrook90.ISO.R_VSMOW²H )
            du_Cᵢ¹⁸_GWAT = in_GWAT/V_GWAT*(LWFBrook90.ISO.δ_to_x.(δ18O_in_GWAT, LWFBrook90.ISO.R_VSMOW¹⁸O) -  Cᵢ¹⁸O_GWAT) # - out_GWAT*(δ18O_out_GWAT - u_δ18O_GWAT) # <- last part == 0
            du_Cᵢ²H_GWAT = in_GWAT/V_GWAT*(LWFBrook90.ISO.δ_to_x.(δ2H_in_GWAT,  LWFBrook90.ISO.R_VSMOW²H ) -  Cᵢ²H_GWAT) # - out_GWAT*(δ2H_out_GWAT  - u_δ2H_GWAT)  # <- last part == 0

            # go back from atom fraction to delta values
            du_δ18O_GWAT  = LWFBrook90.ISO.dxdt_to_dδdt.(du_Cᵢ¹⁸_GWAT, Cᵢ¹⁸O_GWAT, LWFBrook90.ISO.R_VSMOW¹⁸O)
            du_δ2H_GWAT   = LWFBrook90.ISO.dxdt_to_dδdt.(du_Cᵢ²H_GWAT, Cᵢ²H_GWAT, LWFBrook90.ISO.R_VSMOW²H)
        end
        # go back from atom fraction to delta values
        du_δ18O_SWATI .= LWFBrook90.ISO.dxdt_to_dδdt.(du_Cᵢ¹⁸_SWATI, C_¹⁸Oᵏ, LWFBrook90.ISO.R_VSMOW¹⁸O)
        du_δ2H_SWATI  .= LWFBrook90.ISO.dxdt_to_dδdt.(du_Cᵢ²H_SWATI, C_²Hᵏ, LWFBrook90.ISO.R_VSMOW²H)

        # 5) Apply changes explicitly
        # Since we update this in the callback we can't overwrite `du` and let DiffEq.jl do
        # the work, so we need to update `u` instead of `du`
        # Using a simple forward Euler update:
        integrator.u.GWAT.d18O   += integrator.dt * du_δ18O_GWAT     #TODO(bernhard)
        integrator.u.GWAT.d2H    += integrator.dt * du_δ2H_GWAT      #TODO(bernhard)
        integrator.u.SWATI.d18O .+= integrator.dt * du_δ18O_SWATI    #TODO(bernhard)
        integrator.u.SWATI.d2H  .+= integrator.dt * du_δ2H_SWATI     #TODO(bernhard)
        ###
        ### end method from compute_isotope_GWAT_SWATI()


        # 6) Compute average composition of δRWU (weighted-mean of all uptake depths)
        # C¹⁸O_RWU = mean(Cᵢ¹⁸O_TRANI, weights(aux_du_TRANI))
        # C²H_RWU  = mean(Cᵢ²H_TRANI, weights(aux_du_TRANI))
        # # go back from atom fraction to delta values
        # δ¹⁸O_RWU = LWFBrook90.ISO.x_to_δ.(C¹⁸O_RWU, LWFBrook90.ISO.R_VSMOW¹⁸O)
        # δ²H_RWU  = LWFBrook90.ISO.x_to_δ.(C²H_RWU, LWFBrook90.ISO.R_VSMOW²H)
        δ¹⁸O_RWU = LWFBrook90.ISO.x_to_δ(mean(Cᵢ¹⁸O_TRANI, weights(aux_du_TRANI)), LWFBrook90.ISO.R_VSMOW¹⁸O)
        δ²H_RWU  = LWFBrook90.ISO.x_to_δ(mean(Cᵢ²H_TRANI,  weights(aux_du_TRANI)), LWFBrook90.ISO.R_VSMOW²H)

        integrator.u.RWU.d18O = δ¹⁸O_RWU
        integrator.u.RWU.d2H  = δ²H_RWU

        # 7) Update signature of Xylem reservoir
        # # Assume δXylem is a reservoir of volume Vxylem (initial conditions ...) and influx is
        # # sum(aux_du_TRANI) of composition δ_RWU and outflux of same volume but composition δ_Xylem:
        # # V_xylem * dδ_xylem/dt = sum(aux_du_TRANI) * (δ_RWU - δ_xylem)
        # # ==> dδ_xylem/dt = sum(aux_du_TRANI)/V_xylem * (δ_RWU - δ_xylem)
        # # ==>     dx/dt    = A * (B - x) = AB - Ax
        # # ==>     x = (x0 - B)*e^(-At) + B
        # # ==>     δ_xylem = δ_RWU + (δ_xylem - δ_RWU)*exp(-sum(aux_du_TRANI)/V_xylem * t)
        # NOTE that above equations are also correct when some elements of aux_du_TRANI are negative, i.e. outflow from xylem to soil water

        # current state and parameter V_xylem
        u_δ¹⁸O_Xylem = integrator.u.XYLEM.d18O
        u_δ²H_Xylem  = integrator.u.XYLEM.d2H

        # a) Approximate solution (Forward Euler integration over dt)
        # if δ¹⁸O_RWU is NaN, it means aux_du_TRANI is zero
        du_δ¹⁸O_Xylem = ifelse(isnan(δ¹⁸O_RWU), 0.0, sum(aux_du_TRANI) / p_VXYLEM * (δ¹⁸O_RWU - u_δ¹⁸O_Xylem))
        du_δ²H_Xylem  = ifelse(isnan(δ²H_RWU), 0.0, sum(aux_du_TRANI) / p_VXYLEM * (δ²H_RWU  - u_δ²H_Xylem))
        integrator.u.XYLEM.d18O += integrator.dt * du_δ¹⁸O_Xylem
        integrator.u.XYLEM.d2H  += integrator.dt * du_δ²H_Xylem
        # b) Precise solution (Analytical integration over dt)
        # u_δ¹⁸O_Xylem = u_δ18O_RWU + (u_δ¹⁸O_Xylem - u_δ18O_RWU)*exp(-sum(aux_du_TRANI)/p_VXYLEM * integrator.dt) #TODO: this should be using δ_to_x() and back
        # u_δ²H_Xylem  = u_δ2H_RWU  + (u_δ²H_Xylem  - u_δ2H_RWU )*exp(-sum(aux_du_TRANI)/p_VXYLEM * integrator.dt) #TODO: this should be using δ_to_x() and back
    end

    return nothing
end




function LWFBrook90R_check_balance_errors!(integrator)
    if integrator.p.compute_intermediate_quantities
        # Compute daily water balance errors
        u_GWAT     = integrator.u.GWAT.mm
        u_INTS     = integrator.u.INTS.mm
        u_INTR     = integrator.u.INTR.mm
        u_SNOW     = integrator.u.SNOW.mm
        # u_CC     = integrator.u.CC.MJm2
        # u_SNOWLQ = integrator.u.SNOWLQ.mm
        u_SWATI    = integrator.u.SWATI.mm

        # a) Get change in total storages
        # Get old total water volumes and compute new total water volumes
        # old_t          = integrator.uprev.accum.ε_prev_t
        # old_SWAT       = integrator.uprev.accum.ε_prev_totalSWAT
        # old_totalWATER = integrator.uprev.accum.ε_prev_totalWATER
        old_SWAT       = integrator.uprev.accum.totalSWAT
        old_totalWATER = integrator.uprev.accum.new_totalWATER
        new_totalSWAT  = sum(u_SWATI)                                      # total soil water in all layers, mm
        new_totalWATER = u_INTR + u_INTS + u_SNOW + new_totalSWAT + u_GWAT # total water in all compartments, mm

        # b) Get fluxes that caused change in total storages
        # Across outer boundaries of above and belowground system:
        # b1) BALERD_SWAT  = old_SWAT - new_totalSWAT + INFLI - DSFLI - VRFL(NLAYER) - TRANI - SLVP
        # accum_qin corresponds to q(t,0)  in Ireson 2023 eq 16 with additionally sources and sinks:
        accum_qin  = 0 +  (integrator.u.accum.slfl - integrator.u.accum.byfl)
        # accum_qout corresponds to q(t,zN) in Ireson 2023 eq 16 with additionally sources and sinks:
        accum_qout = integrator.u.accum.vrfln +
            integrator.u.accum.dsfl +
            integrator.u.accum.cum_d_tran +
            integrator.u.accum.cum_d_slvp

        # b2) BALERD_total = old_totalWATER - new_totalWATER + PRECD - EVAPD - FLOWD - SEEPD
        accum_qin_total  = integrator.u.accum.cum_d_prec
        accum_qout_total = integrator.u.accum.cum_d_evap + integrator.u.accum.flow + integrator.u.accum.seep

        # c) Compute balance error
        # d) Store balance errors into state vector
        # BALERD_total = old_totalWATER - new_totalWATER + PRECD - EVAPD - FLOWD - SEEPD

        integrator.u.accum.BALERD_SWAT  = (accum_qin - accum_qout)             - (new_totalSWAT  - old_SWAT)       # implicitly assume accum.XXX are reset to zero daily AND this callback is called daily
        integrator.u.accum.BALERD_total = (accum_qin_total - accum_qout_total) - (new_totalWATER - old_totalWATER) # implicitly assume accum.XXX are reset to zero daily AND this callback is called daily
        # alternatively:
        # integrator.u.accum.BALERD_SWAT  = (accum_qin - accum_qout) * (integrator.t - old_t)             - (new_totalSWAT  - old_SWAT)
        # integrator.u.accum.BALERD_total = (accum_qin_total - accum_qout_total) * (integrator.t - old_t) - (new_totalWATER - old_totalWATER)

        # Save current totals into caches to check error for next period
        integrator.u.accum.totalSWAT    = new_totalSWAT
        integrator.u.accum.new_totalWATER   = new_totalWATER
        # integrator.u.accum.ε_prev_t            = integrator.t
        # integrator.u.accum.ε_prev_totalSWAT    = new_totalSWAT
        # integrator.u.accum.ε_prev_totalWATER   = new_totalWATER
        # integrator.u.accum.ε_prev_qcumInSWAT   = accum_qin        # no need to store these
        # integrator.u.accum.ε_prev_qcumOutSWAT  = accum_qout       # no need to store these
        # integrator.u.accum.ε_prev_qcumInWater  = accum_qin_total  # no need to store these
        # integrator.u.accum.ε_prev_qcumOutWater = accum_qout_total # no need to store these
    end

    return nothing
end
