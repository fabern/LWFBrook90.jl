"""
    define_LWFB90_cb()

Generate callback function cb needed for ODE() problem in DiffEq.jl package.

LWFBrook90 updates states INTS, INTR, SNOW, CC, SNOWLQ not continuously but only
once per day. This operator splitting (daily vs continuous update of ODEs) is
implemented by using this callback function which is called once per day.
"""
function define_LWFB90_cb()
    # A) Define updating function
    function LWFBrook90R_update_INTS_INTR_SNOW_CC_SNOWLQ!(integrator)
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
            p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_AGE, p_RELDEN = integrator.p[2]

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

        idx_u_vector_amounts       = integrator.p[1][4][4] # 7:(6+NLAYER)
        idx_u_vector_accumulators  = integrator.p[1][4][5] # 6+0+0+1*NLAYER .+ (1:25)

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

        p_fT_DAYLEN, p_fT_I0HDAY, p_fT_SLFDAY, p_fu_HEIGHT, p_fu_LAI, p_fu_SAI, p_RTLEN,
          p_RPLANT,p_fu_Z0GS, p_fu_Z0C, p_fu_DISPC, p_fu_Z0, p_fu_DISP, p_fu_ZA,
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
        # Combine day and night rates to average daily rate
        (p_fu_PTRAN, p_fu_GEVP, p_fu_PINT, p_fu_GIVP, p_fu_PSLVP, aux_du_TRANI) = # TODO(bernhard): p_fu_PSLVP is unused
            MSBDAYNIGHT_postprocess(FLAG_MualVanGen, NLAYER, p_fu_PTR, p_fu_GER, p_fu_PIR, p_fu_GIR, p_fu_ATRI, p_fT_DAYLEN, p_fu_PGER)
        #* * * * * * * *  E N D   D A Y - N I G H T   L O O P  * * * * * * * * * *

        ####################################################################
        # 1) Update snow accumulation/melt: u_SNOW, u_CC, u_SNOWLQ
        #    and compute fluxes to/from interception storage
        (# compute some fluxes as intermediate results:
        p_fT_SFAL, p_fT_RFAL, p_fu_RNET, p_fu_PTRAN,
        # compute changes in soil water storage:
        aux_du_TRANI, aux_du_SLVP,
        # compute change in interception storage:
        aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP,
        # compute change in snow storage:
        aux_du_RSNO, aux_du_SNVP, aux_du_SMLT,
        # compute updated states:
        u_SNOW, u_CC, u_SNOWLQ) =
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
                   u_CC, u_SNOWLQ, p_fu_PSNVP, p_fu_SNOEN, p_MAXLQF, p_GRDMLT)


        ####################################################################
        # 2) Update states of interception storage over entire precipitation
        #    interval: u_INTS, u_INTR
        u_INTS = u_INTS + (aux_du_SINT - aux_du_ISVP) * p_DTP
        u_INTR = u_INTR + (aux_du_RINT - aux_du_IRVP) * p_DTP

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

        # save intermediate results for use in ODE (function f())
        # integrator.p[3] = [p_fu_RNET, aux_du_SMLT, aux_du_TRANI, aux_du_SLVP]
        integrator.p[3][1] = p_fu_RNET
        integrator.p[3][2] = aux_du_SMLT
        integrator.p[3][3] = aux_du_TRANI
        integrator.p[3][4] = aux_du_SLVP



        ##########################################
        # Accumulate flows to compute daily sums
        # Note that below state variables serve only as accumulator but do not affect
        # the evolution of the system.
        if compute_intermediate_quantities
            # @assert length(idx_u_vector_accumulators) == 25

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
            # timeseries_balerd[IDAY]=BALERD

            # TODO(bernhard): use SavingCallback() for all quantities that have u=... and du=0
            #                 only keep du=... for quantities for which we compute cumulative sums
        end
        ##########################################
    end

    # B) Define callback
    cb_func = PeriodicCallback(LWFBrook90R_update_INTS_INTR_SNOW_CC_SNOWLQ!,  1.0;
                               initial_affect = true,
                               save_positions=(false,false));

    return cb_func
end