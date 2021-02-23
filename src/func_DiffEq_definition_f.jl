function define_DiffEq_timestep_cb()
    # A) Define updating function
    function LWFBrook90R_compute_RHS_and_timestep!(u,t,integrator)
        # if integrator.t>=247.8 && integrator.t<251
        #     @info "Timestep callback. t: $(integrator.t), DTRI: $(integrator.p[4][1]), Integrator.dt: $(get_proposed_dt(integrator)), u_SWATI_1:$(integrator.u[7])"
        # end
        # NOTE: we can make use of those:
        # integrator.t
        # integrator.p
        # integrator.u

        ##################
        # Parse parameters
        ## A) constant parameters:
        (p_DT, NLAYER, IMODEL, compute_intermediate_quantities, Reset,
        p_SWATMX, p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_PSIG, p_KF,
        p_THSAT, p_θr, p_MvGα, p_MvGn, p_MvGl, p_Ksat,

        # FOR MSBITERATE:
        p_QLAYER, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC, p_IMPERV,
        p_LENGTH, p_DSLOPE, p_RHOWG, p_DPSIMX, #TODO(bernhard) p_RHOWG is a global constant
        p_KSAT, p_DRAIN, p_DTIMAX, p_INFRAC, p_DSWMAX,
        p_GSC, p_GSP, p_THICK, p_STONEF,

        p_BYPAR) = integrator.p[1][1]

        # unused are the constant parameters saved in: = integrator.p[1][2]

        ## B) time dependent parameters
        (p_DOY, p_MONTHN, p_SOLRAD, p_TMAX, p_TMIN, p_EA, p_UW, p_PRECIN, p_DTP, p_NPINT, p_MESFL,
        _, _, _, _, _) = integrator.p[2]

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

            p_fT_PREINT = p_PRECIN(t) / p_DTP # (mm/day)
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

        # These were computed in the callback and are kept constant in between two
        # callbacks.
        p_fu_RNET, aux_du_SMLT, aux_du_TRANI, aux_du_SLVP = integrator.p[3]


        ##################
        # Parse states
        u_GWAT     = integrator.u[1]
        #u_INTS     = integrator.u[2]
        #u_INTR     = integrator.u[3]
        #u_SNOW     = integrator.u[4]
        #u_CC       = integrator.u[5]
        #u_SNOWLQ   = integrator.u[6]
        u_SWATI    = integrator.u[7:(7+NLAYER-1)]

        (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
            LWFBrook90Julia.KPT.derive_auxiliary_SOILVAR(u_SWATI, p_SWATMX, p_THSAT,
                 p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_KF,
                 p_θr, p_MvGα, p_MvGn, p_MvGl, p_Ksat,
                 p_PSIG, NLAYER, IMODEL)

        ##################
        # Compute fluxes

        # Bypass fraction of infiltration to each layer
        p_fu_BYFRAC = LWFBrook90Julia.WAT.BYFLFR(
                      NLAYER, p_BYPAR, p_QFPAR, p_QFFC, u_aux_WETNES, p_WETF)

        # Assert time step is smaller than remaining DTRI in order to avoid overshooting a day
        # Variant 1) Derive from time t
        # DTRI = ceil(t; digits=0) - t
        # if (DTRI == 0)
        #     DTRI = p_DTP
        # end
        # Variant 2) Derive from saved DTRI
        DTRI = integrator.p[4][1]
        if (DTRI <= 0)
            DTRI = p_DTP # TODO(bernhard): sometimes DTRI == 0, resulting in division by DTI=0
        end
        if (integrator.opts.adaptive)
            error("DTRI is not working with adaptive solvers from DiffEq.jl")
        end

        # first approximation for iteration time step, time remaining or DTIMAX
        DTI  = min(DTRI, p_DTIMAX)

        # Water movement through soil
        (p_fu_SRFL, p_fu_SLFL, aux_du_DSFLI, aux_du_VRFLI, aux_du_INFLI, aux_du_BYFLI,
        du_NTFLI, du_GWFL, du_SEEP, DTINEW, aux_du_TRANI_corrected, aux_du_SLVP_corrected) =
            MSBITERATE(IMODEL, p_QLAYER,
                    # for SRFLFR:
                    u_SWATI, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC,
                    #
                    p_IMPERV, p_fu_RNET, aux_du_SMLT, NLAYER,
                    p_LENGTH, p_DSLOPE,
                    # for DSLOP:
                    p_RHOWG, u_aux_PSIM, p_THICK, p_STONEF, p_fu_KK,
                    #
                    u_aux_PSITI, p_DPSIMX,
                    # for VERT:
                    p_KSAT,
                    #
                    p_DRAIN, DTI,
                    # for INFLOW:
                    p_INFRAC, p_fu_BYFRAC, aux_du_TRANI, aux_du_SLVP, p_SWATMX,
                    # for ITER:
                    u_aux_θ, u_aux_WETNES,
                    p_DSWMAX, p_THSAT, p_θr, p_BEXP, p_PSIF, p_WETF, p_CHM, p_CHN, p_WETINF, p_MvGα, p_MvGn,
                    # for GWATER:
                    u_GWAT, p_GSC, p_GSP, p_DT)

        # Transport flow (heat, solutes, isotopes, ...)
        # TODO(bernhard): see initial prototype code... (script3)

        # DEBUG
        #@info """t:$(round(t;digits=4)), DTI:$DTINEW, DTRI:$DTRI, NITS:$(integrator.u[7+NLAYER+25]), sum(aux_du_BYFLI):$(sum(aux_du_BYFLI)), sum(aux_du_DSFLI):$(sum(aux_du_DSFLI))""" # , p_fu_SRFL:$p_fu_SRFL, du_GWFL:$du_GWFL
        # END DEBUG

        # save intermediate results for use in ODE (function f())
        integrator.p[4][2] = aux_du_VRFLI[NLAYER]
        integrator.p[4][3] = du_GWFL
        integrator.p[4][4] = du_SEEP
        integrator.p[4][5] = p_fu_SRFL
        integrator.p[4][6] = p_fu_SLFL
        integrator.p[4][7] = du_NTFLI
        integrator.p[4][8] = aux_du_BYFLI
        integrator.p[4][9] = aux_du_DSFLI

        # Force next time step to be: DTINEW
        integrator.dtcache = DTINEW # https://github.com/SciML/OrdinaryDiffEq.jl/issues/1351
        set_proposed_dt!(integrator, DTINEW)
        # Update DTRI:
        integrator.p[4][1] -= DTINEW
        # Update NITS:
        integrator.u[7+NLAYER+25] += 1 # cum_d_nits

        # Apply corrected SLVP and TRANI
        integrator.u[7+NLAYER+ 8] += DTINEW*(aux_du_SLVP_corrected + sum(aux_du_TRANI_corrected))  # cum_d_evap
        integrator.u[7+NLAYER+ 9] += DTINEW*(                        sum(aux_du_TRANI_corrected))
        integrator.u[7+NLAYER+12] += DTINEW*(aux_du_SLVP_corrected)
        integrator.p[3][3] = aux_du_TRANI_corrected
        integrator.p[3][4] = aux_du_SLVP_corrected

        ##################
        # Update soil limited boundary flows during iteration loop
        if Reset == 1
            #error("The case with updated flows (Reset==1) is not implemented here. It was initially implemented in LWFB90V4.jl")

            if (integrator.u[7+NLAYER+25] > 0)
                # starting from the 2nd iteration
                # update the following quantities based on recomputed state variables
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

                # 3) Update: p_fu_RNET, p_fu_PTRAN, aux_du_TRANI, aux_du_SLVP, aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, aux_du_RSNO, aux_du_SNVP, aux_du_SMLT,
                #            u_SNOW, u_CC, u_SNOWLQ
                #    using recomputed:     u_aux_PSIM, p_fu_PSLVP, p_fu_PINT, p_fu_PTRAN, aux_du_TRANI, p_fu_GIVP, p_fu_GEVP
                #    Note that it is
                #    using not recomputed: p_fu_TA, p_fu_LAI, p_fu_SAI, p_fu_PSNVP, p_fu_SNOEN,
                #                          u_INTR, u_SNOW, u_CC, u_SNOWLQ,

                (# compute some fluxes as intermediate results:
                    #p_fT_SFAL, p_fT_RFAL, p_fu_RNET, p_fu_PTRAN,
                    p_fT_SFAL, p_fT_RFAL, _, p_fu_PTRAN, # NOTE: Don't modify these variable as this is the reset run!
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
            end
        end
    end

    # B) Define callback
    cb_func = FunctionCallingCallback(LWFBrook90R_compute_RHS_and_timestep!;
                 func_everystep=true, func_start = true)
    return cb_func
end

""" define_DiffEq_f()\n
    Generates function f (right-hand-side of ODEs) needed for ODE() problem in DiffEq.jl package.
"""
# function define_DiffEq_f()
    function f_LWFBrook90R(du,u,p,t)
        # computes right hand side (du) of states: u_GWAT, u_SWATI

        # NOTE: Within the ufnction we can make use of:
        # t
        # p
        # u

        ##################
        # Parse parameters
        ## A) constant parameters:
        (p_DT, NLAYER, IMODEL, compute_intermediate_quantities, Reset) = p[1][1]

        ## A) solution depenedent parameters (computed in callback):
        (DTRI, # DTRI unused
        aux_du_VRFLI__NLAYER, du_GWFL, du_SEEP, p_fu_SRFL, p_fu_SLFL,
        du_NTFLI,
        aux_du_BYFLI,
        aux_du_DSFLI,) = p[4]

        ##################
        # Update solution:
        # u is a state vector with u[1] = S relative saturation (-)
        # Update GWAT:
        du[1] = aux_du_VRFLI__NLAYER - du_GWFL - du_SEEP
        # Update SWATI for each layer:
        #du[7] = 0
        #du[8] = 0
        #du[9] = 0
        #du[10] = 0
        #du[11] = 0
        du[7:(7+NLAYER-1)] = du_NTFLI[:]

        # Do not modify INTS, INTR, SNOW, CC, SNOWLQ
        # as they are separately modified by the callback.
        # Therfore set their rate of change to 0
        du[2] = 0
        du[3] = 0
        du[4] = 0
        du[5] = 0
        du[6] = 0

        ##########################################
        # Accumulate flows to compute daily sums
        # Note that below state variables serve only as accumulator but do not affect
        # the evolution of the system.
        if compute_intermediate_quantities
            # 1) Either set daily sum if rate is constant throughout precipitation interval: p_DTP*(...)
            # 2) or then set daily sum to zero and use ODE to accumulate flow.
            du[7+NLAYER+0] = 0 # cum_d_prec, was computed in callback
            du[7+NLAYER+1] = 0 # cum_d_rfal, was computed in callback
            du[7+NLAYER+2] = 0 # cum_d_sfal, was computed in callback
            du[7+NLAYER+3] = 0 # cum_d_rint, was computed in callback
            du[7+NLAYER+4] = 0 # cum_d_sint, was computed in callback
            du[7+NLAYER+5] = 0 # cum_d_rsno, was computed in callback
            du[7+NLAYER+6] = 0 # cum_d_rnet, was computed in callback
            du[7+NLAYER+7] = 0 # cum_d_smlt, was computed in callback

            du[7+NLAYER+ 8] = 0 # cum_d_evap,  was computed in callback
            du[7+NLAYER+ 9] = 0 # cum_d_tran,  was computed in callback
            du[7+NLAYER+10] = 0 # cum_d_irvp,  was computed in callback
            du[7+NLAYER+11] = 0 # cum_d_isvp,  was computed in callback
            du[7+NLAYER+12] = 0 # cum_d_slvp,  was computed in callback
            du[7+NLAYER+13] = 0 # cum_d_snvp,  was computed in callback
            du[7+NLAYER+14] = 0 # cum_d_pint,  was computed in callback
            du[7+NLAYER+15] = 0 # cum_d_ptran, was computed in callback
            du[7+NLAYER+16] = 0 # cum_d_pslvp, was computed in callback

            du[7+NLAYER+17] = p_fu_SRFL +
                              sum(aux_du_BYFLI) +
                              sum(aux_du_DSFLI) +
                              du_GWFL # SRFLD + BYFLD + DSFLD + GWFLD # flow
            du[7+NLAYER+18] = du_SEEP                                 # seep
            du[7+NLAYER+19] = p_fu_SRFL                               # srfl
            du[7+NLAYER+20] = p_fu_SLFL                               # slfl
            du[7+NLAYER+21] = sum(aux_du_BYFLI)                       # byfl
            du[7+NLAYER+22] = sum(aux_du_DSFLI)                       # dsfl
            du[7+NLAYER+23] = du_GWFL                                 # gwfl
            du[7+NLAYER+24] = aux_du_VRFLI__NLAYER                    # vrfln

            # du[7+NLAYER+25]= 0 # cum_d_nits, was computed in callback
            # du[7+NLAYER+26]= 0 # cum_d_rthr, was computed in callback
            # du[7+NLAYER+27]= 0 # cum_d_sthr, was computed in callback
            # du[7+NLAYER+28]= 0 # mesfl, was computed in callback
            # balerd[IDAY]=BALERD

            # TODO(bernhard): use SavingCallback() for all quantities that have du=0
            #                 only keep du=... for quantities for which we compute cumulative sums
        end

        return
    end
#     f_func = f_LWFBrook90R
#     return f_func
# end
