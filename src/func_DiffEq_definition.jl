""" define_DiffEq_f()\n
    Generates function f (right-hand-side of ODEs) needed for ODE() probelm in DiffEq.jl package.
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
        (p_DT, NLAYER, IMODEL, compute_intermediate_quantities, Reset,
        p_SWATMX, p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_PSIG, p_KF,
        p_THSAT, p_θr, p_MvGα, p_MvGn, p_MvGl, p_Ksat,

        # FOR MSBITERATE:
        p_QLAYER, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC, p_IMPERV,
        p_LENGTH, p_DSLOPE, p_RHOWG, p_DPSIMX, #TODO(bernhard) p_RHOWG is a global constant
        p_KSAT, p_DRAIN, p_DTIMAX, p_INFRAC, p_DSWMAX,
        p_GSC, p_GSP, p_THICK, p_STONEF,

        p_BYPAR) = p[1][1]

        # unused are the constant parameters saved in: = p[1][2]

        ## B) time dependent parameters
        (p_DOY, p_MONTHN, p_SOLRAD, p_TMAX, p_TMIN, p_EA, p_UW, p_PRECIN, p_DTP, p_NPINT, p_MESFL,
        _, _, _, _, _) = p[2]

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

        DTRI = p_DTP

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
        p_fu_RNET, aux_du_SMLT, aux_du_TRANI, aux_du_SLVP = p[3]


        ##################
        # Parse states
        u_GWAT     = u[1]
        #u_INTS     = u[2]
        #u_INTR     = u[3]
        #u_SNOW     = u[4]
        #u_CC       = u[5]
        #u_SNOWLQ   = u[6]
        u_SWATI    = u[7:(7+NLAYER-1)]

        (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
            LWFBrook90Julia.KPT.derive_auxiliary_SOILVAR(u_SWATI, p_SWATMX, p_THSAT,
                 p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_KF,
                 p_θr, p_MvGα, p_MvGn, p_MvGl, p_Ksat,
                 p_PSIG, NLAYER, IMODEL)

        ##################
        # Update soil limited boundary flows during iteration loop
        if Reset == 1
            error("The case with updated flows (Reset==1) is not implemented here. It was initially implemented in LWFB90V4.jl")
        end

        ##################
        # Compute fluxes

        # Bypass fraction of infiltration to each layer
        p_fu_BYFRAC = LWFBrook90Julia.WAT.BYFLFR(
                      NLAYER, p_BYPAR, p_QFPAR, p_QFFC, u_aux_WETNES, p_WETF)

        # Water movement through soil
        (p_fu_SRFL, p_fu_SLFL, aux_du_DSFLI, aux_du_VRFLI, DTI, aux_du_INFLI, aux_du_BYFLI,
        du_NTFLI, du_GWFL, du_SEEP) =
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
                    p_DRAIN, DTRI, p_DTIMAX,
                    # for INFLOW:
                    p_INFRAC, p_fu_BYFRAC, aux_du_TRANI, aux_du_SLVP, p_SWATMX,
                    # for FDPSIDW:
                    u_aux_WETNES, p_BEXP, p_PSIF, p_WETF, p_CHM, p_CHN, p_MvGα, p_MvGn,
                    # for ITER:
                    p_DSWMAX, p_THSAT, p_θr, u_aux_θ,
                    # for GWATER:
                    u_GWAT, p_GSC, p_GSP, p_DT)

        # Transport flow (heat, solutes, isotopes, ...)
        # TODO(bernhard): see initial prototype code... (script3)

        ##################
        # Update solution:
        # u is a state vector with u[1] = S relative saturation (-)
        # Update GWAT:
        du[1] = aux_du_VRFLI[NLAYER] - du_GWFL - du_SEEP
        # Update SWATI for each layer:
        #du[7] = 0
        #du[8] = 0
        #du[9] = 0
        #du[10] = 0
        #du[11] = 0
        du[7:(7+(NLAYER-1))] = du_NTFLI[:]

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
            du[7+NLAYER+0] = 0 # was computed in callback         # was timeseries_prec
            du[7+NLAYER+1] = 0 # was computed in callback         # was timeseries_evp
            du[7+NLAYER+2] = p_fu_SRFL + sum(aux_du_BYFLI) + sum(aux_du_DSFLI) + du_GWFL # SRFLD + BYFLD + DSFLD + GWFLD # was timeseries_flow
            du[7+NLAYER+3] = 0 # was computed in callback         # was timeseries_rnet
            du[7+NLAYER+4] = 0 # was computed in callback         # was timeseries_irvp
            du[7+NLAYER+5] = 0 # was computed in callback         # was timeseries_isvp
            du[7+NLAYER+6] = 0 # was computed in callback         # was timeseries_ptran
            du[7+NLAYER+7] = 0 # was computed in callback         # was timeseries_pint
            du[7+NLAYER+8] = 0 # was computed in callback         # was timeseries_snvp
            du[7+NLAYER+9] = 0 # was computed in callback         # was timeseries_slvp
            du[7+NLAYER+10]= 0 # was computed in callback         # was timeseries_trand
            du[7+NLAYER+11]= 0 # was computed in callback         # was timeseries_mesfld
            du[7+NLAYER+12]= 0 # was computed in callback         # was timeseries_smltd
            du[7+NLAYER+13]= p_fu_SLFL                            # was timeseries_slfld
            du[7+NLAYER+14]= 0 # was computed in callback         # was timeseries_rfald
            du[7+NLAYER+15]= 0 # was computed in callback         # was timeseries_sfald
            #awat[IDAY]=AWAT # TODO(bernhard): AWAT and ADEF are never overwritten and remain NaN as initialized
            #adef[IDAY]=ADEF # TODO(bernhard): AWAT and ADEF are never overwritten and remain NaN as initialized
            du[7+NLAYER+16]= 0 # was computed in callback         # was timeseries_sintd
            du[7+NLAYER+17]= 0 # was computed in callback         # was timeseries_rintd
            du[7+NLAYER+18]= 0 # was computed in callback         # was timeseries_cum_d_rthr
            du[7+NLAYER+19]= 0 # was computed in callback         # was timeseries_cum_d_sthr
            du[7+NLAYER+20]= 0 # was computed in callback         # was timeseries_rsnod
            # balerd[IDAY]=BALERD

            du[7+NLAYER+21]= p_fu_SRFL                            # was timeseries_SRFLD
            du[7+NLAYER+22]= du_SEEP                              # was timeseries_SEEP
            du[7+NLAYER+23]= du_GWFL                              # was timeseries_GWFL
            du[7+NLAYER+24]= aux_du_VRFLI[NLAYER]                 # was timeseries_VRFLN
            du[7+NLAYER+25]= sum(aux_du_BYFLI)                    # was timeseries_BYFLD
            du[7+NLAYER+26]= sum(aux_du_DSFLI)                    # was timeseries_DSFLD

            # TODO(bernhard): use SavingCallback() for all quantities that have du=0
            #                 only keep du=... for quantities for which we compute cumulative sums
        end

        return
    end
#     f_func = f_LWFBrook90R
#     return f_func
# end

""" define_DiffEq_cb()\n
    Generates callback function cb needed for ODE() probelm in DiffEq.jl package.
    LWFBrook90 updates states INTS, INTR, SNOW, CC, SNOWLQ not continuously but only
    once per day. This operator splitting (daily vs continuous update of ODEs) is
    implemented by using this callback function which is called once per day.
"""
function define_DiffEq_cb()
    # A) Define updating function
    function LWFBrook90R_update_INTS_INTR_SNOW_CC_SNOWLQ!(integrator)
        # NOTE: we can make use of those:
        # integrator.t
        # integrator.p
        # integrator.u

        ############
        ### Compute parameters
        ## A) constant parameters:
        (p_DT, NLAYER, IMODEL, compute_intermediate_quantities, Reset,
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
            p_MESFL, p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_AGE = integrator.p[2]

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

        IDAY = floor(integrator.t) # TODO(bernhard) is just for debug, remove again after

        p_fT_DAYLEN, p_fT_I0HDAY, p_fT_SLFDAY, p_fu_HEIGHT, p_fu_LAI, p_fu_SAI, p_fu_RTLEN,
          p_fu_RPLANT,p_fu_Z0GS, p_fu_Z0C, p_fu_DISPC, p_fu_Z0, p_fu_DISP, p_fu_ZA,
          p_fu_RXYLEM, p_fu_RROOTI, p_fu_ALPHA,
          p_fu_SHEAT,p_fu_SOLRADC, p_fu_TA, p_fu_TADTM, p_fu_TANTM, p_fu_UADTM, p_fu_UANTM,
          p_fT_SNOFRC,p_fu_TSNOW,p_fu_PSNVP, p_fu_ALBEDO,p_fu_RSS, p_fu_SNOEN =
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


        # Calculate average daily rate of potential and actual interception,
        # evaporation, and transpiration by considering weighted average of rate
        # during day and rate during night:
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
        (p_fu_PTRAN, p_fu_GEVP, p_fu_PINT, p_fu_GIVP, p_fu_PSLVP, aux_du_TRANI) = # TODO(bernhard): p_fu_PSLVP is unused
            MSBDAYNIGHT_postprocess(IMODEL, NLAYER, p_fu_PTR, p_fu_GER, p_fu_PIR, p_fu_GIR, p_fu_ATRI, p_fT_DAYLEN, p_fu_PGER, p_DT)
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
            # 1) Either set daily sum if rate is constant throughout precipitation interval: p_DTP*(...)
            # 2) or then set daily sum to zero and use ODE to accumulate flow.
            integrator.u[7+NLAYER+0] = p_DTP*(p_fT_RFAL + p_fT_SFAL)   # RFALD + SFALD
            integrator.u[7+NLAYER+1] = p_DTP*(aux_du_IRVP + aux_du_ISVP + aux_du_SNVP + aux_du_SLVP + sum(aux_du_TRANI))
            integrator.u[7+NLAYER+2] = 0 # du/dt = p_fu_SRFL + sum(aux_du_BYFLI) + sum(aux_du_DSFLI) + du_GWFL
            integrator.u[7+NLAYER+3] = p_DTP*(p_fT_RFAL - aux_du_RINT - aux_du_RSNO) # cum_d_RTHR - RSNOD
            integrator.u[7+NLAYER+4] = p_DTP*(aux_du_IRVP)
            integrator.u[7+NLAYER+5] = p_DTP*(aux_du_ISVP)
            integrator.u[7+NLAYER+6] = p_DTP*(p_fu_PTRAN)
            integrator.u[7+NLAYER+7] = p_DTP*(p_fu_PINT)
            integrator.u[7+NLAYER+8] = p_DTP*(aux_du_SNVP)
            integrator.u[7+NLAYER+9] = p_DTP*(aux_du_SLVP)
            integrator.u[7+NLAYER+10]= p_DTP*(sum(aux_du_TRANI))
            integrator.u[7+NLAYER+11]= p_DTP*(p_MESFL(integrator.t))
            integrator.u[7+NLAYER+12]= p_DTP*(aux_du_SMLT)
            integrator.u[7+NLAYER+13]= 0 # du/dt = p_fu_SLFL
            integrator.u[7+NLAYER+14]= p_DTP*(p_fT_RFAL)
            integrator.u[7+NLAYER+15]= p_DTP*(p_fT_SFAL)
            #timeseries_awat[IDAY]=AWAT # TODO(bernhard): AWAT and ADEF are never overwritten and remain NaN as initialized
            #timeseries_adef[IDAY]=ADEF # TODO(bernhard): AWAT and ADEF are never overwritten and remain NaN as initialized
            integrator.u[7+NLAYER+16]= p_DTP*(aux_du_SINT)
            integrator.u[7+NLAYER+17]= p_DTP*(aux_du_RINT)
            integrator.u[7+NLAYER+18]= p_DTP*(p_fT_RFAL - aux_du_RINT)
            integrator.u[7+NLAYER+19]= p_DTP*(p_fT_SFAL - aux_du_SINT)
            integrator.u[7+NLAYER+20]= p_DTP*(aux_du_RSNO)
            # timeseries_balerd[IDAY]=BALERD
            integrator.u[7+NLAYER+21]= 0 # du/dt = p_fu_SRFL
            integrator.u[7+NLAYER+22]= 0 # du/dt = du_SEEP
            integrator.u[7+NLAYER+23]= 0 # du/dt = du_GWFL
            integrator.u[7+NLAYER+24]= 0 # du/dt = aux_du_VRFLI[NLAYER]
            integrator.u[7+NLAYER+25]= 0 # du/dt = sum(aux_du_BYFLI)
            integrator.u[7+NLAYER+26]= 0 # du/dt = sum(aux_du_DSFLI)

            # TODO(bernhard): use SavingCallback() for all quantities that have u=... and du=0
            #                 only keep du=... for quantities for which we compute cumulative sums
        end
        ##########################################
    end

    # B) Define callback
    cb_func = PeriodicCallback(LWFBrook90R_update_INTS_INTR_SNOW_CC_SNOWLQ!,  1.0;
                               initial_affect = true);
    return cb_func
end

""" define_DiffEq_ODE()\n
    Generates an ODEProblem from DiffEq.jl.\n\n
    # An ODE problem which consists of
    #   - definition of right-hand-side (RHS) function f
    #   - definition of callback function cb
    #   - initial condition of states
    #   - definition of simulation time span
    #   - parameters\n
    Seperate updating of different states (INTS, INTR, SNOW, CC, SNOWLQ are updated once per
    day while GWAT and SWATI are updated continuously) is implemented by means of operator
    splitting using a callback function for the daily updates and a ODE RHS (right hand
    side) for the continuous update.
"""
function define_DiffEq_ODE(u0, tspan, p)

    # Define callback functions
    cb_func = define_DiffEq_cb()

    # swcheck_cb = ContinuousCallback()    #TODO(bernhard) Implement swchek as ContinuousCallback
    # Reset_cb = FunctionCallingCallback() #TODO(bernhard) low priority: implement Reset==1

    # Define ODE problem
    ode = ODEProblem(f_LWFBrook90R,
                     u0,
                     tspan,
                     p,
                     callback=cb_func)
    return ode
end

""" define_DiffEq_u0()\n
    Generates vector u0 needed for ODE() problem in DiffEq.jl package."""
function define_DiffEq_u0(u_GWAT_init,
                            u_INTS_init,
                            u_INTR_init,
                            u_SNOW_init,
                            u_CC_init,
                            u_SNOWLQ_init,
                            u_SWATIinit,
                            compute_intermediate_quantities)

    if compute_intermediate_quantities
    # TODO(bernhard): are these store somewhere else than input_siteparam and pfile_param
        u0 = [u_GWAT_init;
            u_INTS_init;
            u_INTR_init;
            u_SNOW_init;
            u_CC_init;
            u_SNOWLQ_init;
            u_SWATIinit;
            # accumulation variables:
            0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;
            0;0;0;0;0;0]
    else
        u0 = [u_GWAT_init;
                u_INTS_init;
                u_INTR_init;
                u_SNOW_init;
                u_CC_init;
                u_SNOWLQ_init;
                u_SWATIinit]
    end

    return u0
end


""" define_diff_eq_parameters(NLAYER, IMODEL, constant_dt_solver, NOOUTF, Reset, compute_intermediate_quantities,
    pfile_meteo, pfile_siteparam, pfile_param, pfile_soil, pfile_pdur)\n
    Generates vector p needed for ODE() problem in DiffEq.jl package."""
function define_DiffEq_parameters(NLAYER, IMODEL, constant_dt_solver, NOOUTF, Reset, compute_intermediate_quantities,
    pfile_meteo, pfile_siteparam, pfile_param, pfile_soil, pfile_pdur)


    ########
    # 1) Parse pfile inputs:

    # heat flow (unimplemented)
    p_HEAT    = pfile_param["HEAT"]
    p_TopInfT = pfile_soil["TopInfT"]

    # radiation aspect parameters
    p_LAT   = pfile_siteparam["p_LAT"]
    p_GLMAX = pfile_param["GLMAX"]
    p_GLMIN = pfile_param["GLMIN"]
    p_ESLOPE= pfile_param["ESLOPE"]
    p_DSLOPE= pfile_param["DSLOPE"]
    # equivalent slope for radiation calculations
    p_L1, p_L2 = LWFBrook90Julia.SUN.EQUIVSLP(p_LAT, p_ESLOPE, pfile_param["ASPECT"])

    p_SNODEN = pfile_param["SNODEN"]
    p_MXRTLN = pfile_param["MXRTLN"]
    p_MXKPL  = pfile_param["MXKPL"]
    p_CS     = pfile_param["CS"]
    p_Z0S    = pfile_param["Z0S"]
    p_Z0G    = pfile_param["Z0G"]
    p_ZMINH  = pfile_param["ZMINH"]
    p_CZS    = pfile_param["CZS"]
    p_CZR    = pfile_param["CZR"]
    p_HS     = pfile_param["HS"]
    p_HR     = pfile_param["HR"]
    p_LPC    = pfile_param["LPC"]
    p_RTRAD  = pfile_param["RTRAD"]
    p_FXYLEM = pfile_param["FXYLEM"]
    p_WNDRAT = pfile_param["WNDRAT"]
    p_FETCH  = pfile_param["FETCH"]
    p_Z0W    = pfile_param["Z0W"]
    p_ZW     = pfile_param["ZW"]
    p_RSTEMP = pfile_param["RSTEMP"]
    p_LWIDTH = pfile_param["LWIDTH"]
    p_RHOTP  = pfile_param["RHOTP"]
    p_NN     = pfile_param["NN"]
    p_KSNVP  = pfile_param["KSNVP"]
    p_ALBSN  = pfile_param["ALBSN"]
    p_ALB    = pfile_param["ALB"]
    p_RSSA   = pfile_param["RSSA"]
    p_RSSB   = pfile_param["RSSB"]
    p_CCFAC  = pfile_param["CCFAC"]
    p_MELFAC = pfile_param["MELFAC"]
    p_LAIMLT = pfile_param["LAIMLT"]
    p_SAIMLT = pfile_param["SAIMLT"]
    p_GRDMLT = pfile_param["GRDMLT"]
    p_C1     = pfile_param["C1"]
    p_C2     = pfile_param["C2"]
    p_C3     = pfile_param["C3"]
    p_CR     = pfile_param["CR"]
    p_R5     = pfile_param["R5"]
    p_CVPD   = pfile_param["CVPD"]
    p_RM     = pfile_param["RM"]
    p_TL     = pfile_param["TL"]
    p_T1     = pfile_param["T1"]
    p_T2     = pfile_param["T2"]
    p_TH     = pfile_param["TH"]
    p_FSINTL = pfile_param["FSINTL"]
    p_FSINTS = pfile_param["FSINTS"]
    p_CINTSL = pfile_param["CINTSL"]
    p_CINTSS = pfile_param["CINTSS"]
    p_FRINTL = pfile_param["FRINTL"]
    p_FRINTS = pfile_param["FRINTS"]
    p_CINTRL = pfile_param["CINTRL"]
    p_CINTRS = pfile_param["CINTRS"]
    p_MAXLQF = pfile_param["MAXLQF"]
    p_QFPAR  = pfile_param["QFPAR"]
    p_QFFC   = pfile_param["QFFC"]
    p_IMPERV = pfile_param["IMPERV"]
    p_BYPAR  = pfile_param["BYPAR"]
    p_LENGTH = pfile_param["LENGTH"]
    p_DPSIMX = pfile_param["DPSIMX"]
    p_DRAIN  = pfile_param["DRAIN"]
    p_DTIMAX = pfile_param["DTIMAX"]
    p_DSWMAX = pfile_param["DSWMAX"]
    p_GSC    = pfile_param["GSC"]
    p_GSP    = pfile_param["GSP"]

    p_DURATN = pfile_pdur["DURATN"]

    # root density parameters
    p_inirdep    = pfile_param["inirdep"] # for RootGrowth in LWFBrook90Julia
    p_inirlen    = pfile_param["inirlen"] # for RootGrowth in LWFBrook90Julia
    p_rgroper    = pfile_param["rgroper"] # for RootGrowth in LWFBrook90Julia
    p_tini       = pfile_soil["tini"]     # for RootGrowth in LWFBrook90Julia
    p_frelden    = pfile_soil["frelden"]  # for RootGrowth in LWFBrook90Julia

    # unused p_HeatCapOld = pfile_soil["HeatCapOld"]

    # soil water parameters
    # unused u_aux_PSIMinit = pfile_soil["PSIM_init"]
    p_PSICR= pfile_param["PSICR"]
    p_THICK= pfile_soil["THICK"]
    p_STONEF=pfile_soil["STONEF"]
    p_THSAT =pfile_soil["PAR"][!,"θs"]
    if IMODEL == 0
        p_THETAF = pfile_soil["PAR"][!,"θf"]
        p_KF     = pfile_soil["PAR"][!,"kf"]
        p_PSIF   = pfile_soil["PAR"][!,"ψf"]
        p_BEXP   = pfile_soil["PAR"][!,"bexp"]
        p_WETINF = pfile_soil["PAR"][!,"wetinf"]
        p_Kθfc   = fill(NaN, NLAYER)
        p_Ksat   = fill(NaN, NLAYER)
        p_MvGα   = fill(NaN, NLAYER)
        p_MvGn   = fill(NaN, NLAYER)
        p_MvGl   = fill(NaN, NLAYER)
        p_θr     = fill(NaN, NLAYER)
    elseif IMODEL == 1
        p_THETAF = fill(NaN, NLAYER)
        p_KF     = fill(NaN, NLAYER)
        p_PSIF   = fill(NaN, NLAYER)
        p_BEXP   = fill(NaN, NLAYER)
        p_WETINF = fill(NaN, NLAYER)
        p_Kθfc = pfile_soil["PAR"][!,"K(θ_fc)"]
        p_Ksat = pfile_soil["PAR"][!,"Ksat"]
        p_MvGα = pfile_soil["PAR"][!,"α"]
        p_MvGn = pfile_soil["PAR"][!,"n"]
        p_MvGl = pfile_soil["PAR"][!,"tort"]
        p_θr   = pfile_soil["PAR"][!,"θr"]
    else
        error("Unsupported IMODEL: $IMODEL")
    end

    # infiltration parameters INFPAR()
    p_ILAYER = pfile_param["ILAYER"] # number of layers over which infiltration is distributed
    p_INFEXP = pfile_param["INFEXP"]
    p_INFRAC = LWFBrook90Julia.WAT.INFPAR(p_INFEXP, p_ILAYER, p_THICK, NLAYER)


    # soil water parameters
    (p_PSIG,             # gravity potential (kPa),
    p_SWATMX,           # maximum water storage for layer (mm),
    p_WETF,             # wetness at field capacity (dimensionless),
    #p_WETc, # wetness at field capacity (dimensionless) # TODO(bernhard): was unused with Brook90R, check if used in LWFBrook90R
    p_CHM,              # Clapp and Hornberger m (kPa),
    p_CHN,              # Clapp and Hornberger n,
    p_KSAT,             # saturated hydraulic conductivity (mm/d)
    p_PSIF,
    p_THETAF
    # TODO(bernhard): clean up use of p_PSIF and p_THETAF (depending on IMODEL they are
    #                 initialized as NaN. Are they overwritten by SOILPAR?)
    ) = LWFBrook90Julia.KPT.SOILPAR(
            LWFBrook90Julia.CONSTANTS.p_RHOWG,
            p_THICK, # layer thicknesses (mm)
            p_THETAF,# θ at field capacity (-)
            p_THSAT, # θ at saturation == matrix porosity (-)

            p_STONEF,# stone volume fraction, unitless
            p_BEXP,  # exponent for psi-theta relation
            p_KF,    # hydraulic conductivity at field capacity (mm/d)
            p_PSIF,  # ψ at field capacity (kPa)
            p_WETINF,# wetness at dry end of near-saturation range

            p_Kθfc,
            p_PSICR, # minimum plant leaf water potential (MPa)
            p_Ksat, p_MvGl, p_MvGn, p_MvGα, p_θr,

            NLAYER, IMODEL)

    # p_PsiCrit is the ψ value that corresponds to the constant, critical θ value p_ThCrit
    # Note that p_PSICR is different!
    if (IMODEL == 0)
        p_PsiCrit = LWFBrook90Julia.KPT.FPSIMF_CH.(LWFBrook90Julia.CONSTANTS.p_ThCrit./p_THSAT,
                                                   p_PSIF, p_BEXP, p_WET∞, p_WETF, p_CHM, p_CHN)
    elseif (IMODEL == 1)
        p_PsiCrit = LWFBrook90Julia.KPT.FPSIM_MvG.(LWFBrook90Julia.CONSTANTS.p_ThCrit./(p_THSAT .- p_θr),
                                                   p_MvGα, p_MvGn)
    else
        error("Unknown IMODEL!")
    end
    # TODO(bernhard): treat following note:
    # NOTE(bernhard) the difference between p_PSICR and p_PsiCrit:
    # p_PSICR (Brook90): PSICR (Canopy parameter) - minimum plant leaf water
    #    potential, MPa. PSICR is the critical leaf water potential at which stomates
    #    close. BROOK90 assumes that transpiration is limited by potential
    #    transpiration (PTRAN) until water uptake at a plant water potential of PSICR
    #    is less than PTRAN. PSICR can be considered as the water potential at the
    #    turgor-loss point. PSICR varies from -1.5 to -3.0 MPa for most species and is
    #    quite species dependent (Hinckley et al. 1978). This parameter is best
    #    selected from knowledge of the water potential - diffusion resistance relation
    #    for the species involved. [see KPT-SOILPAR] [see EVP-TBYLAYER]
    # p_PsiCrit (LWFBrook90):
    #    Definition: minimum soil matric potential to allow water supply for
    #    evapotranspiration , Hammel, 2001 (p_PsiCrit = f(ThCrit) = FPSIM(ThCrit))
    #    use: if (PsiM < PsiCrit) TRANI = 0 # no transpiration
    #    use: if (PsiM < PsiCrit) SLVP = 0  # no soil evaporation # TODO(bernhard): this seems incorrect

    # source area parameters SRFPAR()
    p_QLAYER = # number of soil layers for SRFL
        pfile_param["QLAYER"]
    p_SWATQX = # maximum water storage for layers 1 through QLAYER (mm)
        sum(p_SWATMX[1:pfile_param["QLAYER"]])
    p_SWATQF = # water storage at field capacity for layers 1 through QLAYER (mm)
        sum(p_THETAF[1:pfile_param["QLAYER"]] .* p_THICK[1:pfile_param["QLAYER"]] .* (1 .- p_STONEF[1:pfile_param["QLAYER"]]))


    ########
    # 2) Define parameters for differential equation:
    # 2a) Constant parameters

    # p_cst_1 for both RHS and CallBack in DiffEq.jl
    p_cst_1 = (constant_dt_solver, NLAYER, IMODEL, compute_intermediate_quantities, Reset,
        p_SWATMX, p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_PSIG, p_KF,
        p_THSAT, p_θr, p_MvGα, p_MvGn, p_MvGl, p_Ksat,

        # FOR MSBITERATE:
        p_QLAYER, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC, p_IMPERV,
        p_LENGTH, p_DSLOPE, LWFBrook90Julia.CONSTANTS.p_RHOWG, p_DPSIMX,
        p_KSAT, p_DRAIN, p_DTIMAX, p_INFRAC, p_DSWMAX,
        p_GSC, p_GSP, p_THICK, p_STONEF,

        # FOR UNIMPLEMENTED HEAT FLOW:
        p_HEAT, p_TopInfT,

        p_BYPAR)

    # p_cst_2 only for CallBack in DiffEq.jl
    p_cst_2 = (p_LAT, p_ESLOPE, p_L1, p_L2,
        p_SNODEN, p_MXRTLN, p_MXKPL, p_CS,
        p_Z0S, p_Z0G,
        p_ZMINH, p_CZS, p_CZR, p_HS, p_HR, p_LPC,
        p_RTRAD, p_FXYLEM,
        p_WNDRAT, p_FETCH, p_Z0W, p_ZW,
        p_RSTEMP,
        LWFBrook90Julia.CONSTANTS.p_CVICE,
        p_LWIDTH, p_RHOTP, p_NN, p_KSNVP,
        p_ALBSN, p_ALB,
        p_RSSA, p_RSSB,
        p_CCFAC, p_MELFAC, p_LAIMLT, p_SAIMLT,

        p_inirdep, p_inirlen, p_rgroper, p_tini, p_frelden,

        LWFBrook90Julia.CONSTANTS.p_WTOMJ, p_C1, p_C2, p_C3, p_CR,
        p_GLMIN, p_GLMAX, p_R5, p_CVPD, p_RM, p_TL, p_T1, p_T2, p_TH,
        p_PSICR, NOOUTF, p_PsiCrit,

        # for MSBPREINT:
        p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS,
        p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS,
        p_DURATN, p_MAXLQF, p_GRDMLT)

    p_cst = (p_cst_1, p_cst_2)

    # 2b) Time varying parameters (e.g. meteorological forcings)
    p_NPINT = pfile_siteparam["p_NPINT"]
    if (isequal(p_NPINT, 1))
        # one precipitation interval per day:
        #   p_DTP = 1.0 day, i.e equal to simulatino interval p_DT = 1.0 day
        #   use precipitation rate (mm/day) = p_PRECIN/1 = p_PRECIN/p_DTP
        p_DTP = 1 # time step for precipitation interval (for PRECIN and MESFL), may be <= 1 d
    else
        # multiple precipitation intervals per day: (not implemented)
        p_DTP = Nothing # time step for precipitation interval (for PRECIN and MESFL), may be <= 1 d
        error("Case with multiple precipitation intervals (using PRECDAT) is not implemented.")
        # This would result in use of e.g. hourly PRECIN and hourly MESFL
    end

    p_DOY_inclRef = (t) -> LWFBrook90Julia.p_DOY(t, pfile_meteo["input_reference_date"])
    p_MONTHN_inclRef = (t) -> LWFBrook90Julia.p_MONTHN(t, pfile_meteo["input_reference_date"])
    p_fT = (p_DOY_inclRef,
            p_MONTHN_inclRef,
            pfile_meteo["p_GLOBRAD"],
            pfile_meteo["p_TMAX"],
            pfile_meteo["p_TMIN"],
            pfile_meteo["p_VAPPRES"],
            pfile_meteo["p_WIND"],
            pfile_meteo["p_PRECIN"],
            p_DTP, p_NPINT,
            pfile_meteo["p_MESFL"],
            pfile_meteo["p_DENSEF"],
            pfile_meteo["p_HEIGHT"],
            pfile_meteo["p_LAI"],
            pfile_meteo["p_SAI"],
            pfile_meteo["p_AGE"])

    # 2c) Time varying "parameters" (depending on state variables)
    #     These need to be exchanged between CallBack and RHS in DiffEq.jl which is why they
    #     can temporarily be saved in the parameter vector to avoid computing them twice

    # Initialize placeholder for parameters that depend on solution and are computed
    p_fu = [NaN, NaN, fill(NaN, NLAYER), NaN] # Use array instead of tuple to be able to mutate

    # 3) Return different types of parameters as a single object
    return (p_cst, p_fT, p_fu)
end