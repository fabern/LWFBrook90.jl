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
            du[7+NLAYER+27]= 0 # was computed in callback         # was timeseries_pslvp

            # TODO(bernhard): use SavingCallback() for all quantities that have du=0
            #                 only keep du=... for quantities for which we compute cumulative sums
        end

        return
    end
#     f_func = f_LWFBrook90R
#     return f_func
# end
