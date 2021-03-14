"""
    define_LWFB90_f()

Generate function f (right-hand-side of ODEs) needed for ODE() problem in DiffEq.jl package.
"""
# function define_LWFB90_f()
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
        p_soil,

        # FOR MSBITERATE:
        p_QLAYER, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC, p_IMPERV,
        p_LENGTH, p_DSLOPE, p_RHOWG, p_DPSIMX, #TODO(bernhard) p_RHOWG is a global constant
        p_DRAIN, p_DTIMAX, p_INFRAC, p_DSWMAX,
        p_GSC, p_GSP,

        p_BYPAR) = p[1][1]
        # unused are the constant parameters saved in: = p[1][2]

        ## B) time dependent parameters
        (p_DOY, p_MONTHN, p_SOLRAD, p_TMAX, p_TMIN, p_EA, p_UW, p_PREC, p_DTP, p_NPINT, p_MESFL,
        _, _, _, _, _) = p[2]

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
            LWFBrook90.KPT.derive_auxiliary_SOILVAR(u_SWATI, p_soil)

        ##################
        # Update soil limited boundary flows during iteration loop
        if Reset == 1
            error("The case with updated flows (Reset==1) is not implemented here. A test is implemented in branch 005.")
        end

        ##################
        # Compute fluxes

        # Bypass fraction of infiltration to each layer
        p_fu_BYFRAC = LWFBrook90.WAT.BYFLFR(
                      NLAYER, p_BYPAR, p_QFPAR, p_QFFC, u_aux_WETNES, p_soil)

        # Water movement through soil
        (p_fu_SRFL, p_fu_SLFL, aux_du_DSFLI, aux_du_VRFLI, DTI, aux_du_INFLI, aux_du_BYFLI,
        du_NTFLI, du_GWFL, du_SEEP) =
            MSBITERATE(IMODEL, NLAYER, p_QLAYER, p_soil,
                    # for SRFLFR:
                    u_SWATI, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC,
                    #
                    p_IMPERV, p_fu_RNET, aux_du_SMLT,
                    p_LENGTH, p_DSLOPE,
                    # for DSLOP:
                    p_RHOWG, u_aux_PSIM, p_fu_KK,
                    #
                    u_aux_PSITI, p_DPSIMX,
                    #
                    p_DRAIN, DTRI, p_DTIMAX,
                    # for INFLOW:
                    p_INFRAC, p_fu_BYFRAC, aux_du_TRANI, aux_du_SLVP,
                    # for FDPSIDW:
                    u_aux_WETNES,
                    # for ITER:
                    p_DSWMAX, u_aux_θ,
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
            du[7+NLAYER+24] = aux_du_VRFLI[NLAYER]                    # vrfln

            # du[7+NLAYER+25]= 0 # cum_d_rthr, was computed in callback
            # du[7+NLAYER+26]= 0 # cum_d_sthr, was computed in callback
            # du[7+NLAYER+27]= 0 # mesfl, was computed in callback
            # balerd[IDAY]=BALERD

            # TODO(bernhard): use SavingCallback() for all quantities that have du=0
            #                 only keep du=... for quantities for which we compute cumulative sums
        end

        return
    end
#     f_func = f_LWFBrook90R
#     return f_func
# end
