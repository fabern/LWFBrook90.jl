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
        p_soil = p[1][1]
        (NLAYER, FLAG_MualVanGen, compute_intermediate_quantities, Reset,
        p_DTP, p_NPINT,

        # FOR MSBITERATE:
        p_QLAYER, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC, p_IMPERV,
        p_LENGTH_SLOPE, p_DSLOPE, p_RHOWG, p_DPSIMAX, #TODO(bernhard) p_RHOWG is a global constant
        p_DRAIN, p_DTIMAX, p_INFRAC, p_DSWMAX,
        p_GSC, p_GSP,

        p_BYPAR) = p[1][2]
        # unused are the constant parameters saved in: = p[1][3]

        ## B) time dependent parameters
        (p_DOY, p_MONTHN, p_SOLRAD, p_TMAX, p_TMIN, p_EA, p_UW, p_PREC,
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
        idx_u_vector_amounts       = p[1][4][4] # 7:(6+NLAYER)
        idx_u_vector_accumulators  = p[1][4][5] # 6+0+0+1*NLAYER .+ (1:25)
        # idx_u_scalar_amounts       = p[1][4][6] # 1:6

        u_GWAT     = u[1]
        #u_INTS     = u[2]
        #u_INTR     = u[3]
        #u_SNOW     = u[4]
        #u_CC       = u[5]
        #u_SNOWLQ   = u[6]
        u_SWATI     = u[idx_u_vector_amounts]

        LWFBrook90.KPT.SWCHEK!(u_SWATI, p_soil.p_SWATMAX, t)

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
            MSBITERATE(FLAG_MualVanGen, NLAYER, p_QLAYER, p_soil,
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
                    p_DRAIN, DTRI, p_DTIMAX,
                    # for INFLOW:
                    p_INFRAC, p_fu_BYFRAC, aux_du_TRANI, aux_du_SLVP,
                    # for FDPSIDW:
                    u_aux_WETNES,
                    # for ITER:
                    p_DSWMAX, u_aux_θ,
                    # for GWATER:
                    u_GWAT, p_GSC, p_GSP)

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
        du[idx_u_vector_amounts] .= du_NTFLI[:]

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
            # @assert length(idx_u_vector_accumulators) == 25

            # 1) Either set daily sum if rate is constant throughout precipitation interval: p_DTP*(...)
            # 2) or then set daily sum to zero and use ODE to accumulate flow.
            6+0+0+1*NLAYER .+ (1:25)
            du[idx_u_vector_accumulators[ 1]] = 0 # cum_d_prec, was computed in callback
            du[idx_u_vector_accumulators[ 2]] = 0 # cum_d_rfal, was computed in callback
            du[idx_u_vector_accumulators[ 3]] = 0 # cum_d_sfal, was computed in callback
            du[idx_u_vector_accumulators[ 4]] = 0 # cum_d_rint, was computed in callback
            du[idx_u_vector_accumulators[ 5]] = 0 # cum_d_sint, was computed in callback
            du[idx_u_vector_accumulators[ 6]] = 0 # cum_d_rsno, was computed in callback
            du[idx_u_vector_accumulators[ 7]] = 0 # cum_d_rnet, was computed in callback
            du[idx_u_vector_accumulators[ 8]] = 0 # cum_d_smlt, was computed in callback

            du[idx_u_vector_accumulators[ 9]] = 0 # cum_d_evap,  was computed in callback
            du[idx_u_vector_accumulators[10]] = 0 # cum_d_tran,  was computed in callback
            du[idx_u_vector_accumulators[11]] = 0 # cum_d_irvp,  was computed in callback
            du[idx_u_vector_accumulators[12]] = 0 # cum_d_isvp,  was computed in callback
            du[idx_u_vector_accumulators[13]] = 0 # cum_d_slvp,  was computed in callback
            du[idx_u_vector_accumulators[14]] = 0 # cum_d_snvp,  was computed in callback
            du[idx_u_vector_accumulators[15]] = 0 # cum_d_pint,  was computed in callback
            du[idx_u_vector_accumulators[16]] = 0 # cum_d_ptran, was computed in callback
            du[idx_u_vector_accumulators[17]] = 0 # cum_d_pslvp, was computed in callback

            du[idx_u_vector_accumulators[18]] = p_fu_SRFL +
                              sum(aux_du_BYFLI) +
                              sum(aux_du_DSFLI) +
                              du_GWFL # SRFLD + BYFLD + DSFLD + GWFLD # flow
            du[idx_u_vector_accumulators[19]] = du_SEEP                                 # seep
            du[idx_u_vector_accumulators[20]] = p_fu_SRFL                               # srfl
            du[idx_u_vector_accumulators[21]] = p_fu_SLFL                               # slfl
            du[idx_u_vector_accumulators[22]] = sum(aux_du_BYFLI)                       # byfl
            du[idx_u_vector_accumulators[23]] = sum(aux_du_DSFLI)                       # dsfl
            du[idx_u_vector_accumulators[24]] = du_GWFL                                 # gwfl
            du[idx_u_vector_accumulators[25]] = aux_du_VRFLI[NLAYER]                    # vrfln

            # du[idx_u_vector_accumulators[26]]= 0 # cum_d_rthr, was computed in callback
            # du[idx_u_vector_accumulators[27]]= 0 # cum_d_sthr, was computed in callback
            # balerd[IDAY]=BALERD

            # TODO(bernhard): use SavingCallback() for all quantities that have du=0
            #                 only keep du=... for quantities for which we compute cumulative sums
        end

        return
    end
#     f_func = f_LWFBrook90R
#     return f_func
# end
