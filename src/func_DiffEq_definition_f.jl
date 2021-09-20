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
        (p_fu_RNET, aux_du_SMLT, aux_du_SLVP) = p[3][1]
        aux_du_TRANI = p[3][2]


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
        u_SWATI     = u[idx_u_vector_amounts] # 0.000002 seconds (1 allocation: 144 bytes)

        LWFBrook90.KPT.SWCHEK!(u_SWATI, p_soil.p_SWATMAX, t)

        (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
            LWFBrook90.KPT.derive_auxiliary_SOILVAR(u_SWATI, p_soil) # 0.000007 seconds (7 allocations: 1008 bytes)

        ##################
        # Update soil limited boundary flows during iteration loop
        if Reset == 1
            error("The case with updated flows (Reset==1) is not implemented here. A test is implemented in branch 005.")
        end

        ##################
        # Compute fluxes

        # Bypass fraction of infiltration to each layer
        p_fu_BYFRAC = LWFBrook90.WAT.BYFLFR(
                      NLAYER, p_BYPAR, p_QFPAR, p_QFFC, u_aux_WETNES, p_soil) # 0.000002 seconds (1 allocation: 144 bytes)

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
                    p_DRAIN, p_DTP, t, p_DTIMAX,
                    # for INFLOW:
                    p_INFRAC, p_fu_BYFRAC, aux_du_TRANI, aux_du_SLVP,
                    # for FDPSIDW:
                    u_aux_WETNES,
                    # for ITER:
                    p_DSWMAX, u_aux_θ,
                    # for GWATER:
                    u_GWAT, p_GSC, p_GSP) # 0.000011 seconds (15 allocations: 2.109 KiB)

        # Transport flow (heat, solutes, isotopes, ...)
        # TODO(bernhard): see initial prototype code... (script3)
        # du_δ2H_GWAT, du_δ2H_SWATI, du_δ2H_INTS, du_δ2H_INTR, du_δ2H_SNOW =
        #     TODO(bernhard)
        du_δ18O_GWAT, du_δ18O_SWATI, du_δ18O_INTS, du_δ18O_INTR, du_δ18O_SNOW =
            compute_isotope_du( # check chart "../docs/src/assets/b90flow.gif"
                # for GWAT:
                du_GWFL, du_SEEP, # aux_du_VRFLI, δ18O_GWAT (out), δ18O_SWATI (in)
                # for SWATI:
                aux_du_VRFLI, δ18O_SWATI, Diffusivity_18O,
                aux_du_TRANI, aux_du_SLVP, # (non-fractionated)
                aux_du_INFLI, δ18O_INFLI,
                # TODO(bernhard): acutally instead of INFLI we need to use SMLT(δ_SNOW) and RNET(δ_PREC)
                # TODO(benrhard): or alternatively compute δ_INFLI = (p_fu_RNET*δ18O_Precip + aux_du_SMLT*δ18O_SNOW)/(p_fu_RNET + aux_du_SMLT)

                # for INTS, INTR, SNOW: (is done in daily time steps (i.e. in the callback function cb()))
                ###
                # TODO(bernhard): also compute δ of EVAP, SEEP and FLOW
                )

        ##################
        # Update solution:
        # u is a state vector with u[1] = S relative saturation (-)
        # Update GWAT:
        du[1] = aux_du_VRFLI[NLAYER] - du_GWFL - du_SEEP

        # Do not modify INTS, INTR, SNOW, CC, SNOWLQ
        # as they are separately modified by the callback.
        # Therfore set their rate of change to 0
        du[2] = 0
        du[3] = 0
        du[4] = 0
        du[5] = 0
        du[6] = 0

        # Update SWATI for each layer:
        du[idx_u_vector_amounts] .= du_NTFLI[:] # 0.000002 seconds (3 allocations: 224 bytes)


        ##########################################
        # Accumulate flows to compute daily sums
        # Note that below state variables serve only as accumulator but do not affect
        # the evolution of the system.
        if compute_intermediate_quantities
            # @assert length(idx_u_vector_accumulators) == 25

            # 1) Either set daily sum if rate is constant throughout precipitation interval: p_DTP*(...)
            # 2) or then set daily sum to zero and use ODE to accumulate flow.
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


        # NOTE(bernharf): To replicate behavior of LWFBrook90R with Reset == True, the following would be needed (but is not planned):
        #                 - with Reset == TRUE, certain quantities are updated after each evaluation of the rhs `f`
        #                 - see the LWFBrook90R code between:
        #                   https://github.com/pschmidtwalter/LWFBrook90R/blob/6f23dc1f6be9e1723b8df5b188804da5acc92e0f/src/md_brook90.f95#L578
        #                   and
        #                   https://github.com/pschmidtwalter/LWFBrook90R/blob/6f23dc1f6be9e1723b8df5b188804da5acc92e0f/src/md_brook90.f95#L703

    end
#     f_func = f_LWFBrook90R
#     return f_func
# end
