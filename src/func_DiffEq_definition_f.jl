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
        p_DOY, p_MONTHN, p_SOLRAD, p_TMAX, p_TMIN, p_EA, p_UW, p_PREC,
            p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_AGE, p_RELDEN,
            p_δ18O_PREC, p_δ2H_PREC, ref_date = p[2]

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
        (δ18O_SLFL, δ2H_SLFL,
            p_fu_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP,
            p_fu_STHR, aux_du_RSNO, aux_du_SNVP,
            aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP, u_SNOW_old) = p[3][1]
        aux_du_TRANI = p[3][2]

        # Pre-allocated caches to save memory allocations
        (u_aux_WETNES,u_aux_PSIM,u_aux_PSITI,u_aux_θ,u_aux_θ_tminus1,p_fu_KK,
            aux_du_DSFLI,aux_du_VRFLI,aux_du_VRFLI_1st_approx,aux_du_INFLI,aux_du_BYFLI, du_NTFLI,
            p_fu_BYFRAC) = p[4][1]

        ##################
        # Parse states
        states = NamedTuple{p[1][4].u0_field_names}(u.x)
        rates  = NamedTuple{p[1][4].u0_field_names}(du.x)

        @assert p[1][4].u0_variable_names == (d18O = 2, d2H = 3)
        u_GWAT      = states[:GWAT  ][1,1]
        # u_INTS      = states[:INTS  ][1,1]
        # u_INTR      = states[:INTR  ][1,1]
        # u_SNOW      = states[:SNOW  ][1,1]
        # u_CC        = states[:CC    ][1,1]
        # u_SNOWLQ    = states[:SNOWLQ][1,1]
        u_SWATI     = states[:SWATI ][:,1]

        # simulate_isotopes = p[1][4].simulate_isotopes
        # if simulate_isotopes
        #     u_δ18O_GWAT   = states[:GWAT ][1,p[1][4].col_idx_d18O]
        #     u_δ2H_GWAT    = states[:GWAT ][1,p[1][4].col_idx_d2H ]
        #     # u_δ18O_INTS = states[:INTS ][1,p[1][4].col_idx_d18O]
        #     # u_δ2H_INTS  = states[:INTS ][1,p[1][4].col_idx_d2H ]
        #     # u_δ18O_INTR = states[:INTR ][1,p[1][4].col_idx_d18O]
        #     # u_δ2H_INTR  = states[:INTR ][1,p[1][4].col_idx_d2H ]
        #     # u_δ18O_SNOW = states[:SNOW ][1,p[1][4].col_idx_d18O]
        #     # u_δ2H_SNOW  = states[:SNOW ][1,p[1][4].col_idx_d2H ]
        #     # u_δ18O_RWU  = states[:RWU  ][1,p[1][4].col_idx_d18O]
        #     # u_δ2H_RWU   = states[:RWU  ][1,p[1][4].col_idx_d2H ]
        #     # u_δ18O_XYL  = states[:XYLEM][1,p[1][4].col_idx_d18O]
        #     # u_δ2H_XYL   = states[:XYLEM][1,p[1][4].col_idx_d2H ]
        #     u_δ18O_SWATI  = states[:SWATI][:,p[1][4].col_idx_d18O]
        #     u_δ2H_SWATI   = states[:SWATI][:,p[1][4].col_idx_d2H ]
        # end

        LWFBrook90.KPT.SWCHEK!(u_SWATI, p_soil.p_SWATMAX, t)

        u_aux_θ_tminus1 .= u_aux_θ #TODO(bernhard): this does not seem to be correctly updated
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
        # If we compute scalar transport as ODEs, the du's need to be computed here in the f-function
        # However, if we compute scalar transport as separate step from the flow equation, we can
        #     compute them separately in a callback function that is called after each time step
        # if simulate_isotopes
        #     EffectiveDiffusivity_18O = 0  # TODO(bernhard): compute correct values using eq 33 Zhou-2021-Environ_Model_Softw.pdf
        #     EffectiveDiffusivity_2H  = 0  # TODO(bernhard): compute correct values using eq 33 Zhou-2021-Environ_Model_Softw.pdf
        #     δ18O_INFLI = δ18O_SLFL
        #     δ2H_INFLI  = δ2H_SLFL

        #     du_δ18O_GWAT, du_δ2H_GWAT, du_δ18O_SWATI, du_δ2H_SWATI =
        #         compute_isotope_du_GWAT_SWATI(
        #             # for GWAT:
        #             u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT,
        #             # for SWATI:
        #             du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI,  # (non-fractionating)
        #             aux_du_SLVP, p_fu_TADTM, p_EA(t), p_δ2H_PREC(t), p_δ18O_PREC(t), u_aux_WETNES, # (fractionating)
        #             u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, 0, 0) #EffectiveDiffusivity_18O, EffectiveDiffusivity_2H)
        # end

        ##################
        # Update solution:
        # u is a state vector with u[1] = S relative saturation (-)
        # Update GWAT:
        rates[:GWAT ][1,1] = aux_du_VRFLI[NLAYER] - du_GWFL - du_SEEP

        # Do not modify INTS, INTR, SNOW, CC, SNOWLQ
        # as they are separately modified by the callback.
        # Therfore set their rate of change to 0
        rates[:INTS  ][1,1] = 0
        rates[:INTR  ][1,1] = 0
        rates[:SNOW  ][1,1] = 0
        rates[:CC    ][1,1] = 0
        rates[:SNOWLQ][1,1] = 0

        # Update SWATI for each layer:
        rates[:SWATI ][:,1] .= du_NTFLI # 0.000002 seconds (3 allocations: 224 bytes)

        # save intermediate results from flow calculation into a cache
        # for efficient use in transport calculation (in callbacks)
        p[3][3]      .=  [du_NTFLI  aux_du_VRFLI  aux_du_DSFLI  aux_du_INFLI  u_aux_WETNES]
        p[3][4][1:2] .=  [du_GWFL, du_SEEP]

        ##########################################
        # Accumulate flows to compute daily sums
        # Note that below state variables serve only as accumulator but do not affect
        # the evolution of the system.
        if compute_intermediate_quantities
            # @assert length(p[1][4].row_idx_accum) == 25

            # 1) Either set daily sum if rate is constant throughout precipitation interval: p_DTP*(...)
            # 2) or then set daily sum to zero and use ODE to accumulate flow.
            rates[:accum][1]
            rates[:accum][ 1] = 0 # cum_d_prec ][:,was computed in callback
            rates[:accum][ 2] = 0 # cum_d_rfal ][:,was computed in callback
            rates[:accum][ 3] = 0 # cum_d_sfal ][:,was computed in callback
            rates[:accum][ 4] = 0 # cum_d_rint ][:,was computed in callback
            rates[:accum][ 5] = 0 # cum_d_sint ][:,was computed in callback
            rates[:accum][ 6] = 0 # cum_d_rsno ][:,was computed in callback
            rates[:accum][ 7] = 0 # cum_d_rnet ][:,was computed in callback
            rates[:accum][ 8] = 0 # cum_d_smlt ][:,was computed in callback

            rates[:accum][ 9] = 0 # cum_d_evap ][:, was computed in callback
            rates[:accum][10] = 0 # cum_d_tran ][:, was computed in callback
            rates[:accum][11] = 0 # cum_d_irvp ][:, was computed in callback
            rates[:accum][12] = 0 # cum_d_isvp ][:, was computed in callback
            rates[:accum][13] = 0 # cum_d_slvp ][:, was computed in callback
            rates[:accum][14] = 0 # cum_d_snvp ][:, was computed in callback
            rates[:accum][15] = 0 # cum_d_pint ][:, was computed in callback
            rates[:accum][16] = 0 # cum_d_ptran ][:,was computed in callback
            rates[:accum][17] = 0 # cum_d_pslvp ][:,was computed in callback

            rates[:accum][18] = p_fu_SRFL +
                              sum(aux_du_BYFLI) +
                              sum(aux_du_DSFLI) +
                              du_GWFL # SRFLD + BYFLD + DSFLD + GWFLD # flow
            rates[:accum][19] = du_SEEP                                 # seep
            rates[:accum][20] = p_fu_SRFL                               # srfl
            rates[:accum][21] = p_fu_SLFL                               # slfl
            rates[:accum][22] = sum(aux_du_BYFLI)                       # byfl
            rates[:accum][23] = sum(aux_du_DSFLI)                       # dsfl
            rates[:accum][24] = du_GWFL                                 # gwfl
            rates[:accum][25] = aux_du_VRFLI[NLAYER]                    # vrfln

            # rates[:accum][26] = 0 # cum_d_rthr ][:,was computed in callback
            # rates[:accum][27] = 0 # cum_d_sthr ][:,was computed in callback
            rates[:accum][28] = 0 # new_SWAT ][:,is computed in callback
            rates[:accum][29] = 0 # new_totalWATER ][:,is computed in callback
            rates[:accum][30] = 0 # BALERD_SWAT ][:,is computed in callback
            rates[:accum][31] = 0 # BALERD_total ][:,is computed in callback

            # TODO(bernhard): use SavingCallback() for all quantities that have du=0
            #                 only keep du=... for quantities for which we compute cumulative sums
        end

        # @infiltrate any(abs.(du_δ2H_SWATI) .> 10)
        # @infiltrate any(abs.(du_δ18O_SWATI) .> 10)
        # @infiltrate any(isnan.(du_δ2H_SWATI))
        # @infiltrate any(isnan.(du_δ18O_SWATI))
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
