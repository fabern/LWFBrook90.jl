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
        @unpack p_soil = p;
        @unpack NLAYER, FLAG_MualVanGen, compute_intermediate_quantities, Reset, p_DTP,
            # FOR MSBITERATE:
            p_QLAYER, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC, p_IMPERV,
            p_LENGTH_SLOPE, p_DSLOPE, p_RHOWG, p_DPSIMAX, #TODO(bernhard) p_RHOWG is a global constant
            p_DRAIN, p_DTIMAX, p_INFRAC, p_DSWMAX, p_GSC, p_GSP, p_BYPAR = p;

        ## B) time dependent parameters
        @unpack p_DOY, p_MONTHN, p_GLOBRAD, p_TMAX, p_TMIN, p_VAPPRES, p_WIND,  p_PREC,
            p_DENSEF, p_HEIGHT, p_LAI, p_SAI, p_fT_RELDEN,
            p_δ18O_PREC, p_δ2H_PREC, REFERENCE_DATE = p;

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
        # @unpack p_fu_RNET, aux_du_SMLT, aux_du_SLVP, aux_du_TRANI = p;
        @unpack p_fu_δ18O_SLFL, p_fu_δ2H_SLFL,
            p_fT_TADTM, p_fu_RNET, aux_du_SMLT, aux_du_SLVP, p_fu_STHR,
            aux_du_RSNO, aux_du_SNVP, aux_du_SINT, aux_du_ISVP, aux_du_RINT, aux_du_IRVP,
            u_SNOW_old, aux_du_TRANI = p

        # Pre-allocated caches to save memory allocations
        @unpack du_GWFL, du_SEEP, du_NTFLI, aux_du_VRFLI, aux_du_DSFLI, aux_du_INFLI, u_aux_WETNES = p;
        @unpack u_aux_WETNES,u_aux_PSIM,u_aux_PSITI,u_aux_θ,u_aux_θ_tminus1,p_fu_KK,
            aux_du_VRFLI_1st_approx, aux_du_BYFLI, p_fu_BYFRAC = p;

        ##################
        # Parse states
        u_GWAT      = u.GWAT.mm
        # u_INTS      = u.INTS.mm
        # u_INTR      = u.INTR.mm
        # u_SNOW      = u.SNOW.mm
        # u_CC        = u.CC.MJm2
        # u_SNOWLQ    = u.SNOWLQ.mm
        u_SWATI     = u.SWATI.mm

        LWFBrook90.KPT.SWCHEK!(u_SWATI, p_soil.p_SWATMAX, t)

        u_aux_θ_tminus1 .= u_aux_θ #TODO(bernhard): this does not seem to be correctly updated
        (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
            LWFBrook90.KPT.derive_auxiliary_SOILVAR(u_SWATI, p_soil) # 0.000007 seconds (7 allocations: 1008 bytes)
        # LWFBrook90.KPT.derive_auxiliary_SOILVAR!(u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK, u_SWATI,  p_soil)

        ##################
        # Update soil limited boundary flows during iteration loop
        if Reset == 1
            error("The case with updated flows (Reset==1) is not implemented here. A test is implemented in branch 005.")
        end

        ##################
        # Compute fluxes

        # Bypass fraction of infiltration to each layer
        p_fu_BYFRAC = fill(NaN, p_soil.NLAYER)
        if (isone(p_BYPAR))
            LWFBrook90.WAT.BYFLFR!(p_fu_BYFRAC, p_QFPAR, p_QFFC, u_aux_WETNES, p_soil) # 0.000002 seconds (1 allocation: 144 bytes)
            # generate bypass flow to avoid saturation
            for i = 1:p_soil.NLAYER
                if (u_aux_WETNES[i] > 0.99)
                    p_fu_BYFRAC[i] = 1
                end
            end
        else
            p_fu_BYFRAC .= 0 # alternatively
        end


        # Water movement through soil (modify DTI time step if needed)
        p_fu_SRFL, p_fu_SLFL, aux_du_DSFLI[:], aux_du_VRFLI[:], aux_du_INFLI[:], aux_du_BYFLI[:], du_NTFLI[:], DTI =
            MSBITERATE(FLAG_MualVanGen, NLAYER, p_QLAYER, p_soil,
                    # for SRFLFR:
                    u_SWATI, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC,
                    #
                    p_IMPERV, p_fu_RNET[1], aux_du_SMLT[1],
                    p_LENGTH_SLOPE, p_DSLOPE,
                    # for DSLOP:
                    p_RHOWG, u_aux_PSIM, p_fu_KK,
                    #
                    u_aux_PSITI, p_DPSIMAX,
                    #
                    p_DRAIN, p_DTP, t, p_DTIMAX,
                    # for INFLOW:
                    p_INFRAC, p_fu_BYFRAC, aux_du_TRANI, aux_du_SLVP[1],
                    # for FDPSIDW:
                    u_aux_WETNES,
                    # for ITER:
                    p_DSWMAX, u_aux_θ) # 0.000011 seconds (15 allocations: 2.109 KiB)

        # groundwater flow and seepage loss
        du_GWFL[1], du_SEEP[1] = LWFBrook90.WAT.GWATER(u_GWAT, p_GSC, p_GSP, aux_du_VRFLI[p_soil.NLAYER]) # [1] as they are vectors and this doesn't allocate

        # Transport flow (heat, solutes, isotopes, ...)
        # If we compute scalar transport as ODEs, the du's need to be computed here in the f-function
        # However, if we compute scalar transport as separate step from the flow equation, we can
        #     compute them separately in a callback function that is called after each time step
        # if simulate_isotopes
        #     EffectiveDiffusivity_18O = 0  # TODO(bernhard): compute correct values using eq 33 Zhou-2021-Environ_Model_Softw.pdf
        #     EffectiveDiffusivity_2H  = 0  # TODO(bernhard): compute correct values using eq 33 Zhou-2021-Environ_Model_Softw.pdf
        #     δ18O_INFLI = p_fu_δ18O_SLFL
        #     δ2H_INFLI  = p_fu_δ2H_SLFL

        #     du_δ18O_GWAT, du_δ2H_GWAT, du_δ18O_SWATI, du_δ2H_SWATI =
        #         compute_isotope_du_GWAT_SWATI(
        #             # for GWAT:
        #             u_GWAT, u_δ18O_GWAT, u_δ2H_GWAT,
        #             # for SWATI:
        #             du_NTFLI, aux_du_VRFLI, aux_du_TRANI, aux_du_DSFLI, aux_du_INFLI, δ18O_INFLI, δ2H_INFLI,  # (non-fractionating)
        #             aux_du_SLVP, p_fT_TADTM, p_VAPPRES(t), p_δ2H_PREC(t), p_δ18O_PREC(t), u_aux_WETNES, # (fractionating)
        #             u_SWATI, u_δ18O_SWATI, u_δ2H_SWATI, 0, 0) #EffectiveDiffusivity_18O, EffectiveDiffusivity_2H)
        # end

        ##################
        # Update solution:
        # u is a state vector with u[1] = S relative saturation (-)
        # Update GWAT:
        du.GWAT.mm   = aux_du_VRFLI[NLAYER] - du_GWFL[1] - du_SEEP[1]

        # Do not modify INTS, INTR, SNOW, CC, SNOWLQ
        # as they are separately modified by the callback.
        # Therfore set their rate of change to 0
        du.INTS.mm   = 0
        du.INTR.mm   = 0
        du.SNOW.mm   = 0
        du.CC.MJm2   = 0
        du.SNOWLQ.mm = 0

        du.XYLEM.mm = 0
        du.RWU.mmday   = 0  # TODO(bernhard): these are rather accumulated quantities than proper state variables
        du.TRANI.mmday .= 0 # TODO(bernhard): these are rather accumulated quantities than proper state variables
        # du.XYLEM.d18O = 0 #TODO(bernhard): remove from here
        # du.RWU.d18O   = 0 #TODO(bernhard): remove from here
        # du.TRANI.d18O = 0 #TODO(bernhard): remove from here

        # Update SWATI for each layer:
        du.SWATI.mm .= du_NTFLI

        # save intermediate results from flow calculation into a cache
        # for efficient use in transport calculation (in callbacks)
        # p[3][3]      .=  [du_NTFLI  aux_du_VRFLI  aux_du_DSFLI  aux_du_INFLI  u_aux_WETNES]
        # p[3][4][1:2] .=  [du_GWFL, du_SEEP]
        # with the method @unpack this is automatically overwritten
        # TODO(bernhard): it seems it is not automatically overwritten. As workaround check:
        # p.du_NTFLI     === du_NTFLI
        # p.aux_du_VRFLI === aux_du_VRFLI
        # p.aux_du_DSFLI === aux_du_DSFLI
        # p.aux_du_INFLI === aux_du_INFLI
        # p.u_aux_WETNES === u_aux_WETNES
        # p.du_GWFL === du_GWFL
        # p.du_SEEP === du_SEEP


        ##########################################
        # Accumulate flows to compute daily sums
        # Note that below state variables serve only as accumulator but do not affect
        # the evolution of the system.
        if compute_intermediate_quantities
            # @assert length(p[1][4].row_idx_accum) == 25

            # 1) Either set daily sum if rate is constant throughout precipitation interval: p_DTP*(...)
            # 2) or then set daily sum to zero and use ODE to accumulate flow.
            du.accum.cum_d_prec     = 0 # was computed in callback
            du.accum.cum_d_rfal     = 0 # was computed in callback
            du.accum.cum_d_sfal     = 0 # was computed in callback
            du.accum.cum_d_rint     = 0 # was computed in callback
            du.accum.cum_d_sint     = 0 # was computed in callback
            du.accum.cum_d_rsno     = 0 # was computed in callback
            du.accum.cum_d_rnet     = 0 # was computed in callback
            du.accum.cum_d_smlt     = 0 # was computed in callback

            du.accum.cum_d_evap     = 0 # was computed in callback
            du.accum.cum_d_tran     = 0 # was computed in callback
            du.accum.cum_d_irvp     = 0 # was computed in callback
            du.accum.cum_d_isvp     = 0 # was computed in callback
            du.accum.cum_d_slvp     = 0 # was computed in callback
            du.accum.cum_d_snvp     = 0 # was computed in callback
            du.accum.cum_d_pint     = 0 # was computed in callback
            du.accum.cum_d_ptran    = 0 # was computed in callback
            du.accum.cum_d_pslvp    = 0 # was computed in callback

            du.accum.flow           = p_fu_SRFL +
                                        sum(aux_du_BYFLI) +
                                        sum(aux_du_DSFLI) +
                                        du_GWFL[1]
            du.accum.seep           = du_SEEP[1]
            du.accum.srfl           = p_fu_SRFL
            du.accum.slfl           = p_fu_SLFL
            du.accum.byfl           = sum(aux_du_BYFLI)
            du.accum.dsfl           = sum(aux_du_DSFLI)
            du.accum.gwfl           = du_GWFL[1]
            du.accum.vrfln          = aux_du_VRFLI[NLAYER]
            # du.accum.cum_d_rthr   = 0 # was computed in callback
            # du.accum.cum_d_sthr   = 0 # was computed in callback
            du.accum.totalSWAT      = 0 # is computed in callback
            du.accum.new_totalWATER = 0 # is computed in callback
            du.accum.BALERD_SWAT    = 0 # is computed in callback
            du.accum.BALERD_total   = 0 # is computed in callback

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
