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
        if Reset != 1
            integrator.u[7+NLAYER+ 9] += DTINEW*(                        sum(aux_du_TRANI_corrected))
            integrator.u[7+NLAYER+12] += DTINEW*(aux_du_SLVP_corrected)
        end
        integrator.p[3][3] = aux_du_TRANI_corrected
        integrator.p[3][4] = aux_du_SLVP_corrected

        ##################
        # Update soil limited boundary flows during iteration loop
        if Reset == 1
            if (integrator.u[7+NLAYER+25] > 0)
                # starting from the 2nd iteration
                # update certain quantities based on recomputed state variables
                apply_reset(integrator, DTINEW)
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
