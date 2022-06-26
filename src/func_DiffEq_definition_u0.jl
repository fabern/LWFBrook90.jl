"""
    define_LWFB90_u0()

Generate vector u0 needed for ODE() problem in DiffEq.jl package.
"""
function define_LWFB90_u0(p, uScalar_initial,
    uSoil_initial_ψM_kPa, uSoil_initial_δ18O_permil, uSoil_initial_δ2H_permil,
    compute_intermediate_quantities;
    simulate_isotopes::Bool = false)

    # amounts:
    u_GWAT_init_mm      = uScalar_initial[1, "u_GWAT_init_mm"]
    u_INTS_init_mm      = uScalar_initial[1, "u_INTS_init_mm"]
    u_INTR_init_mm      = uScalar_initial[1, "u_INTR_init_mm"]
    u_SNOW_init_mm      = uScalar_initial[1, "u_SNOW_init_mm"]
    u_CC_init_MJ_per_m2 = uScalar_initial[1, "u_CC_init_MJ_per_m2"]
    u_SNOWLQ_init_mm    = uScalar_initial[1, "u_SNOWLQ_init_mm"]
    p_soil = p[1][1]
    # u_SWATIinit_mm    = LWFBrook90.KPT.FWETNES(uSoil_initial_ψM_kPa, p_soil) .*
    #                      p_soil.p_SWATMAX
    u_SWATIinit_mm      = LWFBrook90.KPT.FTheta(LWFBrook90.KPT.FWETNES(uSoil_initial_ψM_kPa, p_soil), p_soil) .*
                          p_soil.p_SWATMAX ./ p_soil.p_THSAT # see l.2020: https://github.com/pschmidtwalter/LWFBrook90R/blob/6f23dc1f6be9e1723b8df5b188804da5acc92e0f/src/md_brook90.f95#L2020

    u0 = [u_GWAT_init_mm
        u_INTS_init_mm
        u_INTR_init_mm
        u_SNOW_init_mm
        u_CC_init_MJ_per_m2
        u_SNOWLQ_init_mm
        u_SWATIinit_mm]

    # idx_u_vector_amounts       = p[1][4][4]
    idx_u_vector_accumulators  = p[1][4][5]
    # idx_u_scalar_amounts       = p[1][4][6]

    # isotopes:
    if simulate_isotopes
        # isotopes d18O
        u_GWAT_init_d18O      = uScalar_initial[2, "u_GWAT_init_mm"]
        u_INTS_init_d18O      = uScalar_initial[2, "u_INTS_init_mm"]
        u_INTR_init_d18O      = uScalar_initial[2, "u_INTR_init_mm"]
        u_SNOW_init_d18O      = uScalar_initial[2, "u_SNOW_init_mm"]
        # u_CC_init_MJ_per_m2_d18O = uScalar_initial[2, "u_CC_init_MJ_per_m2"] # state has no isotopic signature
        # u_SNOWLQ_init_d18O       = uScalar_initial[2, "u_SNOWLQ_init_mm"]    # state has no isotopic signature
        u_SWATIinit_d18O      = uSoil_initial_δ18O_permil
        u_RWUinit_d18O        = uSoil_initial_δ18O_permil[1]  # TODO(bernhard): start out with same concentration as in first soil layer
        u_XYLEMinit_d18O      = uSoil_initial_δ18O_permil[1]  # TODO(bernhard): start out with same concentration as in first soil layer

        # isotopes d2H
        u_GWAT_init_d2H      = uScalar_initial[3, "u_GWAT_init_mm"]
        u_INTS_init_d2H      = uScalar_initial[3, "u_INTS_init_mm"]
        u_INTR_init_d2H      = uScalar_initial[3, "u_INTR_init_mm"]
        u_SNOW_init_d2H      = uScalar_initial[3, "u_SNOW_init_mm"]
        # u_CC_init_MJ_per_m2_d2H = uScalar_initial[3, "u_CC_init_MJ_per_m2"] # state has no isotopic signature
        # u_SNOWLQ_init_d2H       = uScalar_initial[3, "u_SNOWLQ_init_mm"]    # state has no isotopic signature
        u_SWATIinit_d2H      = uSoil_initial_δ2H_permil
        u_RWUinit_d2H        = uSoil_initial_δ2H_permil[1]  # TODO(bernhard): start out with same concentration as in first soil layer
        u_XYLEMinit_d2H      = uSoil_initial_δ2H_permil[1]  # TODO(bernhard): start out with same concentration as in first soil layer


        u0 = [u0;
            u_GWAT_init_d18O;
            u_INTS_init_d18O;
            u_INTR_init_d18O;
            u_SNOW_init_d18O;
            # u_CC_init_d18O;
            # u_SNOWLQ_init_d18O;
            u_RWUinit_d18O;
            u_XYLEMinit_d18O
            u_SWATIinit_d18O;

            u_GWAT_init_d2H;
            u_INTS_init_d2H;
            u_INTR_init_d2H;
            u_SNOW_init_d2H;
            # u_CC_init_d2H;
            # u_SNOWLQ_init_d2H;
            u_RWUinit_d2H;
            u_XYLEMinit_d2H;
            u_SWATIinit_d2H]

            # Nδ = 6
            # idx_u_vector_isotopes_d18O = length(idx_u_vector_amounts) .+ (1:(Nδ + length(u_SWATIinit_d18O)))
            # idx_u_vector_isotopes_d2H  = length(idx_u_vector_amounts) .+ (1:(Nδ + length(u_SWATIinit_d2H))) .+ Nδ .+ length(u_SWATIinit_d18O)

            #idx_u_vector_isotopes_d18O = 6+p[1][1].NLAYER     .+ (1:(Nδ+p[1][1].NLAYER))
            #idx_u_vector_isotopes_d2H  = 6+Nδ+2*p[1][1].NLAYER .+ (1:(Nδ+p[1][1].NLAYER))
    else
        #TODO(bernhard): In order to set idx_u_vector_isotopes_d18O and idx_u_vector_isotopes_d2H,
        #                we need to solve the follwoing:
        #                How to define a quantity that can return an empty vector when within []
        #                    u0[idx_u_vector_amounts]
        #                    u0[idx_u_vector_isotopes_d18O]
        #                    u0[idx_u_vector_isotopes_d2H]
        #                    u0[1:2]
        #                    u0[1]
        #                    u0[false]
    end

    if compute_intermediate_quantities
        # accumulation variables:
        N_accumulate_variables = 31
        u0 = [u0; fill(0, N_accumulate_variables, 1)]
        # Set idx_u_vector_accumulators
        # idx_u_vector_accumulators = (length(u0)-N_accumulate_variables+1):length(u0)
        #idx_u_vector_accumulators = 6+Nδ+Nδ+3*p[1][1].NLAYER .+ (1:31)

        # initialize terms for balance errors with initial values from u0
        new_SWAT       = sum(u_SWATIinit_mm) # total soil water in all layers, mm
        new_totalWATER = u_INTR_init_mm + u_INTS_init_mm + u_SNOW_init_mm + new_SWAT + u_GWAT_init_mm # total water in all compartments, mm

        u0[idx_u_vector_accumulators[28]] = new_SWAT
        u0[idx_u_vector_accumulators[29]] = new_totalWATER

    else
        #TODO(bernhard): In order to set idx_u_vector_accumulators, we need to solve the follwoing:
        #                How to define a quantity that can return an empty vector when within []
        #                    u0[idx_u_vector_amounts]
        #                    u0[idx_u_vector_isotopes_d18O]
        #                    u0[idx_u_vector_isotopes_d2H]
        #                    u0[1:2]
        #                    u0[1]
        #                    u0[false]
    end

    return u0
end