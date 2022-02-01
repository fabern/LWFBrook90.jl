"""
    define_LWFB90_u0()

Generate vector u0 needed for ODE() problem in DiffEq.jl package.
"""
function define_LWFB90_u0(p, uScalar_initial, uSoil_initial,
    compute_intermediate_quantities)

    # amounts
    u_GWAT_init_mm = uScalar_initial[1, "u_GWAT_init_mm"]
    u_INTS_init_mm = uScalar_initial[1, "u_INTS_init_mm"]
    u_INTR_init_mm = uScalar_initial[1, "u_INTR_init_mm"]
    u_SNOW_init_mm = uScalar_initial[1, "u_SNOW_init_mm"]
    u_CC_init_MJ_per_m2 = uScalar_initial[1, "u_CC_init_MJ_per_m2"]
    u_SNOWLQ_init_mm = uScalar_initial[1, "u_SNOWLQ_init_mm"]
    u_SWATIinit_mm = uSoil_initial

    u0 = [u_GWAT_init_mm
            u_INTS_init_mm
            u_INTR_init_mm
            u_SNOW_init_mm
            u_CC_init_MJ_per_m2
            u_SNOWLQ_init_mm
            u_SWATIinit_mm]

    # idx_u_vector_amounts       = p[1][4][4]
    # idx_u_vector_accumulators  = p[1][4][5]
    # idx_u_scalar_amounts       = p[1][4][6]

    if compute_intermediate_quantities
        # accumulation variables:
        N_accumulate_variables = 25
        u0 = [u0; fill(0, N_accumulate_variables, 1)]
        # Set idx_u_vector_accumulators
        # idx_u_vector_accumulators = (length(u0)-N_accumulate_variables+1):length(u0)
        #idx_u_vector_accumulators = 6+4+4+3*p[1][1].NLAYER .+ (1:25)
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