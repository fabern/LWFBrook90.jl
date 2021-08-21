"""
    define_LWFB90_u0()

Generate vector u0 needed for ODE() problem in DiffEq.jl package.
"""
function define_LWFB90_u0(u_GWAT_init_mm,
                            u_INTS_init_mm,
                            u_INTR_init_mm,
                            u_SNOW_init_mm,
                            u_CC_init_MJ_per_m2,
                            u_SNOWLQ_init_mm,
                            u_SWATIinit,
                            compute_intermediate_quantities)

    if compute_intermediate_quantities
        u0 = [u_GWAT_init_mm;
            u_INTS_init_mm;
            u_INTR_init_mm;
            u_SNOW_init_mm;
            u_CC_init_MJ_per_m2;
            u_SNOWLQ_init_mm;
            u_SWATIinit;
            # accumulation variables:
            0;0;0;0;0;0;0;0;0;0;
            0;0;0;0;0;0;0;0;0;0;
            0;0;0;0;0]
    else
        u0 = [u_GWAT_init_mm;
                u_INTS_init_mm;
                u_INTR_init_mm;
                u_SNOW_init_mm;
                u_CC_init_MJ_per_m2;
                u_SNOWLQ_init_mm;
                u_SWATIinit]
    end

    return u0
end