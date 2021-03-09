"""
    define_LWFB90_u0()

Generate vector u0 needed for ODE() problem in DiffEq.jl package.
"""
function define_LWFB90_u0(u_GWAT_init,
                            u_INTS_init,
                            u_INTR_init,
                            u_SNOW_init,
                            u_CC_init,
                            u_SNOWLQ_init,
                            u_SWATIinit,
                            compute_intermediate_quantities)

    if compute_intermediate_quantities
    # TODO(bernhard): are these store somewhere else than input_siteparam and pfile_param
        u0 = [u_GWAT_init;
            u_INTS_init;
            u_INTR_init;
            u_SNOW_init;
            u_CC_init;
            u_SNOWLQ_init;
            u_SWATIinit;
            # accumulation variables:
            0;0;0;0;0;0;0;0;0;0;
            0;0;0;0;0;0;0;0;0;0;
            0;0;0;0;0]
    else
        u0 = [u_GWAT_init;
                u_INTS_init;
                u_INTR_init;
                u_SNOW_init;
                u_CC_init;
                u_SNOWLQ_init;
                u_SWATIinit]
    end

    return u0
end