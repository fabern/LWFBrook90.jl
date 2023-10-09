# a) Workaround for stability check with NaN in state vector
# variant a1) also check d18O and d2H, but make sure only to check those that are never NaN
# [compartment[:, :] for compartment in ode.u0.x]
# function unstable_check_function(dt,u,p,t) any(isnan, u.x) end # select only amounts in each compartment
# TODO: unimplemented. Do this after switching to ComponentArrays.jl
# variant a2) only check certain states that are not containing NaNs (δ values might be set to NaN
function unstable_check_function(dt,u::ComponentVector,p,t) # select only amounts in the first 12 compartments, no concentrations
    # any(isnan.([u.GWAT.mm, u.INTS.mm, u.INTR.mm, u.SNOW.mm, u.CC.MJm2, u.SNOWLQ.mm, u.SWATI.mm, u.RWU.mmday, u.XYLEM.mm, u.TRANI.mmday]))
    # # NOTE: `aux`` is missing on purpose
    # any(isnan, u.SWATI) # Do it for SWATI (both amt and concentration)
    any(isnan, u.SWATI.mm) # Do it for SWATI (only amt)
end

# b) Workaround for norm with NaN in state vector
function norm_to_use(u::Real,           t) norm(u)                                end
# function norm_to_use(u::Array, t) DiffEqBase.ODE_DEFAULT_NORM(u[:,1], t) end
# function norm_to_use(u::ArrayPartition, t) DiffEqBase.ODE_DEFAULT_NORM(u.x[1], t) end
# function norm_to_use(u::ArrayPartition, t) DiffEqBase.ODE_DEFAULT_NORM([u[1,1],u[2,1],u[3,1],u[4,1],u[5,1],u[6,1],u[7,1],u[8,1],u[9,1],u[10,1],    #u[11,1], # Not for aux...
#                                                                         u[12,1]], t) end
function norm_to_use(u::ComponentVector, t)
    # DiffEqBase.ODE_DEFAULT_NORM(u.SWATI, t)
    # NOTE: `aux`` is missing on purpose
    DiffEqBase.ODE_DEFAULT_NORM(reduce(ComponentArrays.vcat, # specifying Pkg needed until https://github.com/jonniedie/ComponentArrays.jl/issues/222 solved
                                    [u.GWAT.mm,
                                    u.INTS.mm,
                                    u.INTR.mm,
                                    u.SNOW.mm,
                                    u.CC.MJm2,
                                    u.SNOWLQ.mm,
                                    u.SWATI.mm,
                                    u.RWU.mmday,
                                    u.XYLEM.mm,
                                    u.TRANI.mmday,
                                    u.aux.θ, #u.aux.ψ, u.aux.K,
                                    u.accum]),
        t)
    # DiffEqBase.ODE_DEFAULT_NORM(vcat(u.GWAT.mm,
    #                                  u.INTS.mm,
    #                                  u.INTR.mm,
    #                                  u.SNOW.mm,
    #                                  u.CC.MJm2,
    #                                  u.SNOWLQ.mm,
    #                                  u.SWATI.mm,
    #                                  u.RWU.mmday,
    #                                  u.XYLEM.mm,
    #                                  u.TRANI.mmday,
    #                                  u.accum),
    # t)
end
