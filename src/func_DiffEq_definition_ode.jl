"""
    solve_LWFB90(u0, tspan, p)

Generates an ODEProblem from DiffEq.jl and solves it over period tspan using initial
conditions u0. RHS f() and callbacks cb() are defined internally to LWFBrook90.jl

An ODE problem which consists of
    - definition of right-hand-side (RHS) function f
    - definition of callback function cb
    - initial condition of states
    - definition of simulation time span
    - parameters

Seperate updating of different states (INTS, INTR, SNOW, CC, SNOWLQ are updated once per
day while GWAT and SWATI are updated continuously) is implemented by means of operator
splitting using a callback function for the daily updates and a ODE RHS (right hand
side) for the continuous update.
"""
function solve_LWFB90(u0, tspan, p)
    # 1) Generate the ode problem and solve it with DiffEq.jl
    cb_func = define_LWFB90_cb() # define callback functions
    @assert !any(ismissing.(u0)) """
    There are missing values in the provided initial conditions `u0`. Please correct!"""

    # Define ODE problem which consists of
    #   - definition of right-hand-side (RHS) function f
    #   - definition of callback function cb
    #   - u0:     initial condition of states
    #   - tspan:  definition of simulation time span
    #   - p:      parameters
    ode_LWFBrook90 =
        ODEProblem(LWFBrook90.f_LWFBrook90R,
                    u0,
                    tspan,
                    p,
                    callback=cb_func)

    solve_LWFB90(ode_LWFBrook90)
end
function solve_LWFB90(ode::ODEProblem)
   # 1) Define workarounds to accomodate a state vector (array) that contains NAs
    # E.g. concentrations are undefined when a compartment is empty (e.g. SNOW = 0mm, δ_SNOW = NA)
    # a) stability check
    # b) norm for adaptive time stepping
    # see function definitions of unstable_check_function and norm_to_use

    # 2) Solve the system
    tspan = ode.tspan
    saveat = tspan[1]:tspan[2]
    @show saveat
    @time sol_LWFBrook90 = solve(ode,
        progress = true,
        saveat = saveat, save_everystep = true,
        # ImplicitEuler(autodiff=false);  # for stiff problems (~6.7s for 2.5days of Hammel_loam-NLayer-103)
        # Rodas4P(autodiff=false);        # for stiff problems (~3.2s for 2.5days of Hammel_loam-NLayer-103)
        # TRBDF2(autodiff=false);
        # Rosenbrock23(autodiff=false);   # for stiff problems   (~2.3s for 2.5days of Hammel_loam-NLayer-103)
        # Tsit5(); abstol = 1e-6, reltol = 1e-4, # Tsit5 recommended for non-stiff problems (~ 12.4s for 2.5days of Hammel_loam-NLayer-103)
        # Tsit5(); abstol = 1e-6, reltol = 1e-5, # Tsit5 recommended for non-stiff problems (~ 20.9s for 2.5days of Hammel_loam-NLayer-103)
        # Tsit5(); # Tsit5 recommended for non-stiff problems (~ 6.4s for 2.5days of Hammel_loam-NLayer-103)
        # AutoTsit5(Rosenbrock23(autodiff=false)); reltol = 1e-5, # recommended for problems of unknown stiffnes: (~5.0s)
        Tsit5(); reltol = 1e-5,
        adaptive = true, internalnorm = LWFBrook90.norm_to_use, # fix adaptivity norm for NAs
        unstable_check = LWFBrook90.unstable_check_function,         # fix instability norm for NAs
        dt    = 1e-3,                            # dt is initial dt, but can be changed adaptively
        dtmax = 60/60/24, # 60min max time step
        # reltol = 1e-5, # abstol = 1e-6, # default: abstol = 1e-6, reltol = 1e-3
        # tstops = sol_working.t, #[0.223, 0.4],
        # tstops = tspan[1]:0.00001:tspan[2], adaptive = false
    );

    # @info sol_LWFBrook90.destats
    @info "Time steps for solving: $(sol_LWFBrook90.destats.naccept) ($(sol_LWFBrook90.destats.naccept) accepted out of $(sol_LWFBrook90.destats.nreject + sol_LWFBrook90.destats.naccept) total)"
    # @show now()

    return sol_LWFBrook90
end



# a)
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

function norm_to_use(u::Real,           t) norm(u)                                end
# function norm_to_use(u::Array, t) DiffEqBase.ODE_DEFAULT_NORM(u[:,1], t) end
# function norm_to_use(u::ArrayPartition, t) DiffEqBase.ODE_DEFAULT_NORM(u.x[1], t) end
# function norm_to_use(u::ArrayPartition, t) DiffEqBase.ODE_DEFAULT_NORM([u[1,1],u[2,1],u[3,1],u[4,1],u[5,1],u[6,1],u[7,1],u[8,1],u[9,1],u[10,1],    #u[11,1], # Not for aux...
#                                                                         u[12,1]], t) end
function norm_to_use(u::ComponentVector, t)
    # DiffEqBase.ODE_DEFAULT_NORM(u.SWATI, t)
    # NOTE: `aux`` is missing on purpose
    DiffEqBase.ODE_DEFAULT_NORM(reduce(vcat, [u.GWAT.mm,
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
                                                u.accum
            ]),
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
