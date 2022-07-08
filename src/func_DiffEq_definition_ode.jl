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

    # 2) Define workarounds to accomodate a state vector (array) that contains NAs
    # E.g. concentrations are undefined when a compartment is empty (e.g. SNOW = 0mm, δ_SNOW = NA)
    # a) stability check
    # b) norm for adaptive time stepping

    # a)
    # simulate_isotopes = p[1][4].simulate_isotopes
    # if simulate_isotopes
    #     # only check certain states that are not containing NaNs (δ values might be set to NaN)
    #     # # unstable_check_function = (dt,u,p,t) -> any(isnan, u[[p[1][4].row_idx_SWATI; p[1][4].row_idx_scalars...], [1,2,3]] ]) # select amts, d18O, d2H
    #     # unstable_check_function = (dt,u,p,t) -> any(isnan, u[[p[1][4].row_idx_SWATI; p[1][4].row_idx_scalars...], 1])         # select only amts of SWAT

    #     unstable_check_function = (dt,u,p,t) -> any(isnan, u[:, 1]) # select only amounts
    # else
    #     # check all states
    #     unstable_check_function = (dt,u,p,t) -> any(isnan,u)
    # end
    function unstable_check_function(dt,u,p,t) any(isnan, u[:, 1]) end # select only amounts

    # b) fix adaptivity by only using a norm that considers amounts (1st column of state array), but no concentrations (2nd and 3rd column of state array)
    # Problem: When adding transport of isotopes adaptive time stepping has difficulties and reduces dt below dtmin.
    # Solution: It turns out this is linked to the norm that the adaptivitiy algorithm uses. It can be overruled byfl
    #           by only using a norm that considers amounts (1st column of state array), but no concentrations (2nd and 3rd column of state array)
    #           The norm needs to be defined once for scalar elements of u and once for the whole array of u (multiple dispatch).
    function norm_to_use(u::Real,  t) norm(u) end
    function norm_to_use(u::Array, t) DiffEqBase.ODE_DEFAULT_NORM(u[:,1], t) end
        # norm_to_use(u0, integrator.t)
        # norm(u0[:,1])
        # methods(norm_to_use)
        # methods(DiffEqBase.ODE_DEFAULT_NORM)

    # 3) Solve the system
    saveat = tspan[1]:tspan[2]
    @time sol_LWFBrook90 = solve(ode_LWFBrook90, Tsit5();
        progress       = true,
        saveat         = saveat,
        unstable_check = unstable_check_function,
        dt             = 1e-3, # dt is initial dt, but can be changed adaptively
        adaptive       = true, internalnorm = norm_to_use) # fix adaptivity norm for NAs

    # @info sol_LWFBrook90.destats
    @info "Time steps for solving: $(sol_LWFBrook90.destats.naccept) ($(sol_LWFBrook90.destats.naccept) accepted out of $(sol_LWFBrook90.destats.nreject + sol_LWFBrook90.destats.naccept) total)"
    @show now()

    return sol_LWFBrook90
end