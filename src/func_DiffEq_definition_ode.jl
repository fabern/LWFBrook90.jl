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

    # a)
    # variant a1) also check d18O and d2H, but make sure only to check those that are never NaN
    # [compartment[:, :] for compartment in ode.u0.x]
    # function unstable_check_function(dt,u,p,t) any(isnan, u.x) end # select only amounts in each compartment
    # TODO: unimplemented. Do this after switching to ComponentArrays.jl


    # variant a2) only check certain states that are not containing NaNs (δ values might be set to NaN
    # init_state = NamedTuple{ode.p[1][4].u0_field_names}(ode.u0.x)
    # [compartment[:, 1] for compartment in ode.u0.x[1]] # compartments are ode.p[1][4].u0_field_names
    @assert ode.p[1][4].u0_variable_names == (d18O = 2, d2H = 3)
    # function unstable_check_function(dt,u,p,t) any(isnan, VectorOfArray([compartment[:, 1] for compartment in u.x])...) end # select only amounts in each compartment
    function unstable_check_function1(dt,u,p,t) # select only amounts in each compartment, no concentrations
        # any.(isnan.(VectorOfArray([compartment[:, 1] for compartment in u.x])))
        # any(isnan, VectorOfArray([compartment[:, 1] for compartment in u.x])...)
        return_value = false
        x = [compartment[:, 1] for compartment in u.x]
        for n in eachindex(x)
            if any(isnan.(x[n]))
                return_value = true
                break
            end
        end
        return return_value
    end
    function unstable_check_function2(dt,u,p,t) # select only amounts in the first 12 compartments, no concentrations
        any(isnan.([u[1,1], u[2,1], u[3,1], u[4,1], u[5,1], u[6,1], u[7,1], u[8,1], u[9,1], u[10,1], u[11,1], u[12,1]]))
    end
    function unstable_check_function3(dt,u,p,t) # check all
        any(isnan.(u[:]))
    end
    # @btime unstable_check_function1(0.1,u0,p,t0)
    # @btime unstable_check_function2(0.1,u0,p,t0)
    # @btime unstable_check_function3(0.1,u0,p,t0)


    # b) fix adaptivity by only using a norm that considers amounts (1st column of state array), but no concentrations (2nd and 3rd column of state array)
    # Problem: When adding transport of isotopes adaptive time stepping has difficulties and reduces dt below dtmin.
    # Solution: It turns out this is linked to the norm that the adaptivitiy algorithm uses. It can be overruled
    #           by only using a norm that considers amounts (1st column of state array), but no concentrations (2nd and 3rd column of state array)
    #           The norm needs to be defined once for scalar elements of u and once for the whole array of u (multiple dispatch).
    function norm_to_use(u::Real,  t) norm(u) end
    # function norm_to_use(u::Array, t) DiffEqBase.ODE_DEFAULT_NORM(u[:,1], t) end
    # function norm_to_use(u::ArrayPartition, t) DiffEqBase.ODE_DEFAULT_NORM(VectorOfArray([compartment[:, 1] for compartment in u.x]), t) end
    function norm_to_use(u::ArrayPartition, t) DiffEqBase.ODE_DEFAULT_NORM([u[1,1],u[2,1],u[3,1],u[4,1],u[5,1],u[6,1],u[7,1],u[8,1],u[9,1],u[10,1],u[11,1],u[12,1]], t) end
        # TEST IT:
        # norm_to_use(ode.u0, 1.0)
        # methods(norm_to_use)
        # methods(DiffEqBase.ODE_DEFAULT_NORM)

    # 2) Solve the system
    tspan = ode.tspan
    saveat = tspan[1]:tspan[2]
    @time sol_LWFBrook90 = solve(ode, Tsit5();
        progress       = true,
        saveat         = saveat,
        unstable_check = unstable_check_function2,
        dt             = 1e-3, # dt is initial dt, but can be changed adaptively
        adaptive       = true, internalnorm = norm_to_use) # fix adaptivity norm for NAs

    # @info sol_LWFBrook90.destats
    @info "Time steps for solving: $(sol_LWFBrook90.destats.naccept) ($(sol_LWFBrook90.destats.naccept) accepted out of $(sol_LWFBrook90.destats.nreject + sol_LWFBrook90.destats.naccept) total)"
    # @show now()

    return sol_LWFBrook90
end