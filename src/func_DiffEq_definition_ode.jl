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
    # function unstable_check_function(dt,u::ArrayPartition,p,t) # select only amounts in the first 12 compartments, no concentrations
    #     # any(isnan.([u[1,1], u[2,1], u[3,1], u[4,1], u[5,1], u[6,1], u[7,1], u[8,1], u[9,1], u[10,1], u[11,1], u[12,1]]))
    #     any(isnan, u.x[1]) | any(isnan, u.x[2]) # 2 is for GWAT, which is not needed for stability
    # end
    function unstable_check_function(dt,u::ArrayPartition,p,t) # select only amounts in the first 12 compartments, no concentrations
        # any(isnan.([u[1,1], u[2,1], u[3,1], u[4,1], u[5,1], u[6,1], u[7,1], u[8,1], u[9,1], u[10,1], u[11,1], u[12,1]]))
        # any(isnan, u.x[1]) | any(isnan, u.x[2]) # 2 is for GWAT, which is not needed for stability
        any(isnan, u.x[7]) # Do it for SWATI (both amt and concentration)
        # any(isnan, reduce(vcat,[u.x[idx][:,1] for idx = 1:12]))
        # any(isnan, reduce(vcat,[u.x[1][:,1],
        #                     u.x[2][:,1],
        #                     u.x[3][:,1],
        #                     u.x[4][:,1],
        #                     u.x[5][:,1],
        #                     u.x[6][:,1],
        #                     u.x[7][:,1],
        #                     u.x[8][:,1],
        #                     u.x[9][:,1],
        #                     # 11 is missing on purpose
        #                     u.x[10][:,1],
        #                     u.x[12][:,1]]))
    end

    # b) fix adaptivity by only using a norm that considers amounts (1st column of state array), but no concentrations (2nd and 3rd column of state array)
    # Problem: When adding transport of isotopes adaptive time stepping has difficulties and reduces dt below dtmin.
    # Solution: It turns out this is linked to the norm that the adaptivitiy algorithm uses. It can be overruled
    #           by only using a norm that considers amounts (1st column of state array), but no concentrations (2nd and 3rd column of state array)
    #           The norm needs to be defined once for scalar elements of u and once for the whole array of u (multiple dispatch).
    # for more background compare:
    # a) https://diffeq.sciml.ai/stable/extras/timestepping/#timestepping and
    # b) https://github.com/SciML/SimpleDiffEq.jl/blob/08decf2e22f49ab2aa9b2995759d60fad709e049/src/tsit5/atsit5.jl#L210

    function norm_to_use(u::Real,           t) norm(u)                                end
    # function norm_to_use(u::ArrayPartition, t) DiffEqBase.ODE_DEFAULT_NORM(u.x[1], t) end
    # function norm_to_use(u::ArrayPartition, t) DiffEqBase.ODE_DEFAULT_NORM([u[1,1],u[2,1],u[3,1],u[4,1],u[5,1],u[6,1],u[7,1],u[8,1],u[9,1],u[10,1],u[11,1],u[12,1]], t) end
    # function norm_to_use(u::ArrayPartition, t) DiffEqBase.ODE_DEFAULT_NORM([u[1,1],u[2,1],u[3,1],u[4,1],u[5,1],u[6,1],u[7,1],u[8,1],u[9,1],u[10,1],    #u[11,1], # Not for aux...
    #                                                                         u[12,1]], t) end
    function norm_to_use(u::ArrayPartition, t)
        DiffEqBase.ODE_DEFAULT_NORM(reduce(vcat,[u.x[1][:,1],
                                                u.x[2][:,1],
                                                u.x[3][:,1],
                                                u.x[4][:,1],
                                                u.x[5][:,1],
                                                u.x[6][:,1],
                                                u.x[7][:,1],
                                                u.x[8][:,1],
                                                u.x[9][:,1],
                                                u.x[10][:,1],
                                                # 11 (:aux) is missing on purpose
                                                u.x[12][:,1]]),
                                    t)
    end
    # function norm_to_use(u::Array, t) DiffEqBase.ODE_DEFAULT_NORM(u[:,1], t) end


    # 2) Solve the system
    tspan = ode.tspan
    saveat = tspan[1]:tspan[2]
    @time sol_LWFBrook90 = solve(ode, Tsit5();
        progress       = true,
        saveat         = saveat, save_everystep = true,
        unstable_check = unstable_check_function,
        # dt = 2/60/24, adaptive = false,
        dt             = 1e-3, # dt is initial dt, but can be changed adaptively
        adaptive       = true, internalnorm = norm_to_use
        ) # fix adaptivity norm for NAs

    # @info sol_LWFBrook90.destats
    @info "Time steps for solving: $(sol_LWFBrook90.destats.naccept) ($(sol_LWFBrook90.destats.naccept) accepted out of $(sol_LWFBrook90.destats.nreject + sol_LWFBrook90.destats.naccept) total)"
    # @show now()

    return sol_LWFBrook90
end