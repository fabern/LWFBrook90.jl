"""
    define_LWFB90_ODE()

Generates an ODEProblem from DiffEq.jl

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
function define_LWFB90_ODE(u0, tspan, p)

    # Define callback functions
    cb_func = define_LWFB90_cb()

    # swcheck_cb = ContinuousCallback()    #TODO(bernhard) Implement swchek as ContinuousCallback
    # Reset_cb = FunctionCallingCallback() #TODO(bernhard) low priority: implement Reset==1

    @assert !any(ismissing.(u0)) """
    There are missing values in the provided initial conditions `u0`. Please correct!"""

    # Define ODE problem
    ode = ODEProblem(LWFBrook90.f_LWFBrook90R,
                    u0,
                    tspan,
                    p,
                    callback=cb_func)


    # Define stability check for numerical solver that accepts NAs as undefined isotopic
    # concentrations (e.g. when a compartment such as SNOW is empty its concentration
    # is undefined)
    # unstable_check_function = (dt,u,p,t) -> false
    simulate_isotopes = p[1][4].simulate_isotopes
    if simulate_isotopes
        # only check certain states that are not containing NaNs (δ values might be set to NaN)
        unstable_check_function = (dt,u,p,t) -> any(isnan, u[[p[1][4].row_idx_SWATI; p[1][4].row_idx_scalars...], 1])         # select only amts
        # unstable_check_function = (dt,u,p,t) -> any(isnan, u[[p[1][4].row_idx_SWATI; p[1][4].row_idx_scalars...], [2,3]] ]) # select amts, d18O, d2H
    else
        # check all states
        unstable_check_function = (dt,u,p,t) -> any(isnan,u)
    end

    return ode, unstable_check_function
end