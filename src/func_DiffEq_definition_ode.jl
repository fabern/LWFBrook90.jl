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

    # Define ODE problem
    ode = ODEProblem(f_LWFBrook90R,
                     u0,
                     tspan,
                     p,
                     callback=cb_func)
    return ode
end