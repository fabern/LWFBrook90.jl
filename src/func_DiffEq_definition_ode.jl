function define_DiffEq_ODE(u0, tspan, p)

    # Define callback functions
    cb_func_daily = define_DiffEq_daily_cb()
    cb_func_every_timestep = define_DiffEq_timestep_cb()

    # swcheck_cb = ContinuousCallback()    #TODO(bernhard) Implement swchek as ContinuousCallback
    # Reset_cb = FunctionCallingCallback() #TODO(bernhard) low priority: implement Reset==1

    # Define ODE problem
    ode = ODEProblem(f_LWFBrook90R,
                     u0,
                     tspan,
                     p,
                     callback=CallbackSet(cb_func_daily, cb_func_every_timestep))
    return ode
end