function define_DiffEq_ODE(u0, tspan, p)

    # Define callback functions
    cb_func = define_DiffEq_cb()

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