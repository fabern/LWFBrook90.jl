"""
    run_example()

Run example simulation located in "/example/BEA2016-reset-FALSE-input" for 100 days
and return a Dict containing the solution (solution
object of DifferentialEquations.jl) and other variables useful for plotting.

## Example:
    using LWFBrook90
    using Plots
    example = LWFBrook90.run_example()

    # Plot scalar solution
    # Using simple plot recipe that interpolates, but without dates
    plot(example["solution"];
        vars = [1, 2, 3, 4, 5, 6],
        label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

    # Plot vector solution
    x = example["solutionDates"]
    y = cumsum(example["thickness"])
    z = example["solution"][7 .+ (0:example["NLAYER"]-1), 1, :]./example["thickness"]
    heatmap(x, y, z,
        yflip = true,
        xlabel = "Date",
        ylabel = "Depth [mm]",
        colorbar_title = "θ [-]")

    # Plot both using plot recipe for LWFBrook90
    using Plots, Measures
    using LWFBrook90
    optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)

    pl_final = LWFBrook90.plotlwfbrook90(sol_LWFBrook90, optim_ticks)
    savefig(pl_final,
        input_prefix * "_plotRecipe_NLAYER" * string(sol_LWFBrook90.prob.p[1][1].NLAYER) * ".png")
"""
function run_example()

    @info """
    LWFBrook90 example is being run. Once it has finished, plotting can be done in the
    following way:

    using LWFBrook90
    using Plots
    example = LWFBrook90.run_example()

    # Plot scalar solution
    # Using simple plot recipe that interpolates, but without dates
    plot(example["solution"];
        vars = [1, 2, 3, 4, 5, 6],
        label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

    # Plot vector solution
    x = example["solutionDates"]
    y = cumsum(example["thickness"])
    z = example["solution"][7 .+ (0:example["NLAYER"]-1), 1, :]./example["thickness"]
    heatmap(x, y, z,
        yflip = true,
        xlabel = "Date",
        ylabel = "Depth [mm]",
        colorbar_title = "θ [-]")

    # Plot both using plot recipe for LWFBrook90
    using Plots, Measures
    using LWFBrook90
    optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)

    pl_final = LWFBrook90.plotlwfbrook90(sol_LWFBrook90, optim_ticks)
    savefig(pl_final,
        input_prefix * "_plotRecipe_NLAYER" * string(sol_LWFBrook90.prob.p[1][1].NLAYER) * ".png")
    """

    # 1a) Read in input data
    input_prefix = "BEA2016-reset-FALSE"
    # input_path = "example/"*input_prefix*"-input/"
    input_path = @__DIR__ # https://stackoverflow.com/a/63021629

    ####################
    (input_meteoveg,
        input_meteoveg_reference_date,
        input_param,
        input_storm_durations,
        input_initial_conditions,
        input_soil_horizons,
        input_soil_discretization,
        simOption_FLAG_MualVanGen) = read_inputData(input_path, input_prefix)
    ####################

    ####################
    # Define solver options
    Reset = false                          # currently only Reset = 0 implemented
    compute_intermediate_quantities = true # Flag whether ODE containes additional quantities than only states

    # Override input file settings
    # Here possibility to check and override dataframes input_[...] manually
    # # E.g:
    # # Soil hydraulic model
    # input_param[1,"NOOUTF"] = true # `true` if outflow from roots prevented, `false` if allowed
    ####################

    ####################
    # Define parameters for differential equation
    uSoil_initial, p = define_LWFB90_p(
        input_meteoveg,
        input_meteoveg_reference_date,
        input_param,
        input_storm_durations,
        input_soil_horizons,
        input_soil_discretization,
        simOption_FLAG_MualVanGen;
        Reset = Reset,
        # soil_output_depths = collect(-0.05:-0.05:-1.1),
        compute_intermediate_quantities = compute_intermediate_quantities)
    ####################

    ####################
    # Define initial states of differential equation
    # state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI
    # Create u0 for DiffEq.jl
    u0 = define_LWFB90_u0(p, input_initial_conditions,
        uSoil_initial,
        compute_intermediate_quantities)
    ####################

    ####################
    # Define ODE problem which consists of
    #   - definition of right-hand-side (RHS) function f
    #   - definition of callback function cb
    #   - u0:     initial condition of states
    #   - tspan:  definition of simulation time span
    #   - p:      parameters

    # Define simulation time span:
    tspan = (0.0, 100.0) # simulate 100 days
    ode_LWFBrook90 = define_LWFB90_ODE(u0, tspan, p)
    ####################

    ####################
    ## Solve ODE:
    sol_LWFBrook90 = solve(ode_LWFBrook90, Tsit5();
        progress = true,
        saveat = tspan[1]:tspan[2], dt = 1e-6, adaptive = true) # dt is initial dt, but adaptive
    ####################


    return (
        Dict(["solution" => sol_LWFBrook90,
        "solutionDates" => LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t,
            input_meteoveg_reference_date),
        "thickness" => sol_LWFBrook90.prob.p[1][1].p_THICK,
        "NLAYER" => sol_LWFBrook90.prob.p[1][1].NLAYER])
    )

end