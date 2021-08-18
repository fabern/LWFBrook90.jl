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
    # Using dates (but not interpolated)
    plot(example["solutionDates"],
        example["solution"][[1,2,3,4,5,6],:]',
        label=["GWAT" "INTS" "INTR" "SNOW" "CC" "SNOWLQ"])

    # Using simple plot recipe that interpolates, but without dates
    plot(example["solution"];
        vars = [1, 2, 3, 4, 5, 6],
        label=["GWAT" "INTS" "INTR" "SNOW" "CC" "SNOWLQ"])

    # Plot vector solution
    x = example["solutionDates"]
    y = cumsum(example["thickness"])
    z = example["solution"][7 .+ (0:example["NLAYER"]-1), :]./example["thickness"]
    heatmap(x, y, z,
        yflip = true,
        xlabel = "Date",
        ylabel = "Depth [mm]",
        colorbar_title = "θ [-]")

"""
function run_example()

    @info """
    LWFBrook90 example is being run. Once it has finished, plotting can be done in the
    following way:

    using LWFBrook90
    using Plots
    example = LWFBrook90.run_example()

    # Plot scalar solution
    # Using dates (but not interpolated)
    plot(example["solutionDates"],
        example["solution"][[1,2,3,4,5,6],:]',
        label=["GWAT" "INTS" "INTR" "SNOW" "CC" "SNOWLQ"])

    # Using simple plot recipe that interpolates, but without dates
    plot(example["solution"];
        vars = [1, 2, 3, 4, 5, 6],
        label=["GWAT" "INTS" "INTR" "SNOW" "CC" "SNOWLQ"])

    # Plot vector solution
    x = example["solutionDates"]
    y = cumsum(example["thickness"])
    z = example["solution"][7 .+ (0:example["NLAYER"]-1), :]./example["thickness"]
    heatmap(x, y, z,
        yflip = true,
        xlabel = "Date",
        ylabel = "Depth [mm]",
        colorbar_title = "θ [-]")

    """

    # 1a) Read in input data
    input_prefix = "BEA2016-reset-FALSE"
    # input_path = "example/"*input_prefix*"-input/"
    input_path = @__DIR__ # https://stackoverflow.com/a/63021629

    ####################
    (input_meteoveg,
        input_meteoveg_reference_date,
        input_param,
        input_pdur,
        input_soil_materials,
        input_soil_nodes) = read_LWFBrook90R_inputData(input_path, input_prefix)
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
    p = define_LWFB90_p(
        input_meteoveg,
        input_meteoveg_reference_date,
        input_param,
        input_pdur,
        input_soil_materials,
        input_soil_nodes;
        Reset = Reset,
        compute_intermediate_quantities = compute_intermediate_quantities)
    ####################

    ####################
    # Define initial states of differential equation
    # state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI

    ######
    # Transform initial value of auxiliary state u_aux_PSIM_init into state u_SWATIinit:
    u_aux_PSIM_init = input_soil_nodes[:,"psiini"]
    if any( u_aux_PSIM_init.> 0)
        error("Initial matrix psi must be negative or zero")
    end
    p_soil = p[1][1]
    u_aux_WETNESinit = LWFBrook90.KPT.FWETNES(u_aux_PSIM_init, p_soil)
    u_SWATIinit      = p_soil.p_SWATMAX ./ p_soil.p_THSAT .* LWFBrook90.KPT.FTheta(u_aux_WETNESinit, p_soil)
    ######

    # Create u0 for DiffEq.jl
    u0 = define_LWFB90_u0(input_param[1,"u_GWAT_init"],
                        input_param[1,"u_INTS_init"],
                        input_param[1,"u_INTR_init"],
                        input_param[1,"u_SNOW_init"],
                        input_param[1,"u_CC_init"],
                        input_param[1,"u_SNOWLQ_init"],
                        u_SWATIinit,
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
    tspan = (0.,  100.) # simulate 100 days
    ode_LWFBrook90 = define_LWFB90_ODE(u0, tspan, p)
    ####################

    ####################
    ## Solve ODE:
    sol_LWFBrook90 = solve(ode_LWFBrook90, Tsit5();
        progress = true,
        saveat = tspan[1]:tspan[2], dt=1e-6, adaptive = true); # dt is initial dt, but adaptive
    ####################


    return (
        Dict(["solution"      => sol_LWFBrook90,
              "solutionDates" => LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t,
                                                                        input_meteoveg_reference_date),
              "thickness" => p_soil.p_THICK,
              "NLAYER"    => p_soil.NLAYER])
        )

end