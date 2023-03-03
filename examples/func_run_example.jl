"""
    run_example()

Run example simulation located in "/examples/BEA2016-reset-FALSE-input" for 100 days
and return a Dict containing the solution (solution
object of DifferentialEquations.jl) and other variables useful for plotting.

## Example:
    using LWFBrook90
    using Plots, Measures
    example = LWFBrook90.run_example()

    ###
    # A) Use in-built plotting function
    optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
    pl_inbuilt = LWFBrook90.ISO.plotisotopes(
        example["solution"], optim_ticks;
        layout = grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4]),
        size=(1000,1400), dpi = 300, leftmargin = 15mm);
    savefig(pl_inbuilt, "Isotopeplots_pl_inbuilt.png")

    ###
    # B) Construct plots yourself using the solution object
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
        input_prefix * "_plotRecipe_NLAYER" * string(sol_LWFBrook90.prob.p.p_soil.NLAYER) * ".png")
"""
function run_example()

    @info """
    LWFBrook90 example is being run. Once it has finished, plotting can be done in the
    following way:

    using LWFBrook90
     using Plots, Measures
    example = LWFBrook90.run_example()

    ###
    # A) Use in-built plotting function
    optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
    pl_inbuilt = LWFBrook90.ISO.plotisotopes(
        example["solution"], optim_ticks;
        layout = grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4]),
        size=(1000,1400), dpi = 300, leftmargin = 15mm);
    savefig(pl_inbuilt, "Isotopeplots_pl_inbuilt.png")

    ###
    # B) Construct plots yourself using the solution object
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
        input_prefix * "_plotRecipe_NLAYER" * string(sol_LWFBrook90.prob.p.p_soil.NLAYER) * ".png")
    """

    # 1a) Read in input data
    input_prefix = "BEA2016-reset-FALSE"
    # input_prefix = "isoBEA2016-reset-FALSE"
    # input_path = "examples/"*input_prefix*"-input/"
    input_path = joinpath(@__DIR__, input_prefix*"-input/") # https://stackoverflow.com/a/63021629

    ####################
    # Define simulation model by reading in system definition and input data
    model = loadSPAC(input_path, input_prefix;
                 simulate_isotopes = contains(input_prefix, "iso"));
    ####################

    ####################
    # Prepare simulation by discretizing spatial domain
    simulation = LWFBrook90.setup(model; tspan = (0,100));

    # Solve ODE:
    LWFBrook90.simulate!(simulation)
    # plot(simulation.ODESolution)
    sol_LWFBrook90 = simulation.ODESolution
    ####################

    return (
        Dict(["solution" => sol_LWFBrook90,
        "solutionDates" => LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t,
            simulation.parametrizedSPAC.reference_date),
        "thickness" => sol_LWFBrook90.prob.p.p_soil.p_THICK,
        "NLAYER" => sol_LWFBrook90.prob.p.p_soil.NLAYER])
    )

end
