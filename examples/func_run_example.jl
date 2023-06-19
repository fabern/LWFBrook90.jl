"""
    run_example()

Run example simulation located in "/examples/BEA2016-reset-FALSE-input" for 100 days
and return a Dict containing the solution (solution
object of DifferentialEquations.jl) and other variables useful for plotting.

Kwargs can be provided and are passed through to loadSPAC().

## Example:
    using LWFBrook90
    example = LWFBrook90.run_example()

    ###
    using Plots, Measures
    # A) Use in-built plotting function
    pl1 = plotamounts(example, :above_and_belowground, :showRWUcentroid)
    pl2 = plotisotopes(example, :d18O)
    pl3 = plotforcingandstates(example)
    savefig(pl1, "Example_pl1.png")
    savefig(pl2, "Example_pl2.png")
    savefig(pl3, "Example_pl3.png")

    ###
    # B) Construct plots yourself using the solution object
    # Plot scalar solution
    # Using simple plot recipe that interpolates, but without dates
    plot(example.ODESolution;
        idxs = [1, 2, 3, 4, 5, 6],
        label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

    # Plot vector solution
    x = example.ODESolution_datetime
    y = cumsum(example.ODESolution.prob.p.p_soil.p_THICK)
    z = get_theta(example)
    heatmap(x, y, z,
        yflip = true,
        xlabel = "Date",
        ylabel = "Depth [mm]",
        c = cgrad(:blues,  rev = false),
        colorbar_title = "θ [m3/m3]\n(of fine soil volume)")#, rightmargin = 10mm)
"""
function run_example(;kwargs...)

    @info """
    LWFBrook90 example is being run. Once it has finished, plotting can be done in the
    following way:

    using LWFBrook90
    example = LWFBrook90.run_example()

    ###
    using Plots, Measures
    # A) Use in-built plotting function
    pl1 = plotamounts(example, :above_and_belowground, :showRWUcentroid)
    pl2 = plotisotopes(example, :d18O)
    pl3 = plotforcingandstates(example)
    savefig(pl1, "Example_pl1.png")
    savefig(pl2, "Example_pl2.png")
    savefig(pl3, "Example_pl3.png")

    ###
    # B) Construct plots yourself using the solution object
    # Plot scalar solution
    # Using simple plot recipe that interpolates, but without dates
    plot(example.ODESolution;
        idxs = [1, 2, 3, 4, 5, 6],
        label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

    # Plot vector solution
    x = example.ODESolution_datetime
    y = cumsum(example.ODESolution.prob.p.p_soil.p_THICK)
    z = get_theta(example)
    heatmap(x, y, z,
        yflip = true,
        xlabel = "Date",
        ylabel = "Depth [mm]",
        c = cgrad(:blues,  rev = false),
        colorbar_title = "θ [m3/m3]\n(of fine soil volume)")#, rightmargin = 10mm)
    """

    # 1a) Read in input data
    input_prefix = "DAV2020-full"
    # input_prefix = "isoBEA2016-reset-FALSE"
    # input_path = "examples/"*input_prefix*"-input/"
    input_path = joinpath(@__DIR__, input_prefix) # https://stackoverflow.com/a/63021629

    ####################
    # Define simulation model by reading in system definition and input data
    model = loadSPAC(input_path, input_prefix; simulate_isotopes = true, kwargs...);
    ####################

    ####################
    # Prepare simulation by discretizing spatial domain
    simulation = LWFBrook90.setup(model; requested_tspan = (0,300));

    # Solve ODE:
    LWFBrook90.simulate!(simulation)
    # plot(simulation.ODESolution)
    sol_LWFBrook90 = simulation.ODESolution
    ####################

    return simulation
    # return (
    #     Dict(["simulation" => simulation,
    #     "solutionDates" => LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t,
    #         simulation.parametrizedSPAC.reference_date),
    #     "thickness" => simulation.ODESolution.prob.p.p_soil.p_THICK,
    #     "NLAYER" => simulation.ODESolution.prob.p.p_soil.NLAYER])
    # )

end
