"""
    run_example()

Run example simulation located "/example/BEA2016-reset-FALSE-input" for 100 days
and return a Dict containing the solution (solution
object of DifferentialEquations.jl) other variables useful for plotting.

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
        ylabel = "Depth",
        colorbar_title = "θ")

"""
function run_example()

    @info """
    LWFBrook90 example is run. Plotting can be done in the following way:

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
        ylabel = "Depth",
        colorbar_title = "θ")

    """

    # 1a) Read in input data
    input_prefix = "BEA2016-reset-FALSE"
    input_path = "example/"*input_prefix*"-input/"

    (input_meteo, input_param, input_siteparam, input_precdat, input_pdur,
    input_soil_materials, input_soil_nodes, input_reference_date) =
        read_LWFBrook90R_inputData(input_path, input_prefix)

    # 1c) Parse loaded/redefined input files
    (pfile_meteo, pfile_param, pfile_siteparam, pfile_precdat, pfile_pdur, pfile_soil) =
        derive_params_from_inputData(
            input_meteo, input_param, input_siteparam, input_precdat, input_pdur,
            input_soil_materials, input_soil_nodes, input_reference_date)

    ####################
    # Define simulation
    # Soil hydraulic model
    IMODEL = pfile_param["IMODEL"] # 0 for Clapp-Hornberger; 1 for Mualem-van Genuchten
    NLAYER = pfile_param["NLAYER"]

    # Define solver options
    NOOUTF    = 1 == pfile_param["NOOUTF"] # 1 if no outflow allowed from roots, otherwise 0
    Reset     = 0 # currently only Reset = 0 implemented

    constant_dt_solver = 1 # [days]
    # Flag whether ODE containes additional quantities than only states
    compute_intermediate_quantities = true
    ####################

    ####################
    # Define parameters for differential equation
    p = define_LWFB90_p(NLAYER, IMODEL, constant_dt_solver,
                        NOOUTF, Reset, compute_intermediate_quantities,
                        pfile_meteo,
                        pfile_siteparam,
                        pfile_param,
                        pfile_soil,
                        pfile_pdur)
    ####################

    ####################
    # Define initial states of differential equation
    # state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI
    u_GWAT_init = pfile_siteparam["u_GWAT_init"]
    u_SNOW_init = pfile_siteparam["u_SNOW_init"]
    u_INTS_init = pfile_param["INTS_init"];
    u_INTR_init = pfile_param["INTR_init"];
    u_CC_init     = 0; # any initial snow has zero liquid water and cold content
    u_SNOWLQ_init = 0; # any initial snow has zero liquid water and cold content

    u_aux_PSIM_init = pfile_soil["PSIM_init"]

    ######
    # Transform initial value of auxiliary state u_aux_PSIM_init into state u_SWATIinit:
    u_SWATIinit  = fill(NaN, NLAYER)
    if any( u_aux_PSIM_init.> 0)
        error("Initial matrix psi must be negative or zero")
    end
    if IMODEL == 0
        error("IMODEL==0 is not implemented to get initial SWATI.")
        # TODO(bernhard): implement this.
    elseif IMODEL == 1
        p_SWATMX = p[1][1][6]
        # TODO(bernhard): this hardcoded index is dangerous in case definition of p vector changes
        p_MvGα   = pfile_soil["PAR"][!,"α"]
        p_MvGn   = pfile_soil["PAR"][!,"n"]
        p_THSAT  = pfile_soil["PAR"][!,"θs"]
        p_θr     = pfile_soil["PAR"][!,"θr"]
        # TODO(bernhard): store above quantities somewhere else than p[1][1][3] and pfile_soil?
        for i = 1:NLAYER
            # Define initial u_SWATI based on input parameter
            u_aux_WETNESinit_i = LWFBrook90.KPT.FWETNES_MvG(u_aux_PSIM_init[i],
                                                                p_MvGα[i],
                                                                p_MvGn[i])
            u_SWATIinit[i]     = p_SWATMX[i]/p_THSAT[i] *
                                LWFBrook90.KPT.FTheta_MvG(u_aux_WETNESinit_i,
                                                                p_THSAT[i],
                                                                p_θr[i])

        end
    end
    ######

    # Create u0 for DiffEq.jl
    u0 = define_LWFB90_u0(u_GWAT_init,
                        u_INTS_init,
                        u_INTR_init,
                        u_SNOW_init,
                        u_CC_init,
                        u_SNOWLQ_init,
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
        saveat = tspan[1]:tspan[2], dt=1e-6, adaptive = true); # dt will be overwritten
    ####################


    return (
        Dict(["solution"      => sol_LWFBrook90,
              "solutionDates" => LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t,
                                                                        input_reference_date),
              "thickness" => pfile_soil["THICK"],
              "NLAYER"    => NLAYER])
        )

end