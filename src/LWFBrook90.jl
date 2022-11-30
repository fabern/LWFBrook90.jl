module LWFBrook90

using SciMLBase       # instead of loading the full DifferentialEquations
using OrdinaryDiffEq  # instead of loading the full DifferentialEquations
using DiffEqCallbacks # instead of loading the full DifferentialEquations
using RecipesBase, PlotUtils, Measures
using LinearAlgebra
using StatsBase: mean, weights
using ComponentArrays
using UnPack: @unpack
using Dates: now

# TODO(bernhard): make sure we have documentation for these exported variables
export SPAC, DiscretizedSPAC, discretize, simulate!
export run_simulation
export Rootden_beta_
export RelativeDaysFloat2DateTime
# TODO(bernhard): make sure we have documentation for these exported variables

# read out results for soil domain variables
export get_δsoil,     get_θ,     get_ψ,  get_W, get_SWATI, get_K
export get_deltasoil, get_theta, get_psi
# read out results for aboveground/scalar variables
export get_aboveground, get_δ, get_delta # get_mm


@doc raw"""
    SPAC

An instance of a soil-plant-atmopsheric continuum model.
"""
Base.@kwdef mutable struct SPAC

    # time related:
    "Reference date to relate internal use of days to real-world dates [DateTime]"
    reference_date
    "Time span of available input data [days since reference date]"
    tspan

    # atmospheric forcing:
    "DataFrame with daily atmospheric variables"          #TODO: interpolate this in time
    meteo_forcing
    "DataFrame with isotopic signatures of precipitation" #TODO: interpolate this in time
    meteo_iso_forcing
    "DataFrame with approximate typical storm duration in hours for each month"
    storm_durations

    # soil description:
    "DataFrame containing description of soil layers/horizons and soil hydraulic parameters"
    soil_horizons # a data.frame with columns: DataFrame(HorizonNr::Int, Upper_m::Float64, Lower_m::Float64, shp::soil_hydr_params)

    # vegetation characteristics:
    "DataFrame with daily vegatation variables (DENSEF, HEIGHT, LAI, SAI)" #TODO: interpolate this in time
    canopy_evolution # e.g. LAI(t), ...
    "TODO: function describing root distribution with depth (f(z_m) with z_m in meters and negative downward). Alternatively path to soil_discretization.csv" #TODO: interpolate this in time
    root_distribution # e.g. β_root

    # initial conditions:
    """
    NamedTuple: (continuousIC.soil, continuousIC.scalar), which conatains the initial
    conditions of the state variables (scalar or related to the soil).
    TODO: for soil it is either a path to soil_discretization.csv or alternatively:
        a set of functions describing initial conditions (ψM, δ18O, δ2H) with depth (f(z_m) with z_m in meters and negative downward)
    """
    continuousIC # e.g. .soil = PSIM, δ2H, δ18O, but also .aboveground = SNOW,INTR,INTS

    # simulation parameter:
    """
        params

    NamedTuple: (), which contains parameter values for the model
    """
    params # e.g. (VXYLEM_mm = 20.0, DISPERSIVITY_mm = 10, Str = "oh")
    """
        solver_options

    NamedTuple: (), which contains some solver options for the model
    """
    solver_options # reset_flag, compute_intermediate_quantities, simulate_isotopes, soil_output_depths
end

@doc raw"""
    DiscretizedSPAC

An discretization of a soil-plant-atmopsheric continuum model ready for simulation.
"""
Base.@kwdef mutable struct DiscretizedSPAC
	# input fields:
	continuous_SPAC::SPAC
	soil_discretization

	# derived fields:
	ODEProblem
	ODESolution
    ODESolution_datetime
end

# input_prefix = "isoBEAdense2010-18-reset-FALSE";
# input_path = "examples/isoBEAdense2010-18-reset-FALSE-input/";
# model = SPAC(input_path, input_prefix;
#              simulate_isotopes = simulate_isotopes);
# # # constructor of DiscretizedSPAC
# continuous_SPAC = model
# discrete_SPAC = discretize(continuous_SPAC);
# discrete_SPAC.ODEProblem.tspan
# simulate!(discrete_SPAC)
# using Plots
# plot(discrete_SPAC.ODESolution)

# Include code before defining modules
include("func_read_inputData.jl") # defines RelativeDaysFloat2DateTime which is used in module ISO

# Define modules
# on modules: https://discourse.julialang.org/t/large-programs-structuring-modules-include-such-that-to-increase-performance-and-readability/29102/5
include("module_CONSTANTS.jl");  # to bring into scope: using .CONSTANTS
include("module_KPT.jl");        using .KPT # using to bring exports into scope here
include("module_WAT.jl");        using .WAT # using to bring exports into scope
include("module_SUN.jl");        # to bring into scope: using .SUN
include("module_PET.jl");        # to bring into scope: using .PET
include("module_SNO.jl");        # to bring into scope: using .SNO
include("module_EVP.jl");        # to bring into scope: using .SNO
include("module_ISO.jl");        using .ISO

include("func_discretize_soil_domain.jl")
include("func_DiffEq_definition_u0.jl")
include("func_DiffEq_definition_p.jl")
include("func_DiffEq_definition_cb.jl")
include("func_DiffEq_definition_f.jl")
include("func_DiffEq_definition_ode.jl")
include("func_MSB_functions.jl")
include("func_postprocess.jl")

include("../examples/func_run_example.jl") # defines RelativeDaysFloat2DateTime

# struct DiscretizedSoilDomain{T <: AbstractVector}
#     # Input fields
#     # Geometry
#     "Vector of depths of upper interface [m]"
#     Upper_m::T
#     "Vector of depths of lower interface [m]"
#     Lower_m::T
#     # Soil Hydraulic properties
#     "Vector of soil hydraulic properties"
#     soil::Vector{AbstractSoilHydraulicParams} # containing p_STONEF, p_THSAT, etc...

#     # Derived fields
#     # Geometry
#     NLAYER::Int
#     p_THICK::T
#     p_SWATMAX::T # TODO(bernhard): should be replaced with θmax (i.e. actually θs)

#     # # Inner constructor:
#     # see constructor of KPT_SOILPAR_Ch1d
# end

# struct SPACSimulation
#     soil::DiscretizedSoilDomain
#     params
#     problem::ODEProblem
#     solution::ODESolution?

#     # Constructors
#     function SPACSimulation(config = "conf.toml")
#         params, u0  init_spac(config)
    #         params = init_parameters(config)
    #         u0 = init_state(config)
    #         tspan = ...

#         input("func_DiffEq_definition_cb.jl")
#         #function compute_rhs!(du, u, p, t)
#         #end
#         input("func_DiffEq_definition_cb.jl")
#         #function callbacks(du, u, p, t)
#         #end

#         problem = ODE(u0, tspan, )
#     end
# end

function LWFBrook90.discretize(continuous_SPAC::SPAC;
                                Δz = nothing,
                                tspan = nothing,
                                soil_output_depths = zeros(Float64, 0), # e.g. # soil_output_depths = [-0.35, -0.42, -0.48, -1.05]
                                Δz_functions = (rootden   = ((Δz_m)->LWFBrook90.Rootden_beta_(0.97, Δz_m = Δz_m)),
                                                PSIM_init = ((Δz_m)->fill(-6.3, length(Δz_m))),
                                                δ18Ο_init = ((Δz_m)->ifelse.(cumsum(Δz_m) .<= 0.2, -13., -10.)),
                                                δ2Η_init  = ((Δz_m)->ifelse.(cumsum(Δz_m) .<= 0.2, -95., -70.)))
        )
    soil_horizons = copy(continuous_SPAC.soil_horizons)
    soil_horizon_extent = (maximum(soil_horizons.Upper_m), minimum(soil_horizons.Lower_m))

    ##########
    # Discretize the model in space as `soil_discretization`
    ### returns `soil_discretization` a DataFrame with columns: (Upper_m, Lower_m)
    ### returns `soil_discretization` a DataFrame with columns: (Rootden_)
    ### returns `soil_discretization` a DataFrame with columns: (uAux_PSIM_init_kPa, u_delta18O_init_permil, u_delta2H_init_permil)

    # Define grid for spatial discretization
    # Note that later on grid will be further refined with observation nodes and required interfaces (if needed)
    if isnothing(Δz)
         # a) either read the discretization from a file `soil_discretization.csv`
         #    including Rootden_, and initial PSIM, delta18O, and delta2H
        use_soil_discretization_csv = true
        path_soil_discretization    = continuous_SPAC.continuousIC.soil
        soil_discretization         = LWFBrook90.read_path_soil_discretization(path_soil_discretization)

        # Assert type stability by executing disallowmissing!
        if (continuous_SPAC.solver_options.simulate_isotopes)
            disallowmissing!(soil_discretization, [:Rootden_, :uAux_PSIM_init_kPa, :u_delta18O_init_permil, :u_delta2H_init_permil])
        else
            disallowmissing!(soil_discretization, [:Rootden_, :uAux_PSIM_init_kPa])
        end

    else
        # b) or use manually defined Δz
        #    and require functions for Rootden_, and initial PSIM, delta18O, and delta2H

        @assert all(Δz .> 0)

        use_soil_discretization_csv = false
        interfaces_m = [0, cumsum(Δz)...]
        soil_discretization = DataFrame(
            Upper_m = -interfaces_m[Not(end)],
            Lower_m = -interfaces_m[Not(1)],
            Rootden_ = NaN,
            uAux_PSIM_init_kPa = NaN,
            u_delta18O_init_permil = NaN,
            u_delta2H_init_permil = NaN)

        # check that functions are provided
        @assert keys(Δz_functions) == (:rootden, :PSIM_init, :δ18Ο_init, :δ2Η_init)
        @assert Δz_functions.rootden   isa Function
        @assert Δz_functions.PSIM_init isa Function
        @assert Δz_functions.δ18Ο_init isa Function
        @assert Δz_functions.δ2Η_init  isa Function
    end

    # Check validity of loaded soil discretization
    @assert soil_horizon_extent[1] ≈ soil_discretization.Upper_m[1] "Spatial domain of soil discretization is not compatible with provided soil horizons"
    if soil_horizon_extent[2] ≈ soil_discretization.Lower_m[end]
        # all okay no need to extend z horizon
    elseif soil_horizon_extent[2] < soil_discretization.Lower_m[end]
        @warn """
        Spatial domain soil discretization is smaller than the provided soil horizon information ($soil_horizon_extent m)
        (lower end of requested soil discr.: $(soil_discretization.Lower_m[end]))
        """
    elseif soil_horizon_extent[2] > soil_discretization.Lower_m[end]
        @warn """
        Spatial domain soil discretization is larger than the provided soil horizon information. Lowest soil horizon will be extended.

        Lowest soil horizon from $(soil_horizons[end,:Upper_m])m to $(soil_horizons[end,:Lower_m])m will be extended down to $(soil_discretization.Lower_m[end])m.
        """
        push!(soil_horizons, soil_horizons[end,:])
        soil_horizons[end,:Lower_m] = soil_discretization.Lower_m[end]
        soil_horizons[end,:Upper_m] = soil_horizons[end-1,:Lower_m]
        soil_horizons[end,:HorizonNr] = soil_horizons[end-1,:HorizonNr] + 1
    end

    # Define initial conditions and root densities
    if use_soil_discretization_csv
        # Nothing needed as columns in `soil_discretization` are exptected to be already filled.

        @assert !(any(isnan.(soil_discretization.Rootden_)))           "Error: 'soil_discretization.Rootden_' contains NaNs."
        @assert !(any(isnan.(soil_discretization.uAux_PSIM_init_kPa))) "Error: 'soil_discretization.uAux_PSIM_init_kPa' contains NaNs."
        @assert typeof(continuous_SPAC.continuousIC.soil) == String && isfile(continuous_SPAC.continuousIC.soil) "Expected 'continuous_SPAC.continuousIC.soil' to be a file containing the initial conditions."
        @assert typeof(continuous_SPAC.root_distribution) == String && isfile(continuous_SPAC.root_distribution) "Expected 'continuous_SPAC.continuousIC.soil' to be a file containing the initial conditions."
        @assert continuous_SPAC.continuousIC.soil == continuous_SPAC.root_distribution "If provided as .csv file, initial conditions and root distribution should refer to the same file."
    else
        soil_discretization = discretize_soil(;
            Δz_m = Δz,
            Rootden_               = Δz_functions.rootden,
            uAux_PSIM_init_kPa     = Δz_functions.PSIM_init,
            u_delta18O_init_permil = Δz_functions.δ18Ο_init,
            u_delta2H_init_permil  = Δz_functions.δ2Η_init)
    end

    # TODO(bernhard): make above code a four step procedure:
        # 1) define a grid resolution Δz (either a) reading in from soil_discretization or b) manually defined vector)
        # 2) get info about initial conditions (either a) reading in from soil_discretization or then or b) manually defined mathematical function  )
        # 2) get info about root distribution (either a) reading in from soil_discretization or then or b) manually defined mathematical function  )
        # 3) map initial condition to discretized grid resolution Δz
        # 3) map root distribution to discretized grid resolution Δz (and interpolate in time giving us a function rootden(z, t))
        # 4) include additional interfaces for observation_nodes and other needed interfaces (such as soil layers)...
            # TODO(bernhard): combine step 4) with 1). -> Requires mapping in 3) regardless of whether we have a) or b)
        # 5) map soil physical properties from soil layers to discretized grid resolution Δz
    ####################

    ####################
    ## Discretize soil parameters and interpolate discretized root distribution
    # Define refinement of grid with soil_output_depths
    soil_params =
        LWFBrook90.refine_soil_discretization(
            soil_horizons,
            soil_discretization,
            soil_output_depths,
            continuous_SPAC.params[:IDEPTH_m], # :IDEPTH_m is unused later on
            continuous_SPAC.params[:QDEPTH_m]) # :QDEPTH_m is unused later on
    refined_soil_discretization = soil_params["refined_soil_discretization"]

    # Interpolate discretized root distribution in time
    p_fT_RELDEN = LWFBrook90.HammelKennel_transient_root_density(;
        timepoints         = continuous_SPAC.meteo_forcing.p_days,
        p_AGE              = continuous_SPAC.canopy_evolution.p_AGE,
        p_INITRDEP         = continuous_SPAC.params[:INITRDEP],
        p_INITRLEN         = continuous_SPAC.params[:INITRLEN],
        p_RGROPER_y        = continuous_SPAC.params[:RGROPER],
        p_RGRORATE_m_per_y = continuous_SPAC.params[:RGRORATE],
        p_THICK               = soil_params["THICK"],
        final_Rootden_profile = soil_params["final_Rootden_"]);
    # TODO(bernhard): document input parameters: INITRDEP, INITRLEN, RGROPER, tini, frelden, MAXLAI, HEIGHT_baseline_m
    # TOOD(bernhard): remove from params: IDEPTH_m, QDEPTH_m, INITRDEP, RGRORATE, INITRDEP, INITRLEN, RGROPER
    # display(heatmap(p_fT_RELDEN', ylabel = "SOIL LAYER", xlabel = "Time (days)", yflip=true, colorbar_title = "Root density"))
    ####################

    ####################
    # Define parameters for differential equation
    p = define_LWFB90_p(continuous_SPAC, soil_params, p_fT_RELDEN)
    # using Plots
    # hline([0; cumsum(p.p_THICK)], yflip = true, xticks = false,
    #     title = "N_layer = "*string(p.NLAYER))
   ####################

    ####################
    # Define state vector u for DiffEq.jl
    # a) allocation of u0
    u0 = define_LWFB90_u0(;simulate_isotopes = continuous_SPAC.solver_options.simulate_isotopes,
                          compute_intermediate_quantities = continuous_SPAC.solver_options.compute_intermediate_quantities,
                          NLAYER = soil_params["NLAYER"])
    ####################

    ####################
    # Define initial states of differential equation
    # state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI
    # Create u0 for DiffEq.jl
    # b) initialization of u0
    init_LWFB90_u0!(;u0=u0, continuous_SPAC=continuous_SPAC, soil_params=soil_params, p_soil=p.p_soil)
    ####################

    ####################
    # Define ODE problem which consists of
    #   - definition of right-hand-side (RHS) function f
    #   - definition of callback function cb
    #   - u0:     initial condition of states
    #   - tspan:  definition of model time span
    #   - p:      parameters

    # Define simulation time span:
    # tspan = (0.,  5.) # simulate 5 days
    # tspan = (0.,  100.) # simulate 100 days # NOTE: KAU bugs in "branch 005-" when at least 3*365
    # tspan = (minimum(continuous_SPAC.meteo_forcing[:,"days"]),
    #         maximum(continuous_SPAC.meteo_forcing[:,"days"])) # simulate all available days
    # tspan = (LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1980,1,1), reference_date),
    #          LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1985,1,1), reference_date)) # simulates selected period

    if isnothing(tspan)
        tspan = continuous_SPAC.tspan
    else
        @warn "Overwriting tspan defined in SPAC $(continuous_SPAC.tspan) with provided value of $tspan"
        tspan = tspan
    end


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
    ####################

    return DiscretizedSPAC(;
        continuous_SPAC     = continuous_SPAC,
        soil_discretization = refined_soil_discretization,
        ODEProblem          = ode_LWFBrook90,
        ODESolution         = nothing,
        ODESolution_datetime= nothing)
end

function simulate!(s::DiscretizedSPAC)
    sol_SPAC = solve_LWFB90(s.ODEProblem);
    s.ODESolution = sol_SPAC;

    @assert (SciMLBase.successful_retcode(sol_SPAC)) "Problem with simulation: Return code of simulation was '$(sol_SPAC.retcode)'"

    # also save datetimes
    s.ODESolution_datetime = LWFBrook90.RelativeDaysFloat2DateTime.(s.ODESolution.t, s.continuous_SPAC.reference_date)
    return nothing
end

############################################################################################
############################################################################################
############################################################################################

@doc raw"""
    run_simulation(args)

Runs a simulation defined by input files within a folder and returns solved DiscretizedSPAC.

The function run_simulation() takes as single argument a vector of two or three strings defining
the input_path and input_prefix of a series of input definition files (and "true"/"false" whether
to simulate isotopes).
The function loads these files, runs the simulation and returns the solution object and input arguments
"""
function run_simulation(args)
    # args = ("test-assets/Hammel-2001/input-files", "Hammel_loam-NLayer-27-RESET=TRUE")
    # args = ["examples/isoBEAdense2010-18-reset-FALSE-input/" "isoBEAdense2010-18-reset-FALSE" "true"]
    # args = ["test-assets/Hammel-2001/input-files-ISO" "Hammel_sand-NLayer-27-RESET=FALSE" "true"]
    # args = ["test-assets/Hammel-2001/input-files-ISO" "Hammel_loam-NLayer-27-RESET=FALSE" "true"]
    # args = ["test/test-assets/Hammel-2001/input-files-ISO" "Hammel_loam-NLayer-27-RESET=FALSE" "false"]
    # args = ["test-assets/Hammel-2001/input-files-ISO" "Hammel_loam-NLayer-27-RESET=FALSE" "true"]
    # args = ["test-assets/Hammel-2001/input-files" "Hammel_sand-NLayer-27-RESET=FALSE" "false"]
    # args = ["test-assets/BEA-2016/input-files/"  "BEA2016-reset-FALSE"  "false"]
    # args = ["../examples/isoBEAdense2010-18-reset-FALSE-input/" "isoBEAdense2010-18-reset-FALSE" "true"]
    # args = ["test-assets/DAV-2020/input-files/"  "DAV_LW1_def"  "false"]
    # @show now()
    # @show args

    @assert length(args) >= 2

    input_path = args[1]
    input_prefix = args[2]
    if (length(args) == 3)
        simulate_isotopes = args[3] == "true"
    else
        simulate_isotopes = false
    end
    # @show simulate_isotopes

    ####################
    # Define simulation model by reading in system definition and input data
    model = SPAC(input_path, input_prefix;
                      simulate_isotopes = simulate_isotopes);
    ####################

    ####################
    # Prepare simulation by discretizing spatial domain
    simulation = LWFBrook90.discretize(model);

    # Solve ODE:
    LWFBrook90.simulate!(simulation)
    # plot(simulation.ODESolution)
    sol_LWFBrook90 = simulation.ODESolution

    # simulation.solver_options.simulate_isotopes
    ####################

    # using Plots
    # using Interpolations: interpolate, extrapolate, NoInterp, Gridded, Constant, Next, Previous, Flat, Throw, scale
    # scatter(input_meteoveg.days[1:10], input_meteoveg.PRECIN[1:10])
    # t_to_eval = 0:0.2:10
    # plot!(t_to_eval, p.p_PREC(t_to_eval), xtick = 1:10)
    # plot!(t_to_eval,extrapolate(interpolate((input_meteoveg.days .- 0.00001, ), input_meteoveg.PRECIN, Gridded(Constant{Previous}())), Throw())(t_to_eval))
    # sol_LWFBrook90 = solve_LWFB90(u0, tspan, p)
    # using Plots, Measures
    # optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
    # pl2 = LWFBrook90.ISO.plotisotopes(
    #     sol_LWFBrook90, optim_ticks;
    #     layout = grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4]),
    #     size=(1000,1400), dpi = 300, leftmargin = 15mm);
    # plot!(pl2, link = :x)

    @info sol_LWFBrook90.destats
    @info "Time steps for solving: $(sol_LWFBrook90.destats.naccept) ($(sol_LWFBrook90.destats.naccept) accepted out of $(sol_LWFBrook90.destats.nreject + sol_LWFBrook90.destats.naccept) total)"
    @show now()

    return (simulation, input_prefix, input_path)
end

end # module
