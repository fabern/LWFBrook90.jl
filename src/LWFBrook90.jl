module LWFBrook90

using OrdinaryDiffEq  # instead of loading the full DifferentialEquations
using DiffEqCallbacks # instead of loading the full DifferentialEquations
using RecipesBase
using LinearAlgebra
using StatsBase: mean, weights
using ComponentArrays
using UnPack: @unpack

using Dates: now
# using Infiltrator

export SPAC, DiscretizedSPAC, discretize, simulate, simulate!
export read_inputData
export discretize_soil, Rootden_beta_
export define_LWFB90_p, define_LWFB90_u0, solve_LWFB90
export KPT_SOILPAR_Mvg1d, KPT_SOILPAR_Ch1d
export RelativeDaysFloat2DateTime, plot_LWFBrook90

export run_simulation, plot_and_save_results, find_indices
export get_auxiliary_variables, get_θ, get_δ, get_ψ, get_δsoil, get_aboveground

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
include("module_KPT.jl");        using .KPT # using to bring exports into scope
include("module_WAT.jl");        using .WAT # using to bring exports into scope
include("module_SUN.jl");        # to bring into scope: using .SUN
include("module_PET.jl");        # to bring into scope: using .PET
include("module_SNO.jl");        # to bring into scope: using .SNO
include("module_EVP.jl");        # to bring into scope: using .SNO
include("module_ISO.jl");        # to bring into scope: using .ISO

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

# struct SoilState
#     head
#     c18o
#     c2h
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
# run(sim::SPACsimulation)
#    sim.solution = solve(sim, )
# end

function LWFBrook90.discretize(continuous_SPAC::SPAC;
                                Δz = nothing, tspan = nothing,
                                soil_output_depths = zeros(Float64, 0) # e.g. # soil_output_depths = [-0.35, -0.42, -0.48, -1.05]
        )
    ##########
    # Discretize the model in space
    ### returns DataFrame with columns: (Upper_m, Lower_m)
    ### returns DataFrame with columns: (Rootden_)
    ### returns DataFrame with columns: (uAux_PSIM_init_kPa, u_delta18O_init_permil, u_delta2H_init_permil)

    # Needed extents of domain
    zspan = (maximum(continuous_SPAC.soil_horizons.Upper_m),
             minimum(continuous_SPAC.soil_horizons.Lower_m))

    # Define grid for spatial discretization
    # Note that later on grid will be further refined with observation nodes and required interfaces (if needed)
    if isnothing(Δz)
         # a) either read the discretization from a file `soil_discretization.csv`
        use_soil_discretization_csv = true
        path_soil_discretization = continuous_SPAC.continuousIC.soil
        soil_discretization = LWFBrook90.read_path_soil_discretization(path_soil_discretization)
    else
        # b) or use manually defined Δz
        use_soil_discretization_csv = false
        interfaces_m = [0, cumsum(Δz)...]
        soil_discretization = DataFrame(
            Upper_m = -interfaces_m[Not(end)],
            Lower_m = -interfaces_m[Not(1)],
            Rootden_ = NaN,
            uAux_PSIM_init_kPa = NaN,
            u_delta18O_init_permil = NaN,
            u_delta2H_init_permil = NaN)
    end
    @assert zspan[1] ≈ soil_discretization.Upper_m[1] "Soil discretization is not compatible with provided soil horizons"
    @assert zspan[2] ≈ soil_discretization.Lower_m[end] "Soil discretization is not compatible with provided soil horizons"

    # Define initial conditions and root densities
    if use_soil_discretization_csv
        # Nothing needed as columns in `soil_discretization` are exptected to be already filled.

        # Just checks:
        @assert !(any(isnan.(soil_discretization.Rootden_))) "Error: 'soil_discretization.Rootden_' contains NaNs."
        @assert !(any(isnan.(soil_discretization.uAux_PSIM_init_kPa))) "Error: 'soil_discretization.uAux_PSIM_init_kPa' contains NaNs."

        # Just assert type stability by executing disallowmissing!
        if (continuous_SPAC.solver_options.simulate_isotopes)
            disallowmissing!(soil_discretization, [:Rootden_, :uAux_PSIM_init_kPa, :u_delta18O_init_permil, :u_delta2H_init_permil])
        else
            disallowmissing!(soil_discretization, [:Rootden_, :uAux_PSIM_init_kPa])
        end

        # TODO:
        @assert typeof(continuous_SPAC.continuousIC.soil) == String && isfile(continuous_SPAC.continuousIC.soil) "Expected 'continuous_SPAC.continuousIC.soil' to be a file containing the initial conditions."
        @assert typeof(continuous_SPAC.root_distribution) == String && isfile(continuous_SPAC.root_distribution) "Expected 'continuous_SPAC.continuousIC.soil' to be a file containing the initial conditions."
        @assert continuous_SPAC.continuousIC.soil == continuous_SPAC.root_distribution
        # TODO: assert that root distribution is consistent, otherwise it needs to be discretized and mapped too....
    else
        @warn "Using hardcoded values for soil_discretization.... TODO(improve interface to allow to modify these hardcoded values)"
        f1 = (Δz_m) -> LWFBrook90.Rootden_beta_(0.97, Δz_m = Δz_m)  # function for root density as f(Δz)
        f2 = (Δz_m) -> fill(-6.3, length(Δz_m))          # function for initial conditions as f(Δz)
        f3 = ifelse.(cumsum(Δz) .<= 0.2, -13., -10.)     # initial conditions as vector
        f4 = ifelse.(cumsum(Δz) .<= 0.2, -95., -70.)     # initial conditions as vector

        # continuous_SPAC.root_distribution = f1
        # continuous_SPAC.continuousIC.soil = f2

        # @assert continuous_SPAC.continuousIC.soil isa Function
        # @assert continuous_SPAC.root_distribution isa Function
        # @error "Currently manual definition of Δz is unsupported."

        soil_discretization = discretize_soil(;
            Δz_m = Δz,
            Rootden_ = f1,
            uAux_PSIM_init_kPa = f2,
            u_delta18O_init_permil = f3,
            u_delta2H_init_permil  = f4)

        # # TO DEVELOP: possibly values of Δz_m:
        # ## for Δz_m in (
        # ##     [0.04, 0.04, 0.12, 0.25, 0.3, 0.35, 0.1],
        # ##     [fill(0.04, 5); fill(0.05, 5); 0.3;            0.35;         0.1], # N=13
        # ##     [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 5); 0.1], # N=21
        # ##     [fill(0.02, 10); fill(0.025, 10); fill(0.03, 10); fill(0.035, 10); 0.1], # N=41
        # ##     [fill(0.02, 60)... ], # N=60, 2cm similar to Pollacco et al. 2022 suggestions
        # ##     )
        # # As a test: subsequently increase resolution of the top layers.
        # Δz_m = [0.04, 0.04, 0.12, 0.25, 0.3, 0.35, 0.1];                           # grid spacing (heterogenous), meter (N=7)
        # # Δz_m = [0.04, 0.04, 0.04, 0.08, 0.25, 0.3, 0.35, 0.1]                    # grid spacing (heterogenous), meter (N=8)
        # # Δz_m = [0.04, 0.04, 0.04, 0.04, 0.04, 0.25, 0.3, 0.35, 0.1]              # grid spacing (heterogenous), meter (N=9)
        # # Δz_m = [fill(0.04, 5); fill(0.05, 5); 0.3; 0.35; 0.1]                    # grid spacing (heterogenous), meter (N=13)
        # # Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); 0.35; 0.1]          # grid spacing (heterogenous), meter (N=17)
        # Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 5); 0.1]; # grid spacing (heterogenous), meter (N=21)
        # if zspan[2] == -1.1
        #     Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 5)]; # grid spacing (heterogenous), meter (N=20)
        # end
    end

    # TODO(bernhard): make above code a three step procedure:
        # 1) define a grid resolution Δz (either a) reading in from soil_discretization or b) manually defined vector)
        # 2) get info about initial conditions (either a) reading in from soil_discretization or then or b) manually defined mathematical function  )
        # 2) get info about root distribution (either a) reading in from soil_discretization or then or b) manually defined mathematical function  )
        # 3) map initial condition to discretized grid resolution Δz
        # 3) map root distribution to discretized grid resolution Δz (and interpolate in time giving us a function rootden(z, t))
        # 4) include additional interfaces for observation_nodes...
            # TODO(bernhard): combine step 4) with 1). -> Requires mapping in 3) regardless of whether we have a) or b)
    ####################

    ####################
    ## Discretize soil parameters and interpolate discretized root distribution
    # Define refinement of grid with soil_output_depths
    soil_discr =
        LWFBrook90.refine_soil_discretization(
            continuous_SPAC.soil_horizons,
            soil_discretization,
            soil_output_depths,
            continuous_SPAC.params[:IDEPTH_m], # :IDEPTH_m is unused later on
            continuous_SPAC.params[:QDEPTH_m], # :QDEPTH_m is unused later on
            continuous_SPAC.params[:INITRDEP], # :INITRDEP is unused later on
            continuous_SPAC.params[:RGRORATE]) # :RGRORATE is unused later on
    ####################

    ####################
    # Define parameters for differential equation
    p = define_LWFB90_p(continuous_SPAC, soil_discr)
    # using Plots
    # hline([0; cumsum(p.p_THICK)], yflip = true, xticks = false,
    #     title = "N_layer = "*string(p.NLAYER))
   ####################

    ####################
    # Define state vector u for DiffEq.jl
    # a) allocation of u0
    u0 = define_LWFB90_u0(;simulate_isotopes = continuous_SPAC.solver_options.simulate_isotopes,
                          compute_intermediate_quantities = continuous_SPAC.solver_options.compute_intermediate_quantities,
                          NLAYER = soil_discr["NLAYER"])
    ####################

    ####################
    # Define initial states of differential equation
    # state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI
    # Create u0 for DiffEq.jl
    # b) initialization of u0
    init_LWFB90_u0!(;u0=u0, continuous_SPAC=continuous_SPAC, soil_discr=soil_discr, p_soil=p.p_soil)
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
        soil_discretization = soil_discretization,
        ODEProblem          = ode_LWFBrook90,
        ODESolution         = nothing)
end

function simulate!(s::DiscretizedSPAC)
    sol_SPAC = solve_LWFB90(s.ODEProblem);
    s.ODESolution = sol_SPAC;
    # return sol_SPAC
    return nothing
end

############################################################################################
############################################################################################
############################################################################################

@doc raw"""
    run_simulation(args)

Runs a simulation defined by input files within a folder and returns solution object.

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
    soil_discretization = discretize_soil(model.continuousIC.soil)
    zspan = (maximum(model.soil_horizons.Upper_m),
             minimum(model.soil_horizons.Lower_m))
    @assert zspan[1] == soil_discretization.Upper_m[1] "Soil discretization is not compatible with provided soil horizons"
    @assert zspan[2] == soil_discretization.Lower_m[end] "Soil discretization is not compatible with provided soil horizons"

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
    @info sol_LWFBrook90.destats
    @info "Time steps for solving: $(sol_LWFBrook90.destats.naccept) ($(sol_LWFBrook90.destats.naccept) accepted out of $(sol_LWFBrook90.destats.nreject + sol_LWFBrook90.destats.naccept) total)"
    # using Plots, Measures
    # optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
    # pl2 = LWFBrook90.ISO.plotisotopes(
    #     sol_LWFBrook90, optim_ticks;
    #     layout = grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4]),
    #     size=(1000,1400), dpi = 300, leftmargin = 15mm);
    # plot!(pl2, link = :x)
    @show now()

    return (sol_LWFBrook90, input_prefix, input_path)
end

end # module
