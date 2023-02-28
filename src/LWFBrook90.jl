module LWFBrook90

using SciMLBase       # instead of loading the full DifferentialEquations
using OrdinaryDiffEq  # instead of loading the full DifferentialEquations
using DiffEqCallbacks # instead of loading the full DifferentialEquations
using RecipesBase, PlotUtils, Measures
using LinearAlgebra
using StatsBase: mean, weights
using ComponentArrays
using UnPack: @unpack
using Dates: now, Day, dayofyear
using Printf: @sprintf
using Interpolations: interpolate, extrapolate, NoInterp, Gridded, Constant, Next, Previous, Flat, Throw, scale, BSpline, linear_interpolation

# TODO(bernhard): make sure we have documentation for these exported variables
export SPAC, DiscretizedSPAC, loadSPAC, setup, simulate!
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

An instance of a soil-plant-atmopsheric continuum model with the following fields:

- `reference_date`: DateTime to relate internal use of numerical days to real-world dates
- `tspan`: Tuple `(0, Int)` Time span of available input data in days since reference date
- `solver_options`: `NamedTuple`: (compute_intermediate_quantities, simulate_isotopes, DTIMAX, DSWMAX, DPSIMAX), containing some solver options for the model

- `soil_discretization`: `DataFrame` with contents from `soil_discretizations.csv`, i.e. containing columns:
        `Upper_m`,`Lower_m`,`Rootden_`,`uAux_PSIM_init_kPa`,`u_delta18O_init_permil`,`u_delta2H_init_permil`

        Either loaded from `soil_discretizations.csv` or then specified as argument to loadSPAC().
        Unused values are NaN.

- `forcing`: `NamedTuple`: (:meteo, :meteo\_iso, :storm\_durations), containing atmospheric forcings
    - `forcing.meteo`: DataFrame with daily atmospheric variables
    - `forcing.meteo_iso`: DataFrame with isotopic signatures of precipitation
    - `forcing.storm_durations`: DataFrame with approximate typical storm duration in hours for each month
-   `pars`: `NamedTuple`: (:params, )
    - `pars.params`: `NamedTuple`: (), containing scalar parameter values for the model
    - `pars.root_distribution`: either:
        - `NamedTuple`: (beta = 0.97, z\_rootMax\_m = nothing) parametrization of root distribution with depth (f(z\_m) with z\_m in meters and negative downward). Alternatively path to soil\_discretization.csv"
        - or String: `"soil_discretization.csv"` meaning that it must be defined in `soil_discretizations.csv`
    - `pars.IC_scalar`: `NamedTuple`: (), containing initial conditions of the scalar state variables
    - `pars.IC_soil`: initial conditions of the state variables (scalar or related to the soil). Either:
        - `NamedTuple`: (PSIM\_init, δ18O\_init, δ2H\_init), containing constant values for the initial conditions
        - [NOT IMPELMENTED:]`NamedTuple`: (PSIM\_init, δ18O\_init, δ2H\_init), containing functions of the initial values with argument Δz
        - or String: `"soil_discretization.csv"` meaning that it must be defined in `soil_discretizations.csv`
    - `pars.canopy_evolution`: canopy parameters (LAI, SAI, DENSEF, HEIGHT). Either:
        - `NamedTuple`: (DENSEF = 100, HEIGHT = 25, SAI = 100, LAI = (DOY\_Bstart, Bduration, DOY\_Cstart, Cduration, LAI\_perc\_BtoC, LAI\_perc\_CtoB)), containing constant values and parameters for LAI\_relative interpolation
        - or DataFrame: with daily values
    - `pars.soil_horizons`: DataFrame containing description of soil layers/horizons and soil hydraulic parameters
"""
Base.@kwdef mutable struct SPAC
    # A) Assumed known: will not be estimated
    # A1) Simulation related (non-physical)
    reference_date
    tspan
    solver_options
    soil_discretization
    # A2) Physical forcing
    forcing # atmospheric forcing:  #TODO: interpolate forcing.meteo and forcing.meteo_iso in time
    # B) Model parameters: might be estimated
    pars
        # pars.params:            e.g. (VXYLEM_mm = 20.0, DISPERSIVITY_mm = 10, Str = "oh")
        # pars.root_distribution: e.g. (beta = 0.97, z_rootMax_m = nothing)
        # pars.IC:                e.g. .soil = PSIM, δ2H, δ18O, but also .scalar = SNOW,INTR,INTS
        # pars.root_distribution = (beta = 0.97, z_rootMax_m = nothing)
        # pars.canopy_evolution:  e.g. LAI(t), ...
        # pars.soil_horizons
end

@doc raw"""
    DiscretizedSPAC

An discretization of a soil-plant-atmopsheric continuum model ready for simulation.
"""
Base.@kwdef mutable struct DiscretizedSPAC
	# input fields:
	parametrizedSPAC::SPAC

	# derived fields:
	ODEProblem
	ODESolution
    ODESolution_datetime
end

# input_prefix = "isoBEAdense2010-18-reset-FALSE";
# input_path = "examples/isoBEAdense2010-18-reset-FALSE-input/";
# model = loadSPAC(input_path, input_prefix;
#              simulate_isotopes = simulate_isotopes);
# # # constructor of DiscretizedSPAC
# parametrizedSPAC = model
# discrete_SPAC = setup(parametrizedSPAC);
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
#         params, u0  init_loadSPAC(config)
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



"""
    function setup(parametrizedSPAC::SPAC;
                   requested_tspan = nothing,
                   soil_output_depths_m::Vector = zeros(Float64, 0))

Takes the definition of SPAC model and discretize to a system of ODEs that can be solved by
the package DifferentialEquations.jl

If needed, the computational grid of the soil is refined to output values at specific depths,
e.g. by doing `setup(SPAC; soil_output_depths_m = [-0.35, -0.42, -0.48, -1.05])`.

The argument `requested_tspan` is a tuple defining the duration of the simulation with a start-
and an end-day relative to the reference date: `setup(SPAC; requested_tspan = (0,150))`.
Note that the reference date given by the `SPAC.reference_date`.
"""
function setup(parametrizedSPAC::SPAC;
                requested_tspan = nothing,
                soil_output_depths_m::Vector = zeros(Float64, 0)) # e.g. # soil_output_depths_m = [-0.35, -0.42, -0.48, -1.05]

    @assert all(soil_output_depths_m .< 0)

    # This function prepares a discretizedSPAC, which is a container for a DifferentialEquations::ODEProblem.
    # A discretizedSPAC stores:
        # Its needed input:
            # - SPAC (that contains arguments `requested_tspan` and `soil_output_depths_m`)
            #        (hence setup(discreteSPAC.parametrizedSPAC) works as expected)
        # And derived fields:
            # - ODEProblem
            # - ODESolution
            # - ODESolution_datetime

    # To prepare the ODE we do:
        # make time dependent parameter
        # refine the soil discretization as needed for soil_output_depths_m or infiltration depths...
        # modify the simulation time span `tspan`

    modifiedSPAC = deepcopy(parametrizedSPAC); # make a copy to put into DiscretizedSPAC

    ##########
    # a) Refine soil disretization to provide all needed output
    # Define grid for spatial discretization
    # refined_Δz, IDEPTH_idx, QDEPTH_idx =
    #     LWFBrook90.refine_soil_discretization(
    #         Δz,
    #         soil_output_depths_m,
    #         modifiedSPAC.pars.params[:IDEPTH_m],
    #         modifiedSPAC.pars.params[:QDEPTH_m])
    refined_soil_discretizationDF, IDEPTH_idx, QDEPTH_idx =
        LWFBrook90.refine_soil_discretization(
            # modifiedSPAC.soil_discretization.Δz,
            modifiedSPAC.soil_discretization.df,
            modifiedSPAC.pars.soil_horizons,
            soil_output_depths_m,
            modifiedSPAC.pars.params[:IDEPTH_m],
            modifiedSPAC.pars.params[:QDEPTH_m])

    # Discretize the model in space as `soil_discretization`
    final_soil_discretizationDF = map_soil_horizons_to_discretization(modifiedSPAC.pars.soil_horizons, refined_soil_discretizationDF)#computational_grid)

    # Update soil_discretization in underlying SPAC model
    modifiedSPAC.soil_discretization = (
        Δz = final_soil_discretizationDF.Upper_m - final_soil_discretizationDF.Lower_m,
        df = final_soil_discretizationDF)

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
    ## c) Derive time evolution of aboveground vegetation based on parameter from SPAC
    canopy_evolution_relative = generate_canopy_timeseries_relative(
        modifiedSPAC.pars.canopy_evolution,
        days = modifiedSPAC.forcing.meteo["p_days"],
        reference_date = modifiedSPAC.reference_date)
    canopy_evolutionDF = make_absolute_from_relative(
                aboveground_relative          = canopy_evolution_relative,
                p_MAXLAI                      = modifiedSPAC.pars.params[:MAXLAI],
                p_SAI_baseline_               = modifiedSPAC.pars.params[:SAI_baseline_],
                p_DENSEF_baseline_            = modifiedSPAC.pars.params[:DENSEF_baseline_],
                p_AGE_baseline_yrs            = modifiedSPAC.pars.params[:AGE_baseline_yrs],
                p_HEIGHT_baseline_m           = modifiedSPAC.pars.params[:HEIGHT_baseline_m])
    ####################

    ####################
    ## d) Interpolate vegetation parameter in time for use as parameters
    # Aboveground: LAI, SAI, DENSEF, HEIGHT, AGE
    vegetation_fT = interpolate_aboveground_veg(canopy_evolutionDF.AboveGround)
    ## Interpolate discretized root distribution in time
        # b) Make root growth module on final discretized soil...
    vegetation_fT["p_RELDEN"] = LWFBrook90.HammelKennel_transient_root_density(;
        timepoints         = modifiedSPAC.forcing.meteo["p_days"],
        AGE_at_timepoints  = vegetation_fT["p_AGE"].(modifiedSPAC.forcing.meteo["p_days"]),
        p_INITRDEP         = modifiedSPAC.pars.params[:INITRDEP],
        p_INITRLEN         = modifiedSPAC.pars.params[:INITRLEN],
        p_RGROPER_y        = modifiedSPAC.pars.params[:RGROPER],
        p_RGRORATE_m_per_y = modifiedSPAC.pars.params[:RGRORATE],
        p_THICK               = 1000*modifiedSPAC.soil_discretization.Δz,
        final_Rootden_profile = modifiedSPAC.soil_discretization.df.Rootden_);
    # TODO(bernhard): document input parameters: INITRDEP, INITRLEN, RGROPER, tini, frelden, MAXLAI, HEIGHT_baseline_m
    # TOOD(bernhard): remove from params: IDEPTH_m, QDEPTH_m, INITRDEP, RGRORATE, INITRDEP, INITRLEN, RGROPER
    # display(heatmap(vegetation_fT["p_RELDEN"]', ylabel = "SOIL LAYER", xlabel = "Time (days)", yflip=true, colorbar_title = "Root density"))
    ####################

    ####################
    # Define parameters for differential equation
    p = define_LWFB90_p(modifiedSPAC, vegetation_fT, IDEPTH_idx, QDEPTH_idx)
    # using Plots
    # hline([0; cumsum(p.p_soil.p_THICK)], yflip = true, xticks = false,
    #     title = "N_layer = "*string(p.NLAYER))
   ####################

    ####################
    # Define state vector u for DiffEq.jl and initial states u0
        # state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI
    # a) allocation of u0
    u0 = define_LWFB90_u0(;simulate_isotopes = modifiedSPAC.solver_options.simulate_isotopes,
                          compute_intermediate_quantities = modifiedSPAC.solver_options.compute_intermediate_quantities,
                          NLAYER = nrow(modifiedSPAC.soil_discretization.df))
    # b) initialization of u0
    init_LWFB90_u0!(;u0=u0, parametrizedSPAC=modifiedSPAC, p_soil=p.p_soil)
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
    # tspan = (minimum(modifiedSPAC.meteo_forcing[:,"days"]),
    #         maximum(modifiedSPAC.meteo_forcing[:,"days"])) # simulate all available days
    # tspan = (LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1980,1,1), reference_date),
    #          LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1985,1,1), reference_date)) # simulates selected period

    if isnothing(requested_tspan)
        tspan_to_use = modifiedSPAC.tspan
    else
        @warn "Overwriting tspan defined in SPAC $(modifiedSPAC.tspan) with provided value of $tspan"
        tspan_to_use = requested_tspan

        # Update tspan in underlying SPAC model
        modifiedSPAC.tspan = requested_tspan
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
                    tspan_to_use,
                    p,
                    callback=cb_func)
    ####################

    return DiscretizedSPAC(;
        parametrizedSPAC    = modifiedSPAC,
        ODEProblem          = ode_LWFBrook90,
        ODESolution         = nothing,
        ODESolution_datetime= nothing)
end

function Base.show(io::IO, mime::MIME"text/plain", discSPAC::DiscretizedSPAC)
    println(io, "Discretized SPAC model: =============== solution was computed: $(!isnothing(discSPAC.ODESolution))")
    Base.show(io, mime, discSPAC.parametrizedSPAC; show_SPAC_title=false)

    # DO WE NEED BELOW EXPLICIT CANOPY EVOLUTION? NO.
    # println(io, "\n===== CANOPY EVOLUTION:===============")
    # println(io, show_avg_and_range(model.canopy_evolution.p_AGE.(model.tspan),    "AGE         (years): "))
    # println(io, show_avg_and_range(model.canopy_evolution.p_DENSEF.itp.itp.coefs, "DENSEF         (°C): "))
    # println(io, show_avg_and_range(model.canopy_evolution.p_HEIGHT.itp.itp.coefs, "HEIGHT          (m): "))
    # println(io, show_avg_and_range(model.canopy_evolution.p_LAI.itp.itp.coefs,    "LAI         (m2/m2): "))
    # println(io, show_avg_and_range(model.canopy_evolution.p_SAI.itp.itp.coefs,    "SAI         (m2/m2): "))
end

function simulate!(s::DiscretizedSPAC)
    sol_SPAC = solve_LWFB90(s.ODEProblem);
    s.ODESolution = sol_SPAC;

    @assert (SciMLBase.successful_retcode(sol_SPAC)) "Problem with simulation: Return code of simulation was '$(sol_SPAC.retcode)'"

    # also save datetimes
    s.ODESolution_datetime = LWFBrook90.RelativeDaysFloat2DateTime.(s.ODESolution.t, s.parametrizedSPAC.reference_date)
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
    model = loadSPAC(input_path, input_prefix;
                      simulate_isotopes = simulate_isotopes);
    ####################

    ####################
    # Prepare simulation by discretizing spatial domain
    simulation = LWFBrook90.setup(model);

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
