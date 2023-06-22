module LWFBrook90

using SciMLBase       # instead of loading the full DifferentialEquations
using OrdinaryDiffEq  # instead of loading the full DifferentialEquations
using DiffEqCallbacks # instead of loading the full DifferentialEquations
using ProgressLogging
using RecipesBase, PlotUtils, Measures
using LinearAlgebra
using StatsBase: mean, weights
using ComponentArrays
using UnPack: @unpack
using Dates: now, Day, dayofyear
using Printf: @sprintf
using Interpolations: interpolate, extrapolate, NoInterp, Gridded, Constant, Next, Previous, Flat, Throw, scale, BSpline, linear_interpolation

# DEVELOPERS: EXPORTED ELEMENTS CONSTITUTE THE API AND SHOULD REMAIN DOCUMENTED.
export SPAC, DiscretizedSPAC, loadSPAC, setup, simulate!, remakeSPAC
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
    - `pars.canopy_evolution`: canopy parameters (LAI, SAI, DENSEF, HEIGHT) as relative in percent. Either:
        - `NamedTuple`: (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel = 100, LAI_rel = (DOY\_Bstart, Bduration, DOY\_Cstart, Cduration, LAI\_perc\_BtoC, LAI\_perc\_CtoB)), containing constant values and parameters for LAI\_relative interpolation
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

include("../examples/func_run_example.jl")

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
    remakeSPAC(discrSPAC::DiscretizedSPAC;
                requested_tspan = nothing,
                soil_output_depths_m::Vector = zeros(Float64, 0),
                kwargs...)
or

    remakeSPAC(parametrizedSPAC::SPAC;
                requested_tspan = nothing,
                soil_output_depths_m::Vector = zeros(Float64, 0),
                kwargs...)

Generates a copy of the provided SPAC or DiscretizedSPAC and modifies all the parameter
that are provided as kwargs. This is useful running the same model with a range of different
parameter.

Possible kwargs are:
- `soil_horizons = (ths_ = 0.4, Ksat_mmday = 3854.9, alpha_per_m = 7.11, gravel_volFrac = 0.1)`
- `LAI_rel = (DOY_Bstart = 115,)`
- `root_distribution = (beta = 0.88, z_rootMax_m = -0.6,)`
- `params = (DRAIN=.33, BYPAR=1, IDEPTH_m=0.67, INFEXP=0.33,
            ALB=0.15, ALBSN=0.7, RSSA=720., PSICR=-1.6, FXYLEM=0.4, MXKPL=16.5, MAXLAI=9.999,
            GLMAX=.00801, R5=235., CVPD=1.9, CINTRL=0.18,)`
- `IC_soil = (PSIM_init_kPa = -3.0, delta18O_init_permil = -15.55, )`
- `IC_scalar = ICscalar_tochange`

When the argument `soil_horizons` contains scalar values, the parameter of the toplayer are set
to this value and parameter of all other layers are modified proportionally. Alternatively,
the argument `soil_horizons` can be provided with a vector or vectors with a value for each soil
horizon. A mix of vectors and scalars is also possible.

"""
function remakeSPAC(discrSPAC::DiscretizedSPAC;
                requested_tspan = nothing,
                soil_output_depths_m::Vector = zeros(Float64, 0),
                kwargs...)
    return remakeSPAC(discrSPAC.parametrizedSPAC;
                requested_tspan = requested_tspan,
                soil_output_depths_m = soil_output_depths_m,
                kwargs...)
end
function remakeSPAC(parametrizedSPAC::SPAC;
                requested_tspan = nothing,
                soil_output_depths_m::Vector = zeros(Float64, 0),
                kwargs...) # kwargs collects all others
    # dump(kwargs) # kwargs contains NamedTuples to be modified
    modifiedSPAC = deepcopy(parametrizedSPAC)
    for curr_change in keys(kwargs)
        @assert values(kwargs)[curr_change] isa NamedTuple ||
            # expect a NamedTuple of values to ovewrite, or in the case of IC_scalar a DataFrame
            (curr_change == :IC_scalar && values(kwargs)[curr_change] isa DataFrame) """
        Argument $curr_change is not a NamedTuple. Append single values with comma, e.g. (beta=12,)
        """
        if (curr_change     == :soil_horizons)
            modifiedSPAC = remake_soil_horizons(    modifiedSPAC, values(kwargs)[curr_change])
        elseif (curr_change == :LAI_rel)
            modifiedSPAC = remake_LAI(              modifiedSPAC, values(kwargs)[curr_change])
        elseif (curr_change == :params)
            modifiedSPAC = remake_params(           modifiedSPAC, values(kwargs)[curr_change])
        elseif (curr_change == :root_distribution)
            modifiedSPAC = remake_root_distribution(modifiedSPAC, values(kwargs)[curr_change])
        elseif (curr_change == :IC_soil)
            modifiedSPAC = remake_IC_soil(          modifiedSPAC, values(kwargs)[curr_change])
        elseif (curr_change == :IC_scalar)
            modifiedSPAC = remake_IC_scalar(        modifiedSPAC, values(kwargs)[curr_change])
        else
            error("Unknown argument provided to remake(): $curr_change")
        end
    end
    return setup(modifiedSPAC;
                requested_tspan = requested_tspan,
                soil_output_depths_m = soil_output_depths_m)
end
function remake_soil_horizons(spac, changesNT)
    shp_names = Dict(:ths_           => :p_THSAT,
                     :Ksat_mmday     => :p_KSAT,
                     :alpha_per_m    => :p_MvGα,
                     :gravel_volFrac => :p_STONEF)
    for (key, val) in zip(keys(changesNT), changesNT)
        @assert key ∈ keys(shp_names) "Unclear how to remake '$key' provided to soil_horizons."
        N_horizons = length(spac.pars.soil_horizons.shp)
        if length(val) == 1
            # if only 1 value: modify toplayer to value and all horizons proportionally
            ratio = val/getproperty(spac.pars.soil_horizons.shp[1], shp_names[key])
            target = fill(NaN, N_horizons)
            for idx in 1:N_horizons
                target[idx] = ratio * getproperty(spac.pars.soil_horizons.shp[idx], shp_names[key])
            end
        else
            # if more than 1 value: modify all the layers accordingly
            @assert length(val) == N_horizons """
                Argument '$key' provided to soil_horizons, must either be a single value (for the toplayer) or a vector of length $(N_horizons) (for each soil horizon).
                """
            target = val
        end

        # update values with target values
        for idx in 1:N_horizons
            setproperty!(spac.pars.soil_horizons.shp[idx], shp_names[key], target[idx])
        end
    end
    return spac
end
function remake_LAI(spac, changesNT)
    @assert spac.pars.canopy_evolution isa NamedTuple """
        LAI evolution can only be modified if the SPAC was generated with a parametrized variant of LAI.
        """
    LAI_names = [:DOY_Bstart,:Bduration ,:DOY_Cstart,:Cduration ,:LAI_perc_BtoC,:LAI_perc_CtoB]
    for (key, val) in zip(keys(changesNT), changesNT)
        @assert key ∈ LAI_names "Unclear how to remake '$key' provided to LAI."
    end
    # create new LAI reusing the old one and overwriting whats defined
    new_LAI_pars = (;spac.pars.canopy_evolution.LAI_rel..., changesNT...) # https://stackoverflow.com/a/60883705
    spac.pars = (;spac.pars...,
                  canopy_evolution = (;spac.pars.canopy_evolution..., LAI_rel = new_LAI_pars))
    return spac
end
function remake_root_distribution(spac, changesNT)
    allowed_names = [:beta, :z_rootMax_m]
    for (key, val) in zip(keys(changesNT), changesNT)
        @assert key ∈ allowed_names "Unclear how to remake '$key' provided to params."
    end
    # create new root_distribution reusing the old one and only overwriting whats defined by changesNT
    new_root_distribution = (;spac.pars.root_distribution..., changesNT...) # https://stackoverflow.com/a/60883705
    spac.pars = (;spac.pars..., root_distribution = new_root_distribution)
    return spac
end
function remake_params(spac, changesNT)
    allowed_names = keys(spac.pars.params)
    for (key, val) in zip(keys(changesNT), changesNT)
        @assert key ∈ allowed_names "Unclear how to remake '$key' provided to params."
    end
    # create new params reusing the old one and only overwriting whats defined by changesNT
    new_params = (;spac.pars.params..., changesNT...) # https://stackoverflow.com/a/60883705
    spac.pars = (;spac.pars..., params = new_params)
    return spac
end
function remake_IC_soil(spac, changesNT)
    allowed_names = keys(spac.pars.IC_soil)
    for (key, val) in zip(keys(changesNT), changesNT)
        @assert key ∈ allowed_names "Unclear how to remake '$key' provided to params."
    end
    # create new IC_soil reusing the old one and only overwriting whats defined by changesNT
    new_IC_soil = (;spac.pars.IC_soil..., changesNT...) # https://stackoverflow.com/a/60883705
    spac.pars = (;spac.pars..., IC_soil = new_IC_soil)
    return spac
end
function remake_IC_scalar(spac, newICsoil_DF)
    @assert all(
        names(newICsoil_DF) .==
            ["u_GWAT_init_mm", "u_INTS_init_mm", "u_INTR_init_mm", "u_SNOW_init_mm", "u_CC_init_MJ_per_m2", "u_SNOWLQ_init_mm"]) """
    Unexpected column names in remake_IC_scalar().
    """
    # overwrite the entire DataFrame
    # create new pars with an overwritten IC_scalar
    spac.pars = (;spac.pars..., IC_scalar = newICsoil_DF)
    return spac
end

"""
    setup(parametrizedSPAC::SPAC;
        requested_tspan = nothing,
        soil_output_depths_m::Vector = zeros(Float64, 0))

Takes the definition of SPAC model and discretize to a system of ODEs that can be solved by
the package DifferentialEquations.jl

If needed, the computational grid of the soil is refined to output values at specific depths,
e.g. by doing `setup(SPAC; soil_output_depths_m = [-0.35, -0.42, -0.48, -1.05])`.

The argument `requested_tspan` is a tuple defining the duration of the simulation either
with a start- and an end-day relative to the reference date: `setup(SPAC; requested_tspan = (0,150))`.
or alternatively as a tuple of 2 DataTime's.
Note that the reference date given by the `SPAC.reference_date` refers to the day 0.
Further note that it is possible to provide a tspan not starting at 0, e.g. (120,150).
In that case, initial conditions are applied tspan[1], but atmospheric forcing is correctly taken into account.
"""
function setup(parametrizedSPAC::SPAC;
                requested_tspan = nothing,
                soil_output_depths_m::Vector = zeros(Float64, 0), # e.g. # soil_output_depths_m = [-0.35, -0.42, -0.48, -1.05]
                # Define what is close enough for soil_output_depths_m
                ε = 0.050   # thickness of layer to be inserted, [m]
                )

    @assert all(soil_output_depths_m .< 0)

    if (!isnothing(requested_tspan) && requested_tspan[1] isa DateTime)
        requested_tspan = LWFBrook90.DateTime2RelativeDaysFloat.(
            requested_tspan, parametrizedSPAC.reference_date)
    end
    if (!isnothing(requested_tspan) && !(requested_tspan[1] ≈ 0))
        @warn """
        Requested time span doesn't start at 0. This is supported and correctly takes into account atmospheric forcing.
        Note, however, that initial conditions are applied to t=$(requested_tspan[1]), i.e. at $(parametrizedSPAC.reference_date + Second(floor(requested_tspan[1] * 24*3600))).
        """
    end
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
            modifiedSPAC.pars.params[:QDEPTH_m];
            ε = ε)
    Δz_refined = refined_soil_discretizationDF.Upper_m - refined_soil_discretizationDF.Lower_m
    # if rootden (and initial conditions were given parametrically redo them):
    modifiedSPAC.pars.root_distribution
    if (modifiedSPAC.pars.root_distribution isa DataFrame)
        # keep as is
    elseif (modifiedSPAC.pars.root_distribution isa NamedTuple)
        overwrite_rootden!(refined_soil_discretizationDF, modifiedSPAC.pars.root_distribution, Δz_refined)
        overwrite_IC!(     refined_soil_discretizationDF, modifiedSPAC.pars.IC_soil, modifiedSPAC.solver_options.simulate_isotopes)
    end

    # Discretize the model in space as `soil_discretization`
    final_soil_discretizationDF = map_soil_horizons_to_discretization(modifiedSPAC.pars.soil_horizons, refined_soil_discretizationDF)#computational_grid)

    # Update soil_discretization in underlying SPAC model
    modifiedSPAC.soil_discretization = (
        Δz = Δz_refined,
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
    # Define ODE problem which consists at the least of
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
        @warn "Overwriting tspan defined in SPAC $(modifiedSPAC.tspan) with provided value of $requested_tspan"
        tspan_to_use = requested_tspan

        # Update tspan in underlying SPAC model
        modifiedSPAC.tspan = requested_tspan
    end

    # Seperate updating of different states (INTS, INTR, SNOW, CC, SNOWLQ are updated once per
    # day while GWAT and SWATI are updated continuously) is implemented by means of operator
    # splitting using a callback function for the daily updates and a ODE RHS (right hand
    # side) for the continuous update.

    cb_func = define_LWFB90_cb() # define callback functions
    @assert !any(ismissing.(u0)) """
    There are missing values in the provided initial conditions `u0`. Please correct!"""

    # Note that we require some workarounds for LWFBrook90.jl.
    # Namely:
    # - to accomodate a state vector (array) that contains NaNs, e.g. isotopic concentrations
    #   that are undefined when a compartment is empty (e.g. SNOW = 0mm, δ_SNOW = NA):
    #   a1) modify `unstable_check` (to not use these NaN)
    #   a2) modify `internalnorm` for adaptive time stepping
    # - to set default solver algorithm
    # - to set default time stepping criteria
    # Note that documentation for solve()-arguemnts can be found at:
    # https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options
    ode_LWFBrook90 =
        ODEProblem(LWFBrook90.f_LWFBrook90R, u0, tspan_to_use, p;
                    callback = cb_func,
                    # below the default value for keyword-arguments for solve(),
                    # note they can still be overwritten when calling solve() (or simulate!())
                    # arguments for solve(), can still be overwritten
                    progress_name = "SPAC simulating...",
                    progress = true,
                    saveat = tspan_to_use[1] : 1 : tspan_to_use[2], # in days -> i.e. daily output
                    save_everystep = true, # saves additionally to saveat
                    alg=Tsit5(), reltol = 1e-5,
                    adaptive = true, internalnorm = LWFBrook90.norm_to_use, # fix adaptivity norm for NAs
                    unstable_check = LWFBrook90.unstable_check_function,         # fix instability norm for NAs
                    dt    = 1e-3,                            # dt is initial dt, but can be changed adaptively
                    dtmax = 60/60/24, # 60min max time step
                    )
    # former kwargs included:
        # reltol = 1e-5, # abstol = 1e-6, # default: abstol = 1e-6, reltol = 1e-3
        # tstops = sol_working.t, #[0.223, 0.4],
        # tstops = tspan[1]:0.00001:tspan[2], adaptive = false
        # ImplicitEuler(autodiff=false);  # for stiff problems (~6.7s for 2.5days of Hammel_loam-NLayer-103)
        # Rodas4P(autodiff=false);        # for stiff problems (~3.2s for 2.5days of Hammel_loam-NLayer-103)
        # TRBDF2(autodiff=false);
        # Rosenbrock23(autodiff=false);   # for stiff problems   (~2.3s for 2.5days of Hammel_loam-NLayer-103)
        # Tsit5(); abstol = 1e-6, reltol = 1e-4, # Tsit5 recommended for non-stiff problems (~ 12.4s for 2.5days of Hammel_loam-NLayer-103)
        # Tsit5(); abstol = 1e-6, reltol = 1e-5, # Tsit5 recommended for non-stiff problems (~ 20.9s for 2.5days of Hammel_loam-NLayer-103)
        # Tsit5(); # Tsit5 recommended for non-stiff problems (~ 6.4s for 2.5days of Hammel_loam-NLayer-103)
        # AutoTsit5(Rosenbrock23(autodiff=false)); reltol = 1e-5, # recommended for problems of unknown stiffnes: (~5.0s)
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


"""
    simulate!(s::DiscretizedSPAC; kwargs...)

Simulates a SPAC model and stores the solution in s.ODESolution.

`kwargs...` are passed through to solve(SciML::ODEProblem; ...) and are
documented under https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options
"""
function simulate!(s::DiscretizedSPAC; assert_retcode = true, kwargs...)
    if (:saveat ∈ keys(kwargs))
        @info """
          Start of simulation at $(now()).
                Saving intermediate results (`saveat=`) between: $(extrema(values(kwargs)[:saveat])) days
        """
    elseif (:saveat ∈ keys(s.ODEProblem.kwargs))
        @info """
          Start of simulation at $(now()).
                Saving intermediate results (`saveat=`) between: $(extrema(values(s.ODEProblem.kwargs)[:saveat])) days
        """
    else
         @info "  Start of simulation at $(now())."
    end

    @time s.ODESolution = solve(s.ODEProblem; kwargs...)

    @info "  Time steps for solving: $(s.ODESolution.destats.naccept) ($(s.ODESolution.destats.naccept) accepted out of $(s.ODESolution.destats.nreject + s.ODESolution.destats.naccept) total)"
    @info "  End of simulation at $(now())."

    if assert_retcode
        @assert (SciMLBase.successful_retcode(s.ODESolution)) "Problem with simulation: Return code of simulation was '$(s.ODESolution.retcode)'"
    end

    # also store datetimes
    s.ODESolution_datetime = LWFBrook90.RelativeDaysFloat2DateTime.(s.ODESolution.t, s.parametrizedSPAC.reference_date)
    return nothing
end


############################################################################################
############################################################################################
############################################################################################

@doc raw"""
    run_simulation(args) (deprecated!)

Runs a simulation defined by input files within a folder and returns solved DiscretizedSPAC.

The function `run_simulation()` takes as single argument a vector of two or three strings defining
the `input_path` and `input_prefix` of a series of input definition files (and `true`/`false` whether
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

    @info simulation.ODESolution.destats
    @info "Time steps for solving: $(simulation.ODESolution.destats.naccept) ($(simulation.ODESolution.destats.naccept) accepted out of $(simulation.ODESolution.destats.nreject + simulation.ODESolution.destats.naccept) total)"
    @show now()

    return (simulation, input_prefix, input_path)
end

end # module
