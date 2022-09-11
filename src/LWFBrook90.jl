module LWFBrook90

using OrdinaryDiffEq  # instead of loading the full DifferentialEquations
using DiffEqCallbacks # instead of loading the full DifferentialEquations
using RecipesBase
using LinearAlgebra
using StatsBase: mean, weights

using Dates: now
# using Infiltrator

export read_inputData
export discretize_soil, Rootden_beta_
export define_LWFB90_p, define_LWFB90_u0, solve_LWFB90
export KPT_SOILPAR_Mvg1d, KPT_SOILPAR_Ch1d
export RelativeDaysFloat2DateTime, plot_LWFBrook90

export run_simulation, plot_and_save_results, find_indices
export get_auxiliary_variables, get_θ, get_δ, get_δsoil

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

include("../examples/BEA2016-reset-FALSE-input/func_run_example.jl") # defines RelativeDaysFloat2DateTime

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
    # args = ["test/test-assets/Hammel-2001/input-files-ISO" "Hammel_loam-NLayer-27-RESET=FALSE" "true"]
    @show now()
    @show args

    @assert length(args) >= 2

    input_path = args[1]
    input_prefix = args[2]
    if (length(args) == 3)
        simulate_isotopes = args[3] == "true"
    else
        simulate_isotopes = false
    end
    @show simulate_isotopes

    (input_meteoveg,
    _input_meteoiso, # NOTE: _ indicates that it is possibly unused (depends on simulate_isotopes)
    input_meteoveg_reference_date,
    input_param,
    input_storm_durations,
    input_initial_conditions,
    input_soil_horizons,
    simOption_FLAG_MualVanGen) = read_inputData(input_path, input_prefix;
                                                simulate_isotopes = simulate_isotopes);

    input_soil_discretization = discretize_soil(input_path, input_prefix);

    Reset = false                          # currently only Reset = 0 implemented
    compute_intermediate_quantities = true # Flag whether ODE containes additional quantities than only states

    (ψM_initial, _δ18O_initial, _δ2H_initial), p = define_LWFB90_p(
        input_meteoveg,
        _input_meteoiso, # NOTE: _ indicates that it is possibly unused (depends on simulate_isotopes)
        input_meteoveg_reference_date,
        input_param,
        input_storm_durations,
        input_soil_horizons,
        input_soil_discretization,
        simOption_FLAG_MualVanGen;
        Reset = Reset,
        # soil_output_depths = collect(-0.05:-0.05:-1.1),
        # soil_output_depths = [-0.1, -0.5, -1.0, -1.5, -1.9],
        compute_intermediate_quantities = compute_intermediate_quantities,
        simulate_isotopes = simulate_isotopes);
    u0, p = define_LWFB90_u0(p, input_initial_conditions,
        ψM_initial, _δ18O_initial, _δ2H_initial, # NOTE: _ indicates that it is possibly unused (depends on simulate_isotopes)
        compute_intermediate_quantities;
        simulate_isotopes = simulate_isotopes);

    tspan = (minimum(input_meteoveg[:, "days"]), maximum(input_meteoveg[:, "days"])) # simulate all available days
    # using Plots
    # using Interpolations: interpolate, extrapolate, NoInterp, Gridded, Constant, Next, Previous, Flat, Throw, scale
    # scatter(input_meteoveg.days[1:10], input_meteoveg.PRECIN[1:10])
    # t_to_eval = 0:0.2:10
    # plot!(t_to_eval, p[2].p_PREC(t_to_eval), xtick = 1:10)
    # plot!(t_to_eval,extrapolate(interpolate((input_meteoveg.days .- 0.00001, ), input_meteoveg.PRECIN, Gridded(Constant{Previous}())), Throw())(t_to_eval))
    sol_LWFBrook90 = solve_LWFB90(u0, tspan, p)
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
