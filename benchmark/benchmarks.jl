using BenchmarkTools
using LWFBrook90
using SciMLBase: solve, successful_retcode

const SUITE = BenchmarkGroup()

const EXAMPLE_PATH = normpath(joinpath(@__DIR__, "..", "examples", "DAV2020-full"))
const EXAMPLE_PREFIX = "DAV2020-full"

function load_benchmark_model(; simulate_isotopes, tmax = 30.0)
    model = loadSPAC(EXAMPLE_PATH, EXAMPLE_PREFIX; simulate_isotopes)
    model.tspan = (0.0, tmax)
    return model
end

function solve_benchmark_model(simulation)
    solution = solve(simulation.ODEProblem; progress = false)
    @assert successful_retcode(solution)
    return solution
end

const WATER_MODEL   = load_benchmark_model(; simulate_isotopes = false)
const ISOTOPE_MODEL = load_benchmark_model(; simulate_isotopes = true)

const WATER_MODEL_300   = load_benchmark_model(; simulate_isotopes = false, tmax = 300.0)
const ISOTOPE_MODEL_300 = load_benchmark_model(; simulate_isotopes = true,  tmax = 300.0)

# MODEL SETUP
SUITE["model setup"]["DAV2020-full, 30 days, water"] =
    @benchmarkable setup($WATER_MODEL) evals = 1 samples = 5 seconds = 60
SUITE["model setup"]["DAV2020-full, 30 days, water and isotopes"] =
    @benchmarkable setup($ISOTOPE_MODEL) evals = 1 samples = 5 seconds = 60

# MODEL SIMULATION
SUITE["model simulation"]["DAV2020-full, 30 days, water"] =
    @benchmarkable solve_benchmark_model(simulation30) setup = (simulation30 = setup($WATER_MODEL)) evals = 1 samples = 5 seconds = 60
SUITE["model simulation"]["DAV2020-full, 30 days, water and isotopes"] =
    @benchmarkable solve_benchmark_model(simulation30) setup = (simulation30 = setup($ISOTOPE_MODEL)) evals = 1 samples = 5 seconds = 60

SUITE["model simulation"]["DAV2020-full, 300 days, water"] =
    @benchmarkable solve_benchmark_model(simulation300) setup = (simulation300 = setup($WATER_MODEL_300)) evals = 1 samples = 5 seconds = 60
SUITE["model simulation"]["DAV2020-full, 300 days, water and isotopes"] =
    @benchmarkable solve_benchmark_model(simulation300) setup = (simulation300 = setup($ISOTOPE_MODEL_300)) evals = 1 samples = 5 seconds = 60
