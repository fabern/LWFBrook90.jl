using BenchmarkTools
using LWFBrook90
using SciMLBase: solve, successful_retcode

const SUITE = BenchmarkGroup()
const EXAMPLE_PATH = normpath(joinpath(@__DIR__, "..", "examples", "DAV2020-full"))
const EXAMPLE_PREFIX = "DAV2020-full"
const BENCHMARK_TSPAN = (0.0, 30.0)

function load_benchmark_model(; simulate_isotopes)
    model = loadSPAC(EXAMPLE_PATH, EXAMPLE_PREFIX; simulate_isotopes)
    model.tspan = BENCHMARK_TSPAN
    return model
end

function solve_benchmark_model(simulation)
    solution = solve(simulation.ODEProblem; progress = false)
    @assert successful_retcode(solution)
    return solution
end

const WATER_MODEL = load_benchmark_model(; simulate_isotopes = false)
const ISOTOPE_MODEL = load_benchmark_model(; simulate_isotopes = true)

simulation_suite = SUITE["model simulation"] = BenchmarkGroup()
simulation_suite["DAV2020-full, 30 days, water"] =
    @benchmarkable solve_benchmark_model(simulation) setup = (
        simulation = setup($WATER_MODEL)
    ) evals = 1 samples = 10 seconds = 60
simulation_suite["DAV2020-full, 30 days, water and isotopes"] =
    @benchmarkable solve_benchmark_model(simulation) setup = (
        simulation = setup($ISOTOPE_MODEL)
    ) evals = 1 samples = 10 seconds = 60
