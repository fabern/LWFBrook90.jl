using BenchmarkTools

include(joinpath(@__DIR__, "benchmarks.jl"))

results = run(SUITE; verbose = true)
display(results)

if haskey(ENV, "BENCHMARK_OUTPUT")
    BenchmarkTools.save(ENV["BENCHMARK_OUTPUT"], results)
end
