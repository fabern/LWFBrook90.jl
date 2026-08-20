# Backwards-compatible entry point; the reusable suite lives in benchmark/.
include(joinpath(@__DIR__, "..", "benchmark", "runbenchmarks.jl"))
