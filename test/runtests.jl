using LWFBrook90
using DataFrames
using Test: @testset, @test, @test_throws, @test_broken, @test_skip, @test_logs
using CSV: File
using Random
using Printf
using Logging

# A macro for timing that also prints out the git commit hash:
macro githash_time(variable)
    quote
        #model = replace(chomp(Base.read(`sysctl hw.model`, String)), "hw.model: " => "")
        model = "amberMBP"
        hash = chomp(Base.read(`git rev-parse --short HEAD`, String))
        print(model*"-git-"*hash*":")
        @time $(esc(variable))
    end
end
# Alternatively: think about measuring performance along the lines discussed in
# https://discourse.julialang.org/t/benchmarking-tests-to-ensure-prs-dont-introduce-regressions/8630/6
# or then https://github.com/maxbennedich/julia-regression-analysis or ...)

Random.seed!(1234)

# cd("test")
include("01-unit-tests.jl")

include("02-integration-tests.jl")

# include("03-regression-tests.jl")

