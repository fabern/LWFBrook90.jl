using LWFBrook90
using SciMLBase
using DataFrames
using Test: @testset, @test, @test_throws, @test_broken, @test_skip, @test_logs
using CSV: File
using Random
using Printf
using Logging
using Dates: today

# Note the git hash-string and repository status (can be optionally used in filenames for plots etc.)
git_status_string = "$(today())-git+"*
            "__"*chomp(Base.read(`git branch --show-current`, String))*
            "__"     *chomp(Base.read(`git rev-parse --short HEAD`, String))*
            ifelse(length(Base.read(`git status --porcelain`, String))==0, "+gitclean","+gitdirty")*
            "__"

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

# A flag that determines if tests are run on a CI system
is_a_CI_system = issubset(["GITHUB_ACTION"], collect(keys(ENV))) # checks if ENV["GITHUB_ACTION"] exists
@show is_a_CI_system


Random.seed!(1234)

# cd("test")

# if !is_a_CI_system; include("00-plot-if-not-CI-system.jl"); end

include("01-unit-tests.jl")

plot_flag = false; high_resolution_flag = false;
# plot_flag = true; high_resolution_flag = false; using Plots, Measures
# plot_flag = true; high_resolution_flag = true; using Plots, Measures # run this line for plotting outside of default testing

include("02-integration-tests.jl")

include("03-regression-tests.jl")
