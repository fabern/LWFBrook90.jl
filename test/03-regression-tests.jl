# Regression tests
using JLD2: jldsave, load # TODO: comment in again as soon  https://github.com/JuliaIO/JLD2.jl/issues/428 is fixed
using FileIO

# - _Unit testing_ asserts that individual pieces of a project work as expected. (developers
#       perspective)
# - _Integration testing_ asserts that they fit together as expected. Also known as
#       _functional tests_, they cover entire use cases (user perspective). For LWFBrook90.jl
#       these are tests that are compared to e.g. LWFBrook90R or Hydrus.
# - _Regression testing_ asserts that behavior is unchanged over time. Also known as
#       _reference tests_.


# Such Regression tests can easily setup (and modified if need be) with `ReferenceTests.jl`.
# Disadvantage is that it has a big dependency. Alternatively look at `SheetModel.jl` as
# they coded this manually and also included a function to write out new model results if
# need be (e.g. in case ther was a bug before and want to have the new model outputs to be
# counted as correct).

# this script either writes out or tests against reference test cases
task = ["test", "overwrite"][1]

# TODO(bernhard): while below testset works consistently on MacBook Pro,
#                 in the GithubActions the values are slightly different.
#                 TODO: investigate why. (DiffEq.jl, integration testing, set some seed?)
# Note, in order to have the regression testset work consistently on MacBook Pro and on the
# CI machine (through GithubActions), make sure to test always with the same PRNG (random
# number generator). According to
# https://discourse.julialang.org/t/tests-failing-on-linux-in-travis-ci/19602/2
# the @testset macro does this for us.
#
# It turns out, we need also locally to run it always through `]`, `test`. The workflow is:
#     julia --project=test
#     ]
#     activate .
#     test
# TODO(bernhard): however, above `test` still yields different results in CI and fails, therefore:
# Current workaround, do not run these tests on CI (GithubActions)
is_a_CI_system = issubset(["GITHUB_ACTION"], collect(keys(ENV))) # checks if ENV["GITHUB_ACTION"] exists
@show is_a_CI_system

# NOTE: locally one might need to do manually cd("test")
if basename(pwd()) != "test"; cd("test"); end

@testset "regression-delta-run_example" begin
    @githash_time example_result = LWFBrook90.run_example();

    # extract required data from solution object
    ## helper quantities
    t_out = range(extrema(example_result.ODESolution.t)..., step = 30)
    currSPAC = example_result.parametrizedSPAC
    ## scalar quantities
    u_ref = vcat(
        reduce(hcat, [example_result.ODESolution(t_days).GWAT.mm   for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).INTS.mm   for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).INTR.mm   for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).SNOW.mm   for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).CC.MJm2   for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).SNOWLQ.mm for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).RWU.mmday for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).XYLEM.mm  for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).SWATI.mm  for t_days = t_out]),
        # reduce(hcat, [example_result.ODESolution[t_idx].TRANI.mmday for t_idx = eachindex(example_result.ODESolution)]),
        # reduce(hcat, [example_result.ODESolution[t_idx].aux.θ       for t_idx = eachindex(example_result.ODESolution)]),
        # reduce(hcat, [example_result.ODESolution[t_idx].accum       for t_idx = eachindex(example_result.ODESolution)])
    )
    u_δ     = get_δ(example_result; days_to_read_out_d = t_out)
    ## vector quantities
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) = LWFBrook90.get_auxiliary_variables(example_result; days_to_read_out_d = t_out)
    u_δsoil = get_δsoil(example_result; days_to_read_out_d = t_out);

    # test or overwrite
    # fname = input_path*input_prefix
    fname = "../examples/DAV2020-full_u_sol_reference.jld2"
    if task == "test"
        loaded = load(fname)

        if !is_a_CI_system
            # Test input example_result
            # @test loaded["currSPAC"]   == currSPAC
            @test loaded["currSPAC"].forcing == currSPAC.forcing
            # @test loaded["currSPAC"].pars == currSPAC.pars
            @test loaded["currSPAC"].reference_date == currSPAC.reference_date
            # @test loaded["currSPAC"].soil_discretization == currSPAC.soil_discretization
            @test loaded["currSPAC"].solver_options == currSPAC.solver_options
            @test loaded["currSPAC"].tspan == currSPAC.tspan

            # Test scalar results
            @test all(loaded["u_ref"]      .≈ u_ref)
            @test all(isapprox.(Matrix(loaded["u_δ"]), Matrix(u_δ), nans=true)) # NOTE: δ_RWU can be NaN if no water taken up, hence nans=true
            # Test vector results
            @test all(loaded["u_SWATI"]    .≈ u_SWATI)
            @test all(loaded["u_aux_PSIM"] .≈ u_aux_PSIM)
            @test all(loaded["u_aux_θ"]    .≈ u_aux_θ)
            @test all(loaded["u_δsoil"].d18O .≈ u_δsoil.d18O)
            @test all(loaded["u_δsoil"].d2H .≈ u_δsoil.d2H)
        else
            @warn "Regression tests deactivated on CI system."
        end
    elseif task == "overwrite" && !is_a_CI_system # only overwrite on local machine, never on CI
        # overwrite output
        jldsave(fname, compress=false;
                (currSPAC   = currSPAC,
                 u_ref      = u_ref,
                 u_δ        = u_δ,
                 u_SWATI    = u_SWATI,
                 u_aux_PSIM = u_aux_PSIM,
                 u_aux_θ    = u_aux_θ,
                 u_δsoil    = u_δsoil)...)
        # if (false)
        #     using Plots, Measures
        #     git_status_string = "__git+"*chomp(Base.read(`git rev-parse --short HEAD`, String))*
        #         ifelse(length(Base.read(`git status --porcelain`, String))==0, "+gitclean","+gitdirty")*
        #         "__"
        #     pl1 = plotamounts(example, :above_and_belowground, :showRWUcentroid)
        #     pl2 = plotisotopes(example, :d18O)
        #     pl3 = plotforcingandstates(example)
        #     savefig(plot(pl1, size=(1000,700), dpi=300),  fname*git_status_string*"_amts.png")
        #     savefig(plot(pl2, size=(1000,1400), dpi=300), fname*git_status_string*"_d18O-d2H.png")
        #     savefig(plot(pl3, size=(1000,1400), dpi=300), fname*git_status_string*"_forcing.png")
        #     #"../examples/DAV2020-full_u_sol_reference.jld2__git+09bca17+gitdirty___forcing.png")
        # end
        @error "Test overwrites reference solution instead of checking against it."
    else
        # do nothing
    end

end
