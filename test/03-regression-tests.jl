# Regression tests
using JLD2

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

if !is_a_CI_system
    @testset "regression-nodelta-run_example" begin

        # run simulation
        @githash_time example_result = LWFBrook90.run_example();
        # amberMBP-git-c4275ee: 0.020843 seconds (202.36 k allocations: 16.009 MiB)
        # amberMBP-git-c4275ee: 0.018723 seconds (202.36 k allocations: 16.010 MiB)
        # amberMBP-git-eae940b: 0.028672 seconds (93.92 k allocations: 12.045 MiB)
        # amberMBP-git-eae940b: @btime: 14.003 ms (93919 allocations: 12.05 MiB)
        # amberMBP-git+61a19ed: 0.039381 seconds (117.87 k allocations: 11.615 MiB) 343 time steps

        # extract required data from solution object
        idx_u_scalar_amounts = example_result["solution"].prob.p[1][4].row_idx_scalars;
        idx_u_vector_amounts = example_result["solution"].prob.p[1][4].row_idx_SWATI;
        idx_u = [idx_u_scalar_amounts...; idx_u_vector_amounts];

        u_ref = example_result["solution"][idx_u,1,:];
        # u_ref_d18O= example_result["solution"][idx_u,2,:]
        # u_ref_d2H = example_result["solution"][idx_u,3,:]

        # test or overwrite
        fname = "../examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_u_sol_reference.jld2"
        if task == "test"
            loaded_u_ref = load(fname, "u_ref");
            @test all(abs.((u_ref  .- loaded_u_ref) ./ (loaded_u_ref  .+ eps(Float64))) .< 1e-3) # adding eps for values where _ref is zero
        elseif task == "overwrite"
            jldsave(fname; u_ref);
        else
            # do nothing
        end

        # using Plots, Measures
        # optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
        # pl_final = LWFBrook90.plotlwfbrook90(example_result["solution"], optim_ticks)
        # git_status_string = "__git+"*chomp(Base.read(`git rev-parse --short HEAD`, String))*
        #   ifelse(length(read(`git status --porcelain`, String))==0, "+gitclean","+gitdirty")*
        #   "__"
        # savefig(pl_final, replace(fname, ".jld2"=>git_status_string*".png"))

    end



    @testset "regression-delta-isoBEAdense2010-18" begin
        # run simulation
        @githash_time sol = run_simulation(
        ["../examples/isoBEAdense2010-18-reset-FALSE-input/" "isoBEAdense2010-18-reset-FALSE" "true"]);
        # amberMBP-git-9d3342b: 59.289101 seconds (333.36 M allocations: 53.804 GiB, 20.36% gc time)
        # amberMBP-git-61a19ed: 89.440982 seconds (244.98 M allocations: 37.225 GiB, 9.13% gc time, 49.97% compilation time)
        # amberMBP-git-da150a8: 44.940979 seconds (176.47 M allocations: 33.808 GiB, 15.84% gc time)
        # amberMBP-git-a1872dd: 50.785075 seconds (83.43 M allocations: 4.962 GiB, 2.99% gc time, 96.61% compilation time) 10244 time steps
        # amberMBP-git-a1872dd: 1.800597 seconds (5.00 M allocations: 875.086 MiB, 18.01% gc time) 10244 time steps

        # sol[1] # solution
        # sol[2] # input_prefix
        # sol[3] # input_path

        # extract required data from solution object
        # (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        #         get_auxiliary_variables(sol[1])
        # u_SWATI
        # u_aux_θ
        u_δ = get_δ(sol[1]);
        SWAT_d18O_ref  = u_δ.SWAT.d18O[:, 1:30:end];
        SWAT_d2H_ref   = u_δ.SWAT.d2H[ :,  1:30:end];
        PREC_d18O_ref  = u_δ.PREC.d18O[1:30:end];
        PREC_d2H_ref   = u_δ.PREC.d2H[ 1:30:end];
        # u_δ.INTR.d18O[1:30:end]
        # u_δ.INTR.d2H[ 1:30:end]
        # u_δ.INTS.d18O[1:30:end]
        # u_δ.INTS.d2H[ 1:30:end]
        # u_δ.SNOW.d18O[1:30:end]
        # u_δ.SNOW.d2H[ 1:30:end]
        # u_δ.GWAT.d18O[1:30:end]
        # u_δ.GWAT.d2H[ 1:30:end]

        # test or overwrite
        fname = sol[3]*sol[2]
        if task == "test"
            loaded_SWAT_d18O_ref = load(fname*"_OUTPUT-SWAT_d18O_reference.jld2", "SWAT_d18O_ref");
            loaded_SWAT_d2H_ref  = load(fname*"_OUTPUT-SWAT_d2H_reference.jld2",  "SWAT_d2H_ref");
            @test all(abs.((SWAT_d18O_ref .- loaded_SWAT_d18O_ref) ./ (loaded_SWAT_d18O_ref .+ eps(Float64))) .< 1e-5) # adding eps for values where _ref is zero
            @test all(abs.((SWAT_d2H_ref  .- loaded_SWAT_d2H_ref ) ./ (loaded_SWAT_d2H_ref  .+ eps(Float64))) .< 1e-5) # adding eps for values where _ref is zero
        elseif task == "overwrite"

            # overwrite output
            jldsave(fname*"_OUTPUT-SWAT_d18O_reference.jld2"; SWAT_d18O_ref);
            jldsave(fname*"_OUTPUT-SWAT_d2H_reference.jld2"; SWAT_d2H_ref);

            # plot
            if (false)
                using Plots, Measures
                optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
                pl1 = LWFBrook90.plotlwfbrook90(
                    sol[1], optim_ticks;
                    layout = grid(3,1, heights=[0.2, 0.4, 0.4]),
                    size=(1000,700), dpi=300, leftmargin = 15mm);
                pl2 = LWFBrook90.ISO.plotisotopes(
                    sol[1], optim_ticks;
                    layout = grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4]),
                    size=(1000,1400), dpi = 300, leftmargin = 15mm);
                git_status_string = "__git+"*chomp(Base.read(`git rev-parse --short HEAD`, String))*
                    ifelse(length(read(`git status --porcelain`, String))==0, "+gitclean","+gitdirty")*
                    "__"
                savefig(plot(pl1, size=(1000,700), dpi=300, link=:x),
                        fname*git_status_string*"_theta-psi.png")
                savefig(plot(pl2, size=(1000,1400), dpi=300, link=:x),
                        fname*git_status_string*"_d18O-d2H.png") #"../examples/isoBEAdense2010-18-reset-FALSE-input/idoBEAdense2010-18-reset-FALSE_OUTPUT-simulated.png")
            end
        else
            # do nothing
        end

    end

end