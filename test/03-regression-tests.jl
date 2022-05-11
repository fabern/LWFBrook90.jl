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
@testset "regression-nodelta-run_example" begin

  # run simulation
  @githash_time example_result = LWFBrook90.run_example()
  # amberMBP-git-c4275ee: 0.020843 seconds (202.36 k allocations: 16.009 MiB)
  # amberMBP-git-c4275ee: 0.018723 seconds (202.36 k allocations: 16.010 MiB)
  # amberMBP-git-eae940b: 0.028672 seconds (93.92 k allocations: 12.045 MiB)
  # amberMBP-git-eae940b: @btime: 14.003 ms (93919 allocations: 12.05 MiB)

  # extract required data from solution object
  idx_u_scalar_amounts = example_result["solution"].prob.p[1][4][6]
  idx_u_vector_amounts = example_result["solution"].prob.p[1][4][4]
  idx_u = [idx_u_scalar_amounts; idx_u_vector_amounts]

  u_ref = reshape(example_result["solution"][idx_u,:,1:1:end], length(idx_u), :)

  # test or overwrite
  fname = "../examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_u_sol_reference.jld2"
  if task == "test"
    loaded_u_ref = load(fname, "u_ref")
    @test all(abs.((u_ref  .- loaded_u_ref) ./ (loaded_u_ref  .+ eps(Float64)) .< 1e-3)) # adding eps for values where _ref is zero
  elseif task == "overwrite"
    jldsave(fname; u_ref)
  else
    # do nothing
  end

# @test_skip example_result["solution"].u[10][idx_u] ≈ [
#     0.0
#     0.0
#     1.1102230246251565e-16
#     0.0
#     0.0
#     0.0
#     27.753841212565565
#     20.90274150788789
#     64.73779777927555
#     98.40498418395806
#     122.79978640636511
#     150.91326550464527
#     1.9148395579592201
#   ]

#   @test_skip example_result["solution"].u[50][idx_u] ≈ [
#     0.0
#     0.0
#     0.0
#     0.0
#     0.0
#     0.0
#     25.198288503906493
#     20.322289691259886
#     58.62009000784807
#     119.05056620262427
#     151.1460867821672
#     150.87187388772858
#     1.9139814657581964
#   ]

#   @test_skip example_result["solution"].u[100][idx_u] ≈ [
#     0.0
#     0.0
#     0.0
#     0.1283689357747232
#     0.0
#     0.00641844678873616
#     25.457910149010225
#     20.362161808370338
#     60.64778413324381
#     112.08682989771916
#     154.74527695278573
#     180.6160884889641
#     1.9137807295860196
#   ]
end



@testset "regression-delta-isoBEAdense2010-18" begin
  # run simulation
  @githash_time sol = run_simulation(
    ["../examples/isoBEAdense2010-18-reset-FALSE-input/" "isoBEAdense2010-18-reset-FALSE" "true"]);
  # amberMBP-git-9d3342b: 59.289101 seconds (333.36 M allocations: 53.804 GiB, 20.36% gc time)
  # sol[1] # solution
  # sol[2] # input_prefix
  # sol[3] # input_path

  # extract required data from solution object
  # (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
  #         get_auxiliary_variables(sol[1])
  # u_SWATI
  # u_aux_θ
  u_δ = get_δ(sol[1])
  SWAT_d18O_ref  = u_δ.SWAT.d18O[:, 1:30:end]
  SWAT_d2H_ref   = u_δ.SWAT.d2H[ :,  1:30:end]
  PREC_d18O_ref  = u_δ.PREC.d18O[1:30:end]
  PREC_d2H_ref   = u_δ.PREC.d2H[ 1:30:end]
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
    loaded_SWAT_d18O_ref = load(fname*"_OUTPUT-SWAT_d18O_reference.jld2", "SWAT_d18O_ref")
    loaded_SWAT_d2H_ref  = load(fname*"_OUTPUT-SWAT_d2H_reference.jld2",  "SWAT_d2H_ref")
    @test all(abs.((SWAT_d18O_ref .- loaded_SWAT_d18O_ref) ./ (loaded_SWAT_d18O_ref .+ eps(Float64))) .< 1e-5) # adding eps for values where _ref is zero
    @test all(abs.((SWAT_d2H_ref  .- loaded_SWAT_d2H_ref ) ./ (loaded_SWAT_d2H_ref  .+ eps(Float64))) .< 1e-5) # adding eps for values where _ref is zero
  elseif task == "overwrite"

    # overwrite output
    jldsave(fname*"_OUTPUT-SWAT_d18O_reference.jld2"; SWAT_d18O_ref)
    jldsave(fname*"_OUTPUT-SWAT_d2H_reference.jld2"; SWAT_d2H_ref)

    # # plot
    # using Plots, Measures
    # optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
    # # pl1 = LWFBrook90.ISO.plotisotopes(sol_LWFBrook90);
    # pl2 = LWFBrook90.ISO.plotisotopes(
    #     sol[1], optim_ticks;
    #     layout = grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4]),
    #     size=(1000,1400), dpi = 300, leftmargin = 15mm);
    # plot!(pl2, link = :x)
    # savefig("../examples/isoBEAdense2010-18-reset-FALSE-input/idoBEAdense2010-18-reset-FALSE_OUTPUT-simulated.png")
  else
    # do nothing
  end

end
