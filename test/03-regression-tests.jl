# Regression tests

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



# TODO(bernhard): while below testset works consistently on MacBook Pro,
#                 in the GithubActions the values are slightly different.
#                 TODO: investigate why. (DiffEq.jl, integration testing, set some seed?)
@testset "run_example" begin
  @githash_time example_result = LWFBrook90.run_example()
  # MacBookPro15,2-git-c4275ee: 0.020843 seconds (202.36 k allocations: 16.009 MiB)
  # MacBookPro15,2-git-c4275ee: 0.018723 seconds (202.36 k allocations: 16.010 MiB)

  idx_of_state_variables = 1:(7+(example_result["NLAYER"]-1))
  @test_skip example_result["solution"].u[10][idx_of_state_variables] ≈ [
    0.0
    0.0
    1.1102230246251565e-16
    0.0
    0.0
    0.0
    27.753841212565565
    20.90274150788789
    64.73779777927555
    98.40498418395806
    122.79978640636511
    150.91326550464527
    1.9148395579592201
  ]

  @test_skip example_result["solution"].u[50][idx_of_state_variables] ≈ [
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    25.198288503906493
    20.322289691259886
    58.62009000784807
    119.05056620262427
    151.1460867821672
    150.87187388772858
    1.9139814657581964
  ]

  @test_skip example_result["solution"].u[100][idx_of_state_variables] ≈ [
    0.0
    0.0
    0.0
    0.1283689357747232
    0.0
    0.00641844678873616
    25.457910149010225
    20.362161808370338
    60.64778413324381
    112.08682989771916
    154.74527695278573
    180.6160884889641
    1.9137807295860196
  ]
end

