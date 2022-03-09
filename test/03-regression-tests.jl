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

@testset "run_example" begin
  @githash_time example_result = LWFBrook90.run_example()
  # amberMBP-git-c4275ee: 0.020843 seconds (202.36 k allocations: 16.009 MiB)
  # amberMBP-git-c4275ee: 0.018723 seconds (202.36 k allocations: 16.010 MiB)
  # amberMBP-git-57c3a2b: 0.514906 seconds (1.18 M allocations: 78.538 MiB, 95.01% compilation time)                 # NOTE: this is done from `]`, test
  # amberMBP-git-57c3a2b: 0.611677 seconds (1.18 M allocations: 78.538 MiB, 10.89% gc time, 84.83% compilation time) # NOTE: this is done from `]`, test

  idx_of_state_variables = 1:(7+(example_result["NLAYER"]-1))
  state_var_1 = example_result["solution"].u[10][idx_of_state_variables]
  state_var_2 = example_result["solution"].u[50][idx_of_state_variables]
  state_var_3 = example_result["solution"].u[100][idx_of_state_variables]

  @test state_var_1 ≈ [
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    27.741639126855034
    21.082040356164992
    65.06233753727214
    97.96628781837019
    122.79977817305011
    150.91326522309896
    1.9148395535036713
  ]

  @test state_var_2 ≈ [
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    25.148014734138382
    20.105198659029693
    58.82834300856744
    120.7426501609414
    149.41134071940718
    150.86433794429922
    1.9139813399391166
  ]

  @test state_var_3 ≈ [
    0.0
    0.0
    0.0
    0.1283689357747232
    0.0
    0.00641844678873616
    25.466726726549723
    20.366251907144935
    60.6927711158876
    112.9299734564963
    154.94552290089902
    180.65652824956382
    1.9137705492299724
  ]
end

