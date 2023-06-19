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
#  or simply julia --project=. test/runtests.jl

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
        # if (true) # Do these manually outside of automatic testing in order not to require Plots pkg
        #     using Plots, Measures
        #     pl1 = plotamounts(example_result, :above_and_belowground, :showRWUcentroid)
        #     pl2 = plotisotopes(example_result, :d18O, :showRWUcentroid)
        #     pl3 = plotforcingandstates(example_result)
        #     savefig(plot(pl1, size=(1000,700), dpi=300),  fname*git_status_string*"_amts.png")
        #     savefig(plot(pl2, size=(1000,1400), dpi=300), fname*git_status_string*"_d18O-d2H.png")
        #     savefig(plot(pl3, size=(1000,1400), dpi=300), fname*git_status_string*"_forcing.png")
        #     #"../examples/DAV2020-full_u_sol_reference.jld2__git+09bca17+gitdirty___forcing.png")
        # end
        error("Test overwrites reference solution instead of checking against it.")
    else
        # do nothing
    end


    if !is_a_CI_system && plot_flag
        fname_illustrations = "out/$git_status_string/TESTSET_DAV2020-regressionTest"
        mkpath(dirname(fname_illustrations))

        pl1 = plotamounts(example_result, :above_and_belowground, :showRWUcentroid)
        pl2 = plotisotopes(example_result, :d18O, :showRWUcentroid)
        pl3 = plotforcingandstates(example_result)
        savefig(plot(pl1, size=(1000,1400), dpi=300),  fname_illustrations*"_amts.png")
        savefig(plot(pl2, size=(1000,700), dpi=300), fname_illustrations*"_d18O-d2H.png")
        savefig(plot(pl3, size=(1000,1400), dpi=300), fname_illustrations*"_forcing.png")
    end

end


@testset "DAV2020modified-FLUXES" begin
    
    # Another simulation testing INFL, BYFL, SRFL, DSFLI, VRFLN among other fluxes
    @githash_time example_result2 = LWFBrook90.run_example(
        simulate_isotopes = true, 
        canopy_evolution = (DENSEF_rel = 100,
                            HEIGHT_rel = 100,
                            SAI_rel    = 100,
                            LAI_rel = (DOY_Bstart = 120,
                                        Bduration  = 21,
                                        DOY_Cstart = 270,
                                        Cduration  = 60,
                                        LAI_perc_BtoC = 100,
                                        LAI_perc_CtoB = 20)),
        Δz_thickness_m = fill(0.02, 40), 
        root_distribution = (beta = 0.98, z_rootMax_m=-0.5),
        IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -9.0, delta2H_init_permil = -11.0), 
        IC_scalar = (
            amount = (u_GWAT_init_mm = 0,
                      u_INTS_init_mm = 0,
                      u_INTR_init_mm = 0,
                      u_SNOW_init_mm = 0,
                      u_CC_init_MJ_per_m2 = 0,
                      u_SNOWLQ_init_mm =  0),
            d18O   = (u_GWAT_init_permil = -13.,
                      u_INTS_init_permil = -13.,
                      u_INTR_init_permil = -13.,
                      u_SNOW_init_permil = -13.),
            d2H    = (u_GWAT_init_permil = -95.,
                      u_INTS_init_permil = -95.,
                      u_INTR_init_permil = -95.,
                      u_SNOW_init_permil = -95.)),
        storm_durations_h = [4, 4, 4, 4, 8, 8, 8, 4, 4, 4, 4, 4]);
    @githash_time example_result3 = remakeSPAC(
        example_result2,
        soil_horizons = (ths_ = 0.4, Ksat_mmday = 3854.9, alpha_per_m = 7.11, gravel_volFrac = 0.1),
        LAI_rel = (DOY_Bstart = 115,),
        root_distribution = (beta = 0.88, z_rootMax_m = -0.6,),
        params = (
                # test flow processes INFL, BYFL, SRFL DSFLI, SEEP/GWFL
                IDEPTH_m=0.50, INFEXP=0.33, 
                BYPAR=1, QDEPTH_m = 0.70, QFFC = 0.2, QFPAR = 0.3,
                IMPERV = 0.01,
                DSLOPE = 10,
                DRAIN=.33, 
                GSC = 0.05, GSP = 0.2,
            ALB=0.15, ALBSN=0.7, RSSA=720., PSICR=-1.6, FXYLEM=0.4, MXKPL=16.5, MAXLAI=9.999,
            GLMAX=.00801, R5=235., CVPD=1.9, CINTRL=0.18,))
    simulate!(example_result3)

    # extract required data from solution object
    # extract specificayll flows: INFL-BYFL-SRFL-DSFLI-VRFLN
    ## helper quantities
    # t_out = range(extrema(example_result3.ODESolution.t)..., step = 30)
    t_out = range(extrema(example_result3.ODESolution.t)..., step = 1)
    t_ref = example_result3.ODESolution.prob.p.REFERENCE_DATE
    d_out = RelativeDaysFloat2DateTime.(t_out, t_ref);
    ## scalar quantities
    # keys(example_result3.ODESolution(1).accum)
    simulated_fluxes = (
        cum_d_prec      = [example_result3.ODESolution(t_days).accum.cum_d_prec      for t_days = t_out],
        cum_d_rfal      = [example_result3.ODESolution(t_days).accum.cum_d_rfal      for t_days = t_out],
        cum_d_sfal      = [example_result3.ODESolution(t_days).accum.cum_d_sfal      for t_days = t_out],
        cum_d_rint      = [example_result3.ODESolution(t_days).accum.cum_d_rint      for t_days = t_out],
        cum_d_sint      = [example_result3.ODESolution(t_days).accum.cum_d_sint      for t_days = t_out],
        cum_d_rthr      = [example_result3.ODESolution(t_days).accum.cum_d_rthr      for t_days = t_out],
        cum_d_sthr      = [example_result3.ODESolution(t_days).accum.cum_d_sthr      for t_days = t_out],
        cum_d_rsno      = [example_result3.ODESolution(t_days).accum.cum_d_rsno      for t_days = t_out],
        cum_d_rnet      = [example_result3.ODESolution(t_days).accum.cum_d_rnet      for t_days = t_out],
        cum_d_smlt      = [example_result3.ODESolution(t_days).accum.cum_d_smlt      for t_days = t_out],
        cum_d_evap      = [example_result3.ODESolution(t_days).accum.cum_d_evap      for t_days = t_out],
        cum_d_tran      = [example_result3.ODESolution(t_days).accum.cum_d_tran      for t_days = t_out],
        cum_d_irvp      = [example_result3.ODESolution(t_days).accum.cum_d_irvp      for t_days = t_out],
        cum_d_isvp      = [example_result3.ODESolution(t_days).accum.cum_d_isvp      for t_days = t_out],
        cum_d_slvp      = [example_result3.ODESolution(t_days).accum.cum_d_slvp      for t_days = t_out],
        cum_d_snvp      = [example_result3.ODESolution(t_days).accum.cum_d_snvp      for t_days = t_out],
        cum_d_pint      = [example_result3.ODESolution(t_days).accum.cum_d_pint      for t_days = t_out],
        cum_d_ptran     = [example_result3.ODESolution(t_days).accum.cum_d_ptran     for t_days = t_out],
        # cum_d_pslvp     = [example_result3.ODESolution(t_days).accum.cum_d_pslvp     for t_days = t_out],
        flow            = [example_result3.ODESolution(t_days).accum.flow            for t_days = t_out],
        seep            = [example_result3.ODESolution(t_days).accum.seep            for t_days = t_out],
        srfl            = [example_result3.ODESolution(t_days).accum.srfl            for t_days = t_out],
        slfl            = [example_result3.ODESolution(t_days).accum.slfl            for t_days = t_out],
        byfl            = [example_result3.ODESolution(t_days).accum.byfl            for t_days = t_out],
        dsfl            = [example_result3.ODESolution(t_days).accum.dsfl            for t_days = t_out],
        gwfl            = [example_result3.ODESolution(t_days).accum.gwfl            for t_days = t_out],
        vrfln           = [example_result3.ODESolution(t_days).accum.vrfln           for t_days = t_out],
        totalSWAT       = [example_result3.ODESolution(t_days).accum.totalSWAT       for t_days = t_out],
        new_totalWATER  = [example_result3.ODESolution(t_days).accum.new_totalWATER  for t_days = t_out],
        BALERD_SWAT     = [example_result3.ODESolution(t_days).accum.BALERD_SWAT     for t_days = t_out],
        BALERD_total    = [example_result3.ODESolution(t_days).accum.BALERD_total    for t_days = t_out],
        )
    
    # test or overwrite
    # fname = input_path*input_prefix
    fname = "../examples/DAV2020-full-modified_fluxes_reference.jld2"
    if task == "test"
        loaded = load(fname)
        @test isapprox(loaded["cum_d_prec"], simulated_fluxes.cum_d_prec, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_rfal"], simulated_fluxes.cum_d_rfal, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_sfal"], simulated_fluxes.cum_d_sfal, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_rint"], simulated_fluxes.cum_d_rint, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_sint"], simulated_fluxes.cum_d_sint, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_rthr"], simulated_fluxes.cum_d_rthr, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_sthr"], simulated_fluxes.cum_d_sthr, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_rsno"], simulated_fluxes.cum_d_rsno, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_rnet"], simulated_fluxes.cum_d_rnet, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_smlt"], simulated_fluxes.cum_d_smlt, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_evap"], simulated_fluxes.cum_d_evap, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_tran"], simulated_fluxes.cum_d_tran, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_irvp"], simulated_fluxes.cum_d_irvp, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_isvp"], simulated_fluxes.cum_d_isvp, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_slvp"], simulated_fluxes.cum_d_slvp, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_snvp"], simulated_fluxes.cum_d_snvp, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_pint"], simulated_fluxes.cum_d_pint, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["cum_d_ptran"], simulated_fluxes.cum_d_ptran, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        # @test isapprox(loaded["cum_d_pslvp"], simulated_fluxes.cum_d_pslvp, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["flow"], simulated_fluxes.flow, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["seep"], simulated_fluxes.seep, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["srfl"], simulated_fluxes.srfl, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["slfl"], simulated_fluxes.slfl, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["byfl"], simulated_fluxes.byfl, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["dsfl"], simulated_fluxes.dsfl, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["gwfl"], simulated_fluxes.gwfl, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["vrfln"], simulated_fluxes.vrfln, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["totalSWAT"], simulated_fluxes.totalSWAT, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["new_totalWATER"], simulated_fluxes.new_totalWATER, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["BALERD_SWAT"], simulated_fluxes.BALERD_SWAT, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
        @test isapprox(loaded["BALERD_total"], simulated_fluxes.BALERD_total, atol = 1e-4) # TODO: somehow we need atol. Why does a simple ≈ not work?
    elseif task == "overwrite" && !is_a_CI_system # only overwrite on local machine, never on CI
        # overwrite output
        jldsave(fname, compress=false; (simulated_fluxes   = simulated_fluxes)...)
        error("Test overwrites reference solution instead of checking against it.")
    else
        # do nothing
    end



    if !is_a_CI_system && plot_flag
        fname_illustrations = "out/$git_status_string/TESTSET_DAV2020modified-regressionTest-FLUXES"
        mkpath(dirname(fname_illustrations))

        pl_fluxes = plot(
            plot(d_out, [simulated_fluxes.cum_d_prec      loaded["cum_d_prec"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_prec", ),
            plot(d_out, [simulated_fluxes.cum_d_rfal      loaded["cum_d_rfal"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_rfal", ),
            plot(d_out, [simulated_fluxes.cum_d_sfal      loaded["cum_d_sfal"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_sfal", ),
            plot(d_out, [simulated_fluxes.cum_d_rint      loaded["cum_d_rint"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_rint", ),
            plot(d_out, [simulated_fluxes.cum_d_sint      loaded["cum_d_sint"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_sint", ),
            plot(d_out, [simulated_fluxes.cum_d_rthr      loaded["cum_d_rthr"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_rthr", ),
            plot(d_out, [simulated_fluxes.cum_d_sthr      loaded["cum_d_sthr"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_sthr", ),
            plot(d_out, [simulated_fluxes.cum_d_rsno      loaded["cum_d_rsno"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_rsno", ),
            plot(d_out, [simulated_fluxes.cum_d_rnet      loaded["cum_d_rnet"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_rnet", ),
            plot(d_out, [simulated_fluxes.cum_d_smlt      loaded["cum_d_smlt"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_smlt", ),
            plot(d_out, [simulated_fluxes.cum_d_evap      loaded["cum_d_evap"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_evap", ),
            plot(d_out, [simulated_fluxes.cum_d_tran      loaded["cum_d_tran"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_tran", ),
            plot(d_out, [simulated_fluxes.cum_d_irvp      loaded["cum_d_irvp"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_irvp", ),
            plot(d_out, [simulated_fluxes.cum_d_isvp      loaded["cum_d_isvp"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_isvp", ),
            plot(d_out, [simulated_fluxes.cum_d_slvp      loaded["cum_d_slvp"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_slvp", ),
            plot(d_out, [simulated_fluxes.cum_d_snvp      loaded["cum_d_snvp"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_snvp", ),
            plot(d_out, [simulated_fluxes.cum_d_pint      loaded["cum_d_pint"]],     legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_pint", ),
            plot(d_out, [simulated_fluxes.cum_d_ptran     loaded["cum_d_ptran"]],    legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_ptran", ),
            # plot(d_out, [simulated_fluxes.cum_d_pslvp     loaded["cum_d_pslvp"]],    legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.cum_d_pslvp", ), # PSLVP was removed from output,
            plot(d_out, [simulated_fluxes.flow            loaded["flow"]],           legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.flow", ),
            plot(d_out, [simulated_fluxes.seep            loaded["seep"]],           legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.seep", ),
            plot(d_out, [simulated_fluxes.srfl            loaded["srfl"]],           legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.srfl", ),
            plot(d_out, [simulated_fluxes.slfl            loaded["slfl"]],           legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.slfl", ),
            plot(d_out, [simulated_fluxes.byfl            loaded["byfl"]],           legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.byfl", ),
            plot(d_out, [simulated_fluxes.dsfl            loaded["dsfl"]],           legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.dsfl", ),
            plot(d_out, [simulated_fluxes.gwfl            loaded["gwfl"]],           legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.gwfl", ),
            plot(d_out, [simulated_fluxes.vrfln           loaded["vrfln"]],          legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.vrfln", ),
            plot(d_out, [simulated_fluxes.totalSWAT       loaded["totalSWAT"]],      legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.totalSWAT", ),
            plot(d_out, [simulated_fluxes.new_totalWATER  loaded["new_totalWATER"]], legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.new_totalWATER", ),
            plot(d_out, [simulated_fluxes.BALERD_SWAT     loaded["BALERD_SWAT"]],    legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.BALERD_SWAT", ),
            plot(d_out, [simulated_fluxes.BALERD_total    loaded["BALERD_total"]],   legend = :topright, linestyle = [:solid :dot], label = ["current code" "reference simulation"], title = "accum.BALERD_total", ),
            size=(2000,1000), layout = (7,5)
            )
        savefig(plot(pl_fluxes, size=(2000,1400), dpi=300),  fname_illustrations*"_fluxes.png")
    end

end