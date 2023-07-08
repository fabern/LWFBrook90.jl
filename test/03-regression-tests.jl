# Regression tests
using JLD2: jldsave, load # TODO: comment in again as soon  https://github.com/JuliaIO/JLD2.jl/issues/428 is fixed
# using FileIO
import FileIO

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
    # TODO: replace get_δsoil(...) by get_soil_(:δ18O, ...)
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
        error("Test overwrites reference solution instead of checking against it.")
    else
        # do nothing
    end


    if !is_a_CI_system && plot_flag
        fname_illustrations = "out/$git_status_string/TESTSET_DAV2020-regressionTest"
        mkpath(dirname(fname_illustrations))

        pl1 = plotamounts(example_result, :above_and_belowground, :showRWUcentroid)
        pl2 = plotisotopes(example_result, :d18O, (d18O = (-16, -6), d2H = (-125, -40)), :showRWUcentroid)
        pl3 = plotforcingandstates(example_result)
        savefig(plot(pl1, size=(1000,1400), dpi=300),  fname_illustrations*"_amts.png")
        savefig(plot(pl2, size=(1000,700), dpi=300), fname_illustrations*"_d18O-d2H.png")
        savefig(plot(pl3, size=(1000,1400), dpi=300), fname_illustrations*"_forcing.png")
    end

end

function get_daily_soilFluxes(simulation)
    t_out = range(extrema(simulation.ODESolution.t)..., step = 1)
    t_ref = simulation.ODESolution.prob.p.REFERENCE_DATE
    d_out = RelativeDaysFloat2DateTime.(t_out, t_ref);

    # extract specifically flows: INFL-BYFL-SRFL-DSFLI-VRFLN
    # keys(simulation.ODESolution(1).accum)
    simulated_fluxes = (
        cum_d_prec      = [simulation.ODESolution(t_days).accum.cum_d_prec      for t_days = t_out],
        cum_d_rfal      = [simulation.ODESolution(t_days).accum.cum_d_rfal      for t_days = t_out],
        cum_d_sfal      = [simulation.ODESolution(t_days).accum.cum_d_sfal      for t_days = t_out],
        cum_d_rint      = [simulation.ODESolution(t_days).accum.cum_d_rint      for t_days = t_out],
        cum_d_sint      = [simulation.ODESolution(t_days).accum.cum_d_sint      for t_days = t_out],
        cum_d_rthr      = [simulation.ODESolution(t_days).accum.cum_d_rthr      for t_days = t_out],
        cum_d_sthr      = [simulation.ODESolution(t_days).accum.cum_d_sthr      for t_days = t_out],
        cum_d_rsno      = [simulation.ODESolution(t_days).accum.cum_d_rsno      for t_days = t_out],
        cum_d_rnet      = [simulation.ODESolution(t_days).accum.cum_d_rnet      for t_days = t_out],
        cum_d_smlt      = [simulation.ODESolution(t_days).accum.cum_d_smlt      for t_days = t_out],
        cum_d_evap      = [simulation.ODESolution(t_days).accum.cum_d_evap      for t_days = t_out],
        cum_d_tran      = [simulation.ODESolution(t_days).accum.cum_d_tran      for t_days = t_out],
        cum_d_irvp      = [simulation.ODESolution(t_days).accum.cum_d_irvp      for t_days = t_out],
        cum_d_isvp      = [simulation.ODESolution(t_days).accum.cum_d_isvp      for t_days = t_out],
        cum_d_slvp      = [simulation.ODESolution(t_days).accum.cum_d_slvp      for t_days = t_out],
        cum_d_snvp      = [simulation.ODESolution(t_days).accum.cum_d_snvp      for t_days = t_out],
        cum_d_pint      = [simulation.ODESolution(t_days).accum.cum_d_pint      for t_days = t_out],
        cum_d_ptran     = [simulation.ODESolution(t_days).accum.cum_d_ptran     for t_days = t_out],
        # cum_d_pslvp     = [simulation.ODESolution(t_days).accum.cum_d_pslvp     for t_days = t_out],
        flow            = [simulation.ODESolution(t_days).accum.flow            for t_days = t_out],
        seep            = [simulation.ODESolution(t_days).accum.seep            for t_days = t_out],
        srfl            = [simulation.ODESolution(t_days).accum.srfl            for t_days = t_out],
        slfl            = [simulation.ODESolution(t_days).accum.slfl            for t_days = t_out],
        byfl            = [simulation.ODESolution(t_days).accum.byfl            for t_days = t_out],
        dsfl            = [simulation.ODESolution(t_days).accum.dsfl            for t_days = t_out],
        gwfl            = [simulation.ODESolution(t_days).accum.gwfl            for t_days = t_out],
        vrfln           = [simulation.ODESolution(t_days).accum.vrfln           for t_days = t_out],
        StorageSWAT     = [simulation.ODESolution(t_days).accum.StorageSWAT     for t_days = t_out],
        StorageWATER    = [simulation.ODESolution(t_days).accum.StorageWATER    for t_days = t_out],
        BALERD_SWAT     = [simulation.ODESolution(t_days).accum.BALERD_SWAT     for t_days = t_out],
        BALERD_total    = [simulation.ODESolution(t_days).accum.BALERD_total    for t_days = t_out],
        )
    return d_out, simulated_fluxes
end
function plot_simulated_fluxes_vs_reference(simulated_fluxes, reference, d_out; labels = ["current code" "reference simulation"], kwargs...)
    plot(
        plot(d_out, [simulated_fluxes.cum_d_prec      reference["cum_d_prec"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_prec", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_rfal      reference["cum_d_rfal"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_rfal", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_sfal      reference["cum_d_sfal"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_sfal", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_rint      reference["cum_d_rint"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_rint", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_sint      reference["cum_d_sint"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_sint", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_rthr      reference["cum_d_rthr"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_rthr", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_sthr      reference["cum_d_sthr"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_sthr", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_rsno      reference["cum_d_rsno"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_rsno", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_rnet      reference["cum_d_rnet"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_rnet", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_smlt      reference["cum_d_smlt"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_smlt", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_evap      reference["cum_d_evap"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_evap", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_tran      reference["cum_d_tran"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_tran", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_irvp      reference["cum_d_irvp"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_irvp", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_isvp      reference["cum_d_isvp"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_isvp", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_slvp      reference["cum_d_slvp"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_slvp", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_snvp      reference["cum_d_snvp"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_snvp", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_pint      reference["cum_d_pint"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_pint", kwargs...),
        plot(d_out, [simulated_fluxes.cum_d_ptran     reference["cum_d_ptran"]],    linestyle = [:solid :dot], label = labels, title = "accum.cum_d_ptran", kwargs...),
        # plot(d_out, [simulated_fluxes.cum_d_pslvp     reference["cum_d_pslvp"]],    linestyle = [:solid :dot], label = labels, title = "accum.cum_d_pslvp", kwargs...), # PSLVP was removed from output,
        plot(d_out, [simulated_fluxes.flow            reference["flow"]],           linestyle = [:solid :dot], label = labels, title = "accum.flow", kwargs...),
        plot(d_out, [simulated_fluxes.seep            reference["seep"]],           linestyle = [:solid :dot], label = labels, title = "accum.seep", kwargs...),
        plot(d_out, [simulated_fluxes.srfl            reference["srfl"]],           linestyle = [:solid :dot], label = labels, title = "accum.srfl", kwargs...),
        plot(d_out, [simulated_fluxes.slfl            reference["slfl"]],           linestyle = [:solid :dot], label = labels, title = "accum.slfl", kwargs...),
        plot(d_out, [simulated_fluxes.byfl            reference["byfl"]],           linestyle = [:solid :dot], label = labels, title = "accum.byfl", kwargs...),
        plot(d_out, [simulated_fluxes.dsfl            reference["dsfl"]],           linestyle = [:solid :dot], label = labels, title = "accum.dsfl", kwargs...),
        plot(d_out, [simulated_fluxes.gwfl            reference["gwfl"]],           linestyle = [:solid :dot], label = labels, title = "accum.gwfl", kwargs...),
        plot(d_out, [simulated_fluxes.vrfln           reference["vrfln"]],          linestyle = [:solid :dot], label = labels, title = "accum.vrfln", kwargs...),
        plot(d_out, [simulated_fluxes.StorageSWAT     reference["StorageSWAT"]],    linestyle = [:solid :dot], label = labels, title = "accum.StorageSWAT", kwargs...),
        plot(d_out, [simulated_fluxes.StorageWATER    reference["StorageWATER"]],   linestyle = [:solid :dot], label = labels, title = "accum.StorageWATER", kwargs...),
        plot(d_out, [simulated_fluxes.BALERD_SWAT     reference["BALERD_SWAT"]],    linestyle = [:solid :dot], label = labels, title = "accum.BALERD_SWAT", kwargs...),
        plot(d_out, [simulated_fluxes.BALERD_total    reference["BALERD_total"]],   linestyle = [:solid :dot], label = labels, title = "accum.BALERD_total", kwargs...)
    )
end
function test_fluxes_comparison(simulated_fluxes, reference)
        @test isapprox(reference["cum_d_prec"],     simulated_fluxes.cum_d_prec,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_rfal"],     simulated_fluxes.cum_d_rfal,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_sfal"],     simulated_fluxes.cum_d_sfal,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_rint"],     simulated_fluxes.cum_d_rint,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_sint"],     simulated_fluxes.cum_d_sint,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_rthr"],     simulated_fluxes.cum_d_rthr,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_sthr"],     simulated_fluxes.cum_d_sthr,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_rsno"],     simulated_fluxes.cum_d_rsno,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_rnet"],     simulated_fluxes.cum_d_rnet,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_smlt"],     simulated_fluxes.cum_d_smlt,     atol = 1e-4, rtol = 1e-4)
        @test_skip isapprox(reference["cum_d_evap"],     simulated_fluxes.cum_d_evap,     atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
        @test_skip isapprox(reference["cum_d_tran"],     simulated_fluxes.cum_d_tran,     atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
        @test isapprox(reference["cum_d_irvp"],     simulated_fluxes.cum_d_irvp,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_isvp"],     simulated_fluxes.cum_d_isvp,     atol = 1e-4, rtol = 1e-4)
        @test_skip isapprox(reference["cum_d_slvp"],     simulated_fluxes.cum_d_slvp,     atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
        @test isapprox(reference["cum_d_snvp"],     simulated_fluxes.cum_d_snvp,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_pint"],     simulated_fluxes.cum_d_pint,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_ptran"],    simulated_fluxes.cum_d_ptran,    atol = 1e-4, rtol = 1e-4)
        # @test isapprox(reference["cum_d_pslvp"],  simulated_fluxes.cum_d_pslvp,    atol = 1e-4, rtol = 1e-4)
        @test_skip isapprox(reference["flow"],           simulated_fluxes.flow,           atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
        @test isapprox(reference["seep"],           simulated_fluxes.seep,           atol = 1e-4, rtol = 1e-4)
        @test_skip isapprox(reference["srfl"],           simulated_fluxes.srfl,           atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
        @test_skip isapprox(reference["slfl"],           simulated_fluxes.slfl,           atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
        @test_skip isapprox(reference["byfl"],           simulated_fluxes.byfl,           atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
        @test isapprox(reference["dsfl"],           simulated_fluxes.dsfl,           atol = 1e-4, rtol = 1e-4) #
        @test_skip isapprox(reference["gwfl"],           simulated_fluxes.gwfl,           atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
        @test_skip isapprox(reference["vrfln"],          simulated_fluxes.vrfln,          atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
        @test isapprox(reference["StorageSWAT"],    simulated_fluxes.StorageSWAT,    atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["StorageWATER"],   simulated_fluxes.StorageWATER,   atol = 1e-4, rtol = 1e-4)
        @test_skip isapprox(reference["BALERD_SWAT"],    simulated_fluxes.BALERD_SWAT,    atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
        @test_skip isapprox(reference["BALERD_total"],   simulated_fluxes.BALERD_total,   atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
end

@testset "Oversaturation-infiltration-FLUXES" begin
    model2 = loadSPAC("../examples/infiltrationSaturationINFEXP1/", "infiltrationSaturationINFEXP1";
        simulate_isotopes = false,
        Δz_thickness_m = [fill(0.05, 42);],
        IC_soil = (PSIM_init_kPa = -6.3, delta18O_init_permil = -13.0, delta2H_init_permil = -95.0),
        root_distribution = (beta = 0.90, z_rootMax_m = nothing));
    simulation2  = setup(model2; requested_tspan = (0,300));
    simulate!(simulation2);
    simulation2_withVariousFlows = remakeSPAC(simulation2, params = (
        # activate all flow processes INFL, BYFL, SRFL DSFLI, SEEP/GWFL
        BYPAR=1, QDEPTH_m = 0.40, QFFC = 0.2, QFPAR = 0.3,
        # IMPERV = 0.01,
        # DSLOPE = 10,
        # DRAIN=.33,
        GSC = 0.05, GSP = 0.2));
    simulate!(simulation2_withVariousFlows);

    d_out2,  simulated_fluxes2  = get_daily_soilFluxes(simulation2);
    d_out2b, simulated_fluxes2b = get_daily_soilFluxes(simulation2_withVariousFlows);

    fname2  = "../examples/INFEXP1_fluxes_reference.jld2"
    fname2b = "../examples/INFEXP1-modified_fluxes_reference.jld2"
    if task == "test"
        loaded2 = load(fname2)
        loaded2b = load(fname2b)

        test_fluxes_comparison(simulated_fluxes2, loaded2)
        test_fluxes_comparison(simulated_fluxes2b, loaded2b)
    elseif task == "overwrite" && !is_a_CI_system # only overwrite on local machine, never on CI
        # overwrite output
        jldsave(fname2, compress=false; (simulated_fluxes2   = simulated_fluxes2)..., d_out2 = d_out2)
        jldsave(fname2b, compress=false; (simulated_fluxes2b   = simulated_fluxes2b)..., d_out2b = d_out2b)
        error("Test overwrites reference solution instead of checking against it.")
    else
        # do nothing
    end

    if !is_a_CI_system && plot_flag
        fname_illustrations = "out/$git_status_string/TESTSET_Oversaturation-infiltration-FLUXES"
        mkpath(dirname(fname_illustrations))

        pl_fluxes2 = plot_simulated_fluxes_vs_reference(simulated_fluxes2, loaded2, loaded2["d_out2"]);
        plot!(legend = :topleft, size=(2000,1000), layout = (7,5))
        savefig(plot(pl_fluxes2, size=(2000,1000), dpi=300),  fname_illustrations*"_fluxes_noBYFL_regressionTest.png")

        pl_fluxes2b = plot_simulated_fluxes_vs_reference(simulated_fluxes2b, loaded2b, loaded2b["d_out2b"]);
        plot!(legend = :topleft, size=(2000,1000), layout = (7,5))
        savefig(plot(pl_fluxes2b, size=(2000,1000), dpi=300),  fname_illustrations*"_fluxes_withBYFL_regressionTest.png")

        pl_fluxes2vs2b = plot_simulated_fluxes_vs_reference(
            simulated_fluxes2,
            Dict(String(k) => v for (k,v) in pairs(simulated_fluxes2b)),
            d_out2b,
            labels = ["Without BYFL" "With BYFL"]);
        plot!(legend = :topleft, size=(2000,1000), layout = (7,5))
        savefig(plot(pl_fluxes2vs2b, size=(2000,1000), dpi=300),  fname_illustrations*"_fluxes_BYFL_comparison.png")
        pl_comparison = plot(plotamounts(simulation2, title = "Without BYFL"),
            plotamounts(simulation2_withVariousFlows, title = "With BYFL"),
            layout = (1,2), size=(1460,1400))
        savefig(plot(pl_comparison, size=(2000,1000), dpi=300),  fname_illustrations*"_fluxes_BYFL_comparison2.png")
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
    d_out, simulated_fluxes = get_daily_soilFluxes(example_result3);

    # test or overwrite
    # fname = input_path*input_prefix
    fname = "../examples/DAV2020-full-modified_fluxes_reference.jld2"
    if task == "test"
        loaded = load(fname)
        test_fluxes_comparison(simulated_fluxes, loaded)
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

        pl_fluxes = plot_simulated_fluxes_vs_reference(simulated_fluxes, loaded, d_out);
        plot!(legend = :topright, size=(2000,1000), layout = (7,5))
        savefig(plot(pl_fluxes, size=(2000,1400), layout = (7,5), dpi=300),  fname_illustrations*"_fluxes.png")
    end

end