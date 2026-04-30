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

function get_some_states_to_compare(example_result)
    ## helper quantities
    t_out = range(extrema(example_result.ODESolution.t)..., step = 30)
    ## scalar quantities
    df_aboveground_names = ["GWAT" "INTS" "INTR" "SNOW" "CC" "SNOWLQ" "RWU" "XYLEM"]
    mat_aboveground = vcat(
        reduce(hcat, [example_result.ODESolution(t_days).GWAT.mm   for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).INTS.mm   for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).INTR.mm   for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).SNOW.mm   for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).CC.MJm2   for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).SNOWLQ.mm for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).RWU.mmday for t_days = t_out]),
        reduce(hcat, [example_result.ODESolution(t_days).XYLEM.mm  for t_days = t_out]),
        # reduce(hcat, [example_result.ODESolution.u[t_idx].TRANI.mmday for t_idx = eachindex(example_result.ODESolution)]),
        # reduce(hcat, [example_result.ODESolution.u[t_idx].aux.θ       for t_idx = eachindex(example_result.ODESolution)]),
        # reduce(hcat, [example_result.ODESolution.u[t_idx].accum       for t_idx = eachindex(example_result.ODESolution)])
    )
    u_mm = hcat(DataFrame(time = t_out), DataFrame(permutedims(mat_aboveground), df_aboveground_names[:]))
    u_δ     = innerjoin(
        get_states(example_result, days_to_read_out_d = t_out),
        get_fluxes(example_result, days_to_read_out_d = t_out),
        on = [:dates])[:,[:dates,:PREC_d18O,:GWAT_d18O,:INTS_d18O,:INTR_d18O,:SNOW_d18O,:RWU_d18O,:XYLEM_d18O,:PREC_d2H,:GWAT_d2H,:INTS_d2H,:INTR_d2H,:SNOW_d2H,:RWU_d2H,:XYLEM_d2H]]
    ## vector quantities
    u_belowground = get_soil_([:SWATI, :W, :ψ, :θ, :K, :δ18O, :δ2H], example_result;
        days_to_read_out_d = t_out)
    u_belowground = u_belowground[:,Not(r"(K_)|(W_)")]

    return u_mm, u_δ, u_belowground
    # return t_out, u_ref_names, mat_aboveground, u_δ, u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_θ, p_fu_KK, u_δsoil_d18O, u_δsoil_d2H
end

function plot_simulated_states_vs_reference(
    u_mm, u_δ, u_belowground,
    loaded_u_mm, loaded_u_δ, loaded_u_belowground)
    # t_out = u_mm.time
    u_ref_names = names(u_mm[:,Not(:time)])
    u_δ = u_δ
    u_ref         = permutedims(Matrix(u_mm[:,Not(:time)]))
    u_SWATI       = Matrix(permutedims(select(u_belowground, r"SWATI_")));
    # u_aux_WETNES  = Matrix(permutedims(select(u_belowground, r"W_")));
    u_aux_PSIM    = Matrix(permutedims(select(u_belowground, r"ψ_")));
    u_aux_θ       = Matrix(permutedims(select(u_belowground, r"θ_")));
    # p_fu_KK       = Matrix(permutedims(select(u_belowground, r"K_")));
    u_δsoil_d18O  = select(u_belowground, r"δ18O_") # loaded["u_δsoil"].d18O
    u_δsoil_d2H   = select(u_belowground, r"δ2H_" ) # loaded["u_δsoil"].d2H

    loaded_u_ref         = permutedims(Matrix(loaded_u_mm[:,Not(:time)]))
    loaded_u_SWATI       = Matrix(permutedims(select(loaded_u_belowground, r"SWATI_")));
    # loaded_u_aux_WETNES  = Matrix(permutedims(select(loaded_u_belowground, r"W_")));
    loaded_u_aux_PSIM    = Matrix(permutedims(select(loaded_u_belowground, r"ψ_")));
    loaded_u_aux_θ       = Matrix(permutedims(select(loaded_u_belowground, r"θ_")));
    # loaded_p_fu_KK       = Matrix(permutedims(select(loaded_u_belowground, r"K_")));
    loaded_u_δsoil_d18O  = select(loaded_u_belowground, r"δ18O_") # loaded["u_δsoil"].d18O
    loaded_u_δsoil_d2H   = select(loaded_u_belowground, r"δ2H_" ) # loaded["u_δsoil"].d2H

    d_out = u_belowground.time

    function plot_comparison(A, B_ref; d_out=1:size(A,1), nams = "")
        labels = if (nams=="")
            try names(A);
            catch
                fill("", size(A,2))
            end
        else
            nams
        end
        pl = Plots.plot(d_out, Matrix(A), label = reshape(labels, 1, :), linestyle = :solid)
        [Plots.plot!(pl, d_out, col, label = ifelse(it==1, "reference simulation",""), linestyle = :dash, linecolor = :black)
            for (it, col) in pairs(eachcol(Matrix(B_ref)))]
        return pl
    end

    Plots.plot(
            plot_comparison(u_ref',          loaded_u_ref',          d_out=d_out, nams = u_ref_names),
            plot_comparison(u_δ[:, r"d18O"], loaded_u_δ[:, r"d18O"], d_out=d_out),
            plot_comparison(u_δ[:, r"d2H"],  loaded_u_δ[:, r"d2H"],  d_out=d_out),

            plot_comparison(u_SWATI',        loaded_u_SWATI',        d_out=d_out, nams = ["SWATI_$l" for _ in [1], l in (1:5)]),
            plot_comparison(u_aux_PSIM',     loaded_u_aux_PSIM',     d_out=d_out, nams = ["PSIM_$l"  for _ in [1], l in (1:5)]),
            plot_comparison(u_aux_θ',        loaded_u_aux_θ',        d_out=d_out, nams = ["θ_$l"     for _ in [1], l in (1:5)]),
            plot_comparison(u_δsoil_d18O,    loaded_u_δsoil_d18O, d_out=d_out, nams = names(u_δsoil_d18O)),
            plot_comparison(u_δsoil_d2H,     loaded_u_δsoil_d2H,  d_out=d_out, nams = names(u_δsoil_d2H)),
            layout = (2,4), size = (1200,800))
end

# Helper: serialise an array of structs to a DataFrame (one row per element, one column per field)
function spac_structs_to_df(arr)
    fields = fieldnames(typeof(arr[1]))
    DataFrame([field => [getfield(s, field) for s in arr] for field in fields]...)
end

# Helper: convert a flat NamedTuple of scalars to a single-row DataFrame.
function nt_to_df(nt::NamedTuple)
    DataFrame([k => [v] for (k, v) in pairs(nt)]...)
end

# Save all SPAC fields (except forcing) to CSV files for regression testing.
# Using CSV removes the JLD2 brittleness when type definitions change.
function save_spac_to_csvs(currSPAC, fname_base)
    # Scalars and strings
    write(fname_base * "_spac_scalars.csv", DataFrame(
        reference_date   = [string(currSPAC.reference_date)],
        tspan_start      = [first(currSPAC.tspan)],
        tspan_end        = [last(currSPAC.tspan)],
        root_distribution = [string(currSPAC.pars.root_distribution)],
        IC_soil          = [currSPAC.pars.IC_soil],
    ))
    # NamedTuples of scalar parameters → single-row DataFrames
    write(fname_base * "_spac_params.csv",         nt_to_df(currSPAC.pars.params))
    write(fname_base * "_spac_solver_options.csv", nt_to_df(currSPAC.solver_options))
    # Fields that are already DataFrames → write directly
    write(fname_base * "_spac_IC_scalar.csv",        currSPAC.pars.IC_scalar)
    write(fname_base * "_spac_canopy_evolution.csv", currSPAC.pars.canopy_evolution)
    # Soil discretization
    write(fname_base * "_spac_soil_discretization.csv",
        currSPAC.soil_discretization.df[:, Not(:shp)])
    write(fname_base * "_spac_Δz.csv", DataFrame(Δz = currSPAC.soil_discretization.Δz))
    write(fname_base * "_spac_soil_horizons_shp.csv",
        spac_structs_to_df(currSPAC.pars.soil_horizons.shp))
    write(fname_base * "_spac_soil_discretization_shp.csv",
        spac_structs_to_df(currSPAC.soil_discretization.df[:, :shp]))
end

# Test all SPAC fields against CSV reference files (robust to type-definition changes).
function test_spac_from_csvs(currSPAC, fname_base)
    @testset "test_spac_comparison" begin
        # Scalars and strings
        loaded_scalars = read(fname_base * "_spac_scalars.csv", DataFrame)
        @test loaded_scalars.reference_date[1]    == currSPAC.reference_date
        @test loaded_scalars.tspan_start[1]       == first(currSPAC.tspan)
        @test loaded_scalars.tspan_end[1]         == last(currSPAC.tspan)
        @test loaded_scalars.root_distribution[1] == string(currSPAC.pars.root_distribution) # e.g. string((beta = 0.97, z_rootMax_m = nothing))
        @test loaded_scalars.IC_soil[1]           == currSPAC.pars.IC_soil

        # NamedTuples of scalar parameters
        @test read(fname_base * "_spac_params.csv",         DataFrame) == nt_to_df(currSPAC.pars.params)
        @test read(fname_base * "_spac_solver_options.csv", DataFrame) == nt_to_df(currSPAC.solver_options)

        # DataFrame fields
        @test read(fname_base * "_spac_IC_scalar.csv",        DataFrame) == currSPAC.pars.IC_scalar
        @test read(fname_base * "_spac_canopy_evolution.csv", DataFrame) == currSPAC.pars.canopy_evolution

        # Soil discretization
        @test read(fname_base * "_spac_soil_discretization.csv", DataFrame) ==
              currSPAC.soil_discretization.df[:, Not(:shp)]
        @test read(fname_base * "_spac_Δz.csv", DataFrame).Δz == currSPAC.soil_discretization.Δz

        loaded_shp = read(fname_base * "_spac_soil_horizons_shp.csv", DataFrame)
        @test loaded_shp == spac_structs_to_df(currSPAC.pars.soil_horizons.shp)

        loaded_discr_shp = read(fname_base * "_spac_soil_discretization_shp.csv", DataFrame)
        @test loaded_discr_shp == spac_structs_to_df(currSPAC.soil_discretization.df[:, :shp])
    end
end

function test_states_comparison(u_mm, u_δ, u_belowground,
    loaded_u_mm, loaded_u_δ, loaded_u_belowground)

    # t_out = u_mm.time
    # u_ref_names = names(u_mm)
    u_δ = u_δ
    u_ref         = permutedims(Matrix(u_mm[:,Not(:time)]))
    u_SWATI       = Matrix(permutedims(select(u_belowground, r"SWATI_")));
    # u_aux_WETNES  = Matrix(permutedims(select(u_belowground, r"W_")));
    u_aux_PSIM    = Matrix(permutedims(select(u_belowground, r"ψ_")));
    u_aux_θ       = Matrix(permutedims(select(u_belowground, r"θ_")));
    # p_fu_KK       = Matrix(permutedims(select(u_belowground, r"K_")));
    u_δsoil_d18O  = select(u_belowground, r"δ18O_") # loaded["u_δsoil"].d18O
    u_δsoil_d2H   = select(u_belowground, r"δ2H_" ) # loaded["u_δsoil"].d2H

    loaded_u_ref         = permutedims(Matrix(loaded_u_mm[:,Not(:time)]))
    loaded_u_SWATI       = Matrix(permutedims(select(loaded_u_belowground, r"SWATI_")));
    # loaded_u_aux_WETNES  = Matrix(permutedims(select(loaded_u_belowground, r"W_")));
    loaded_u_aux_PSIM    = Matrix(permutedims(select(loaded_u_belowground, r"ψ_")));
    loaded_u_aux_θ       = Matrix(permutedims(select(loaded_u_belowground, r"θ_")));
    # loaded_p_fu_KK       = Matrix(permutedims(select(loaded_u_belowground, r"K_")));
    loaded_u_δsoil_d18O  = select(loaded_u_belowground, r"δ18O_") # loaded["u_δsoil"].d18O
    loaded_u_δsoil_d2H   = select(loaded_u_belowground, r"δ2H_" ) # loaded["u_δsoil"].d2H

    @testset "test_states_comparison" begin
        # Test scalar states
        compare_scalar = (A,B; nans = false) -> all(isapprox.(Matrix(A), Matrix(B); nans))
        @test compare_scalar(loaded_u_ref, u_ref)
        @test compare_scalar(loaded_u_δ[:, Not(1)],   u_δ[:, Not(1)]; nans=true) # Not(1) removes :time or :dates

        # Test vector states
        compare_vector = (A,B) -> all(isapprox.(Matrix(A), Matrix(B)))
        @test compare_vector(loaded_u_SWATI,       u_SWATI)
        @test compare_vector(loaded_u_aux_PSIM,    u_aux_PSIM)
        @test compare_vector(loaded_u_aux_θ,       u_aux_θ)
        @test compare_vector(loaded_u_δsoil_d18O, u_δsoil_d18O)
        @test compare_vector(loaded_u_δsoil_d2H,  u_δsoil_d2H)
    end
end

@testset "regression-delta-run_example" begin
    @githash_time example_result = LWFBrook90.run_example();

    # extract required data from solution object
    u_mm, u_δ, u_belowground = get_some_states_to_compare(example_result);

    # test or overwrite
    fname = "../examples/DAV2020-full_u_sol_reference.jld2"
    loaded_u_mm          = read(replace(fname, ".jld2"=>"u_mm.csv"), DataFrame)
    loaded_u_δ           = read(replace(fname, ".jld2"=>"u_δ.csv"), DataFrame)
    loaded_u_belowground = read(replace(fname, ".jld2"=>"u_belowground.csv"), DataFrame)

    currSPAC = example_result.parametrizedSPAC;
    fname_base = replace(fname, ".jld2" => "")
    if task == "test"
        test_states_comparison(u_mm, u_δ, u_belowground,
            loaded_u_mm, loaded_u_δ, loaded_u_belowground)
        test_spac_from_csvs(currSPAC, fname_base)
    elseif task == "overwrite" && !is_a_CI_system # only overwrite on local machine, never on CI
        # overwrite output
        fname_illustrations = "out/$git_status_string/TESTSET_DAV2020-regressionTest"
        mkpath(dirname(fname_illustrations))
        pl_comparison = plot_simulated_states_vs_reference(
            u_mm, u_δ, u_belowground,
            loaded_u_mm, loaded_u_δ, loaded_u_belowground)

        savefig(pl_comparison,fname_illustrations*"_regression_test_overwritten.png")
        save_spac_to_csvs(currSPAC, fname_base)
        write(replace(fname, ".jld2"=>"u_mm.csv"), u_mm)
        write(replace(fname, ".jld2"=>"u_δ.csv"), u_δ)
        write(replace(fname, ".jld2"=>"u_belowground.csv"), u_belowground)
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
        savefig(Plots.plot(pl1, size=(1000,1400), dpi=300),  fname_illustrations*"_amts.png")
        savefig(Plots.plot(pl2, size=(1000,700), dpi=300), fname_illustrations*"_d18O-d2H.png")
        savefig(Plots.plot(pl3, size=(1000,1400), dpi=300), fname_illustrations*"_forcing.png")
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
        cum_d_tran      = [simulation.ODESolution(t_days).accum.cum_d_tran      for t_days = t_out],
        cum_d_irvp      = [simulation.ODESolution(t_days).accum.cum_d_irvp      for t_days = t_out],
        cum_d_isvp      = [simulation.ODESolution(t_days).accum.cum_d_isvp      for t_days = t_out],
        cum_d_slvp      = [simulation.ODESolution(t_days).accum.cum_d_slvp      for t_days = t_out],
        cum_d_snvp      = [simulation.ODESolution(t_days).accum.cum_d_snvp      for t_days = t_out],
        cum_d_pint      = [simulation.ODESolution(t_days).accum.cum_d_pint      for t_days = t_out],
        cum_d_ptran     = [simulation.ODESolution(t_days).accum.cum_d_ptran     for t_days = t_out],
        # cum_d_pslvp     = [simulation.ODESolution(t_days).accum.cum_d_pslvp     for t_days = t_out],
        evap            = [simulation.ODESolution(t_days).accum.evap      for t_days = t_out],
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

    return hcat(DataFrame(time = d_out), DataFrame(simulated_fluxes))
end
function get_daily_soilFluxes_new(simulation)
    
    out = simulation.saved_values;
    t_out = out.t;
    t_ref = simulation.ODESolution.prob.p.REFERENCE_DATE;
    d_out = RelativeDaysFloat2DateTime.(t_out, t_ref);

    # format to match order of current fluxes output
    simulated_fluxes = (
        cum_d_prec      = vcat(0, [out.saveval[t_days].accum.cum_d_prec      for t_days = 1:(length(t_out)-1)]),
        cum_d_rfal      = vcat(0, [out.saveval[t_days].accum.cum_d_rfal      for t_days = 1:(length(t_out)-1)]),
        cum_d_sfal      = vcat(0, [out.saveval[t_days].accum.cum_d_sfal      for t_days = 1:(length(t_out)-1)]),
        cum_d_rint      = vcat(0, [out.saveval[t_days].accum.cum_d_rint      for t_days = 1:(length(t_out)-1)]),
        cum_d_sint      = vcat(0, [out.saveval[t_days].accum.cum_d_sint      for t_days = 1:(length(t_out)-1)]),
        cum_d_rthr      = vcat(0, [out.saveval[t_days].accum.cum_d_rthr      for t_days = 1:(length(t_out)-1)]),
        cum_d_sthr      = vcat(0, [out.saveval[t_days].accum.cum_d_sthr      for t_days = 1:(length(t_out)-1)]),
        cum_d_rsno      = vcat(0, [out.saveval[t_days].accum.cum_d_rsno      for t_days = 1:(length(t_out)-1)]),
        cum_d_rnet      = vcat(0, [out.saveval[t_days].accum.cum_d_rnet      for t_days = 1:(length(t_out)-1)]),
        cum_d_smlt      = vcat(0, [out.saveval[t_days].accum.cum_d_smlt      for t_days = 1:(length(t_out)-1)]),
        cum_d_tran      = vcat(0, [out.saveval[t_days].accum.cum_d_tran      for t_days = 1:(length(t_out)-1)]),
        cum_d_irvp      = vcat(0, [out.saveval[t_days].accum.cum_d_irvp      for t_days = 1:(length(t_out)-1)]),
        cum_d_isvp      = vcat(0, [out.saveval[t_days].accum.cum_d_isvp      for t_days = 1:(length(t_out)-1)]),
        cum_d_slvp      = vcat(0, [out.saveval[t_days].accum.cum_d_slvp      for t_days = 1:(length(t_out)-1)]),
        cum_d_snvp      = vcat(0, [out.saveval[t_days].accum.cum_d_snvp      for t_days = 1:(length(t_out)-1)]),
        cum_d_pint      = vcat(0, [out.saveval[t_days].accum.cum_d_pint      for t_days = 1:(length(t_out)-1)]),
        cum_d_ptran     = vcat(0, [out.saveval[t_days].accum.cum_d_ptran     for t_days = 1:(length(t_out)-1)]),
        evap            = vcat(0, [out.saveval[t_days].accum.evap            for t_days = 1:(length(t_out)-1)]),
        flow            = [out.saveval[t_days].accum.flow            for t_days = 1:length(t_out)],
        seep            = [out.saveval[t_days].accum.seep            for t_days = 1:length(t_out)],
        srfl            = [out.saveval[t_days].accum.srfl            for t_days = 1:length(t_out)],
        slfl            = [out.saveval[t_days].accum.slfl            for t_days = 1:length(t_out)],
        byfl            = [out.saveval[t_days].accum.byfl            for t_days = 1:length(t_out)],
        dsfl            = [out.saveval[t_days].accum.dsfl            for t_days = 1:length(t_out)],
        gwfl            = [out.saveval[t_days].accum.gwfl            for t_days = 1:length(t_out)],
        vrfln           = [out.saveval[t_days].accum.vrfln           for t_days = 1:length(t_out)],
        StorageSWAT     = vcat(out.saveval[1].accum.StorageSWAT,  [out.saveval[t_days].accum.StorageSWAT     for t_days = 1:(length(t_out)-1)]),
        StorageWATER    = vcat(out.saveval[1].accum.StorageWATER, [out.saveval[t_days].accum.StorageWATER    for t_days = 1:(length(t_out)-1)]),
        BALERD_SWAT     = vcat(out.saveval[1].accum.BALERD_SWAT,  [out.saveval[t_days].accum.BALERD_SWAT     for t_days = 1:(length(t_out)-1)]),
        BALERD_total    = vcat(out.saveval[1].accum.BALERD_total, [out.saveval[t_days].accum.BALERD_total    for t_days = 1:(length(t_out)-1)]),
        )

    return hcat(DataFrame(time = d_out), DataFrame(simulated_fluxes))
end
function plot_simulated_fluxes_vs_reference(simulated_fluxes_arg, reference_arg; labels = ["current code" "reference simulation"], kwargs...)
    d_out = simulated_fluxes_arg.time
    simulated_fluxes = NamedTuple([k=>v for (k,v) in pairs(eachcol(simulated_fluxes_arg))])
    reference = Dict([String(k)=>v for (k,v) in pairs(eachcol(reference_arg))])

    evap_key = "evap" ∈ keys(reference) ? "evap" : "cum_d_evap"
    Plots.plot(
        Plots.plot(d_out, [simulated_fluxes.cum_d_prec      reference["cum_d_prec"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_prec", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_rfal      reference["cum_d_rfal"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_rfal", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_sfal      reference["cum_d_sfal"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_sfal", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_rint      reference["cum_d_rint"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_rint", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_sint      reference["cum_d_sint"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_sint", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_rthr      reference["cum_d_rthr"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_rthr", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_sthr      reference["cum_d_sthr"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_sthr", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_rsno      reference["cum_d_rsno"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_rsno", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_rnet      reference["cum_d_rnet"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_rnet", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_smlt      reference["cum_d_smlt"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_smlt", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_tran      reference["cum_d_tran"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_tran", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_irvp      reference["cum_d_irvp"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_irvp", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_isvp      reference["cum_d_isvp"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_isvp", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_slvp      reference["cum_d_slvp"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_slvp", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_snvp      reference["cum_d_snvp"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_snvp", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_pint      reference["cum_d_pint"]],     linestyle = [:solid :dot], label = labels, title = "accum.cum_d_pint", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.cum_d_ptran     reference["cum_d_ptran"]],    linestyle = [:solid :dot], label = labels, title = "accum.cum_d_ptran", kwargs...),
        # Plots.plot(d_out, [simulated_fluxes.cum_d_pslvp     reference["cum_d_pslvp"]],    linestyle = [:solid :dot], label = labels, title = "accum.cum_d_pslvp", kwargs...), # PSLVP was removed from output,
        # Plots.plot(d_out, [simulated_fluxes.evap            reference["evap"]],           linestyle = [:solid :dot], label = labels, title = "accum.evap => fate", kwargs...),
        # Plots.plot(d_out, [simulated_fluxes.evap            reference["cum_d_evap"]],     linestyle = [:solid :dot], label = labels, title = "accum.evap => fate", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.evap            reference[evap_key]],         linestyle = [:solid :dot], label = labels, title = "accum.evap => fate", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.flow            reference["flow"]],           linestyle = [:solid :dot], label = labels, title = "accum.flow => fate", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.seep            reference["seep"]],           linestyle = [:solid :dot], label = labels, title = "accum.seep => fate", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.srfl            reference["srfl"]],           linestyle = [:solid :dot], label = labels, title = "accum.srfl", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.slfl            reference["slfl"]],           linestyle = [:solid :dot], label = labels, title = "accum.slfl", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.byfl            reference["byfl"]],           linestyle = [:solid :dot], label = labels, title = "accum.byfl", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.dsfl            reference["dsfl"]],           linestyle = [:solid :dot], label = labels, title = "accum.dsfl", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.gwfl            reference["gwfl"]],           linestyle = [:solid :dot], label = labels, title = "accum.gwfl", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.vrfln           reference["vrfln"]],          linestyle = [:solid :dot], label = labels, title = "accum.vrfln", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.StorageSWAT     reference["StorageSWAT"]],    linestyle = [:solid :dot], label = labels, title = "accum.StorageSWAT", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.StorageWATER    reference["StorageWATER"]],   linestyle = [:solid :dot], label = labels, title = "accum.StorageWATER", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.BALERD_SWAT     reference["BALERD_SWAT"]],    linestyle = [:solid :dot], label = labels, title = "accum.BALERD_SWAT", kwargs...),
        Plots.plot(d_out, [simulated_fluxes.BALERD_total    reference["BALERD_total"]],   linestyle = [:solid :dot], label = labels, title = "accum.BALERD_total", kwargs...)
    )
end
function test_fluxes_comparison(simulated_fluxes_arg, reference_arg)
    simulated_fluxes = NamedTuple([k=>v for (k,v) in pairs(eachcol(simulated_fluxes_arg))])
    reference = Dict([String(k)=>v for (k,v) in pairs(eachcol(reference_arg))])

    @testset "test_fluxes_comparison" begin
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
        @test isapprox(reference["evap"],      simulated_fluxes.evap,     atol = 1e-4, rtol = 1e-4) # TODO: somehow this does not work on CI?
        @test isapprox(reference["cum_d_tran"],     simulated_fluxes.cum_d_tran,     atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
        @test isapprox(reference["cum_d_irvp"],     simulated_fluxes.cum_d_irvp,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_isvp"],     simulated_fluxes.cum_d_isvp,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_slvp"],     simulated_fluxes.cum_d_slvp,     atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
        @test isapprox(reference["cum_d_snvp"],     simulated_fluxes.cum_d_snvp,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_pint"],     simulated_fluxes.cum_d_pint,     atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["cum_d_ptran"],    simulated_fluxes.cum_d_ptran,    atol = 1e-4, rtol = 1e-4)
        # @test isapprox(reference["cum_d_pslvp"],  simulated_fluxes.cum_d_pslvp,    atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["flow"],           simulated_fluxes.flow,           atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
        @test isapprox(reference["seep"],           simulated_fluxes.seep,           atol = 1e-4, rtol = 1e-4)
        @test isapprox(reference["srfl"],           simulated_fluxes.srfl,           atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
        @test isapprox(reference["slfl"],           simulated_fluxes.slfl,           atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
        @test isapprox(reference["byfl"],           simulated_fluxes.byfl,           atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
        @test isapprox(reference["dsfl"],           simulated_fluxes.dsfl,           atol = 1e-4, rtol = 1e-4) #
        @test isapprox(reference["gwfl"],           simulated_fluxes.gwfl,           atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
        @test isapprox(reference["vrfln"],          simulated_fluxes.vrfln,          atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
        @test isapprox(reference["StorageSWAT"],    simulated_fluxes.StorageSWAT,    atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
        @test isapprox(reference["StorageWATER"],   simulated_fluxes.StorageWATER,   atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
        @test isapprox(reference["BALERD_SWAT"],    simulated_fluxes.BALERD_SWAT,    atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
        @test isapprox(reference["BALERD_total"],   simulated_fluxes.BALERD_total,   atol = 1e-1, rtol = 1e-1) # TODO: somehow this does not work on CI?
    end
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

    df_simulatedFluxes2 = get_daily_soilFluxes(simulation2);
    df_simulatedFluxes2b = get_daily_soilFluxes(simulation2_withVariousFlows);

    df_simulatedFluxes2c = get_daily_soilFluxes_new(simulation2);
    df_simulatedFluxes2d = get_daily_soilFluxes_new(simulation2_withVariousFlows);

    fname2  = "../examples/INFEXP1_fluxes_referencedf.csv"
    fname2b = "../examples/INFEXP1-modified_fluxes_referencedf.csv"
    loadeddf2  = read(fname2, DataFrame)
    loadeddf2b = read(fname2b, DataFrame)
    if task == "test"
        test_fluxes_comparison(df_simulatedFluxes2,  loadeddf2)
        test_fluxes_comparison(df_simulatedFluxes2b, loadeddf2b)
        #test_fluxes_comparison(df_simulatedFluxes2c, loadeddf2) # new format
        #test_fluxes_comparison(df_simulatedFluxes2d, loadeddf2b)# new format
        test_fluxes_comparison(df_simulatedFluxes2, df_simulatedFluxes2c) # check that new format matches old format
        test_fluxes_comparison(df_simulatedFluxes2b, df_simulatedFluxes2d) # check that new format matches old format
    elseif task == "overwrite" && !is_a_CI_system # only overwrite on local machine, never on CI
        # overwrite output
        write(fname2,  df_simulatedFluxes2)
        write(fname2b, df_simulatedFluxes2b)
        error("Test overwrites reference solution instead of checking against it.")
    else
        # do nothing
    end

    if !is_a_CI_system && plot_flag
        fname_illustrations = "out/$git_status_string/TESTSET_Oversaturation-infiltration-FLUXES"
        mkpath(dirname(fname_illustrations))

        pl_fluxes2 = plot_simulated_fluxes_vs_reference(df_simulatedFluxes2, loadeddf2);
        Plots.plot!(legend = :topleft, size=(2000,1000), layout = (7,5))
        savefig(Plots.plot(pl_fluxes2, size=(2000,1000), dpi=300),  fname_illustrations*"_fluxes_noBYFL_regressionTest.png")

        pl_fluxes2b = plot_simulated_fluxes_vs_reference(df_simulatedFluxes2b, loadeddf2b);
        Plots.plot!(legend = :topleft, size=(2000,1000), layout = (7,5))
        savefig(Plots.plot(pl_fluxes2b, size=(2000,1000), dpi=300),  fname_illustrations*"_fluxes_withBYFL_regressionTest.png")

        pl_fluxes2vs2b = plot_simulated_fluxes_vs_reference(
            df_simulatedFluxes2,df_simulatedFluxes2b,labels = ["Without BYFL" "With BYFL"]);
        Plots.plot!(legend = :topleft, size=(2000,1000), layout = (7,5))
        savefig(Plots.plot(pl_fluxes2vs2b, size=(2000,1000), dpi=300),  fname_illustrations*"_fluxes_BYFL_comparison.png")
        pl_comparison = Plots.plot(plotamounts(simulation2, title = "Without BYFL"),
            plotamounts(simulation2_withVariousFlows, title = "With BYFL"),
            layout = (1,2), size=(1460,1400))
        savefig(Plots.plot(pl_comparison, size=(2000,1000), dpi=300),  fname_illustrations*"_fluxes_BYFL_comparison2.png")
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
    df_simulatedFluxes = get_daily_soilFluxes(example_result3);
    df_simulatedFluxes_new = get_daily_soilFluxes_new(example_result3);

    # test or overwrite
    fname = "../examples/DAV2020-full-modified_fluxes_referencedf.csv"
    loadeddf = read(fname, DataFrame)
    if task == "test"
        test_fluxes_comparison(df_simulatedFluxes, loadeddf)
        test_fluxes_comparison(df_simulatedFluxes_new, loadeddf)
        test_fluxes_comparison(df_simulatedFluxes_new, df_simulatedFluxes) # check that new format matches old format
    elseif task == "overwrite" && !is_a_CI_system # only overwrite on local machine, never on CI
        # overwrite output
        write(fname, df_simulatedFluxes)
        error("Test overwrites reference solution instead of checking against it.")
    else
        # do nothing
    end



    if !is_a_CI_system && plot_flag
        fname_illustrations = "out/$git_status_string/TESTSET_DAV2020modified-regressionTest-FLUXES"
        mkpath(dirname(fname_illustrations))

        pl_fluxes = plot_simulated_fluxes_vs_reference(df_simulatedFluxes, loadeddf);
        Plots.plot!(legend = :topright, size=(2000,1000), layout = (7,5))
        savefig(Plots.plot(pl_fluxes, size=(2000,1400), layout = (7,5), dpi=300),  fname_illustrations*"_fluxes.png")
    end

end