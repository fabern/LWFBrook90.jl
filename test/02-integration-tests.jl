# Integration tests

# - _Unit testing_ asserts that individual pieces of a project work as expected. (developers
#       perspective)
# - _Integration testing_ asserts that they fit together as expected. Also known as
#       _functional tests_, they cover entire use cases (user perspective). For LWFBrook90.jl
#       these are tests that are compared to e.g. LWFBrook90R or Hydrus.
# - _Regression testing_ asserts that behavior is unchanged over time. Also known as
#       _reference tests_.

include("fct-helpers-for-integration-tests.jl")

# NOTE: locally, i.e. not on CI system, one might need to do manually cd("test")
if basename(pwd()) != "test"; cd("test"); end


# @testset "BEA-2016-θ-ψ-aboveground-states" begin
# @testset "DAV-2020-θ-ψ-aboveground-states" begin
# # source: https://stackoverflow.com/a/63871951
@testset "trial1" begin
    soil_horizons = LWFBrook90.init_soil("test-assets/DAV-2020/input-files/DAV_LW1_def_soil_horizons.csv")
    @test soil_horizons isa DataFrame
end
@testset "trial2" begin
    res = loadSPAC("test-assets/DAV-2020/input-files/", "DAV_LW1_def"; simulate_isotopes = false)
    @test res isa SPAC
end
@testset "loadSPAC()-parametrization of root distribution and IC" begin
    # input_path = "test-assets/Hammel-2001/input-files-ISO"
    # input_prefix = "Hammel_loam-NLayer-27-RESET=FALSE"
    input_prefix = "isoBEAdense2010-18-reset-FALSE";
    input_path = "../examples/isoBEAdense2010-18-reset-FALSE-input/";
    simulate_isotopes = true

    model = loadSPAC(input_path, input_prefix; simulate_isotopes = simulate_isotopes);
    @test model.pars.root_distribution == "soil_discretization.csv"
    @test model.pars.IC_soil == "soil_discretization.csv"
    @test length(model.soil_discretization.Δz) == 21
    model = loadSPAC(input_path, input_prefix; simulate_isotopes = simulate_isotopes,
                    root_distribution = (beta = 0.90, z_rootMax_m = nothing));
    @test model.pars.root_distribution == (beta = 0.90, z_rootMax_m = nothing)
    @test model.pars.IC_soil == "soil_discretization.csv"
    model = loadSPAC(input_path, input_prefix; simulate_isotopes = simulate_isotopes,
                    IC_soil           = (PSIM_init_kPa = -12., delta18O_init_permil = -22., delta2H_init_permil = -88.));
    @test all(model.soil_discretization.df.uAux_PSIM_init_kPa .== -12.)
    model = loadSPAC(input_path, input_prefix; simulate_isotopes = simulate_isotopes,
                    root_distribution = (beta = 0.90, z_rootMax_m = nothing),
                    IC_soil           = (PSIM_init_kPa = -12., delta18O_init_permil = -22., delta2H_init_permil = -88.));
    model
    model = loadSPAC(input_path, input_prefix; simulate_isotopes = simulate_isotopes,
                    Δz_thickness_m = [0.5, 0.5],
                    root_distribution = (beta = 0.90, z_rootMax_m = nothing),
                    IC_soil           = (PSIM_init_kPa = -12., delta18O_init_permil = -22., delta2H_init_permil = -88.));
    @test nrow(model.soil_discretization.df) == 2
    @test model.soil_discretization.df.Rootden_ ≈ [0.019896924495853598, 0.00010254427616865014]
    @test model.soil_discretization.df.uAux_PSIM_init_kPa == [-12., -12.]
    @test model.soil_discretization.df.u_delta18O_init_permil == [-22., -22.]
    @test model.soil_discretization.df.u_delta2H_init_permil  == [-88., -88.]
    @test_throws ErrorException loadSPAC(input_path, input_prefix; simulate_isotopes = simulate_isotopes,
                    Δz_thickness_m = [0.5, 0.5]);
    @test_throws ErrorException loadSPAC(input_path, input_prefix; simulate_isotopes = simulate_isotopes,
                    Δz_thickness_m = [0.5, 0.5],
                    IC_soil           = (PSIM_init_kPa = -12., delta18O_init_permil = -22., delta2H_init_permil = -88.));
    @test_throws ErrorException loadSPAC(input_path, input_prefix; simulate_isotopes = simulate_isotopes,
                    Δz_thickness_m = [0.5, 0.5],
                    root_distribution = (beta = 0.90, z_rootMax_m = nothing));
end

@testset "loadSPAC()-manual_soil_discretization" begin
    input_prefix = "isoBEAdense2010-18-reset-FALSE";
    input_path = "../examples/isoBEAdense2010-18-reset-FALSE-input/";
    simulate_isotopes = true
    Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 5); 0.1]; # grid spacing (heterogenous), meter (N=21)

    # TODO: make this work without parametric root distribution (use)
    @test_throws "no parametric root_distribution provided" loadSPAC(input_path, input_prefix; Δz_thickness_m = Δz_m)
    @test_throws "no parametric soil initial conditions" loadSPAC(input_path, input_prefix; Δz_thickness_m = Δz_m, root_distribution = (beta = 0.98))
    # TODO: once it works without parametric root and IC, below two tests can be simplified by removign the two arguments

    # Check extension
    @test loadSPAC(input_path, "isoBEAdense2010-18-reset-FALSE";
        root_distribution = (beta = 0.98, ), IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -9.0, delta2H_init_permil = -11.0),
        # test is about these added 10cm below soil_horizon:
        Δz_thickness_m = [Δz_m; 0.1]).soil_discretization.df.Lower_m[end] ≈ -1.3
    @test loadSPAC(input_path, "isoBEAdense2010-18-reset-FALSE";
        root_distribution = (beta = 0.98, ), IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -9.0, delta2H_init_permil = -11.0),
        # test is about these added 10cm below soil_horizon:
        Δz_thickness_m = [Δz_m; 100]).soil_discretization.df.Lower_m[end] ≈ -101.2

    @test_logs (:warn, r"soil horizon will be extended") match_mode=:any loadSPAC(input_path, "isoBEAdense2010-18-reset-FALSE";
        root_distribution = (beta = 0.98, ), IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -9.0, delta2H_init_permil = -11.0),
        # test is about these added 10cm below soil_horizon:
        Δz_thickness_m = [Δz_m; 0.7]);

    model = loadSPAC(
        input_path, input_prefix;
        simulate_isotopes = simulate_isotopes,
        root_distribution = (beta = 0.90, z_rootMax_m = nothing));

    @test_throws AssertionError simulation = setup(model; soil_output_depths_m = [0.02, 0.42])
    @test_throws "Requested soil_output_depths_m (additional layers)" simulation = setup(model; soil_output_depths_m = [-0.03, -0.11, -0.112])
end

@testset "$site-θ-ψ-aboveground-states" for site in ["BEA-2016" "DAV-2020"]
    out_figure_string = "test-assets/$site/out"
    folder_with_sim_input_and_ref_output = "test-assets/$site"
    input_prefix = ifelse(site == "BEA-2016", "BEA2016-reset-FALSE", "DAV_LW1_def")
    NLAYERBASE = ifelse(site == "BEA-2016", 7, 5)

    # # Check the RMSE of θ in simulations is below a limit
    sim, ref_NLAYER7, ref_NLAYER14, ref_NLAYER21, ref_NLAYER70, depth_to_read_out_mm =
        prepare_θψAboveground_from_sim_and_ref(folder_with_sim_input_and_ref_output,input_prefix; NLAYERBASE=NLAYERBASE);
    # # Use sensible accuracy values to compare the two solutions (e.g. θ of 0.02, and ψ of 1 kPa)
    # # (e.g. RMSE(reference, simulated) < [hardcoded_value])
    # # Floating point accuracy is used by regression tests.
    # # (Or could be used, once I compare with a very high-fidelity solution, of which I known
    # #  that my code should converge to. E.g. a high-resolution Hydrus simulation.)

    # Compare with LWFBrook90R as reference solution
    # θ:
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"θ_")], ref_NLAYER7.θ) < 0.02
    # ψ:
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"ψ_")], ref_NLAYER7.ψ) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.24,5)
    # "GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)", (not done for: "CC (MJ/m2)" "SNOWLQ (mm)"]):
    @test RMS_differences(sim.above[:,[:time,:GWAT, :INTS, :INTR, :SNOW]],
                            ref_NLAYER7.above[Not(end),[:time, :gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)

    # Note that below we compare a NLAYER7 LWFBrook90.jl solution, with finer resolved
    # LWFBrook90R solutions. It is therefore normal, that the uncertainty can increase...
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"θ_")], ref_NLAYER14.θ) < 0.03
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"ψ_")], ref_NLAYER14.ψ) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.65,5.0)
    @test RMS_differences(sim.above[:,[:time,:GWAT,:INTS,:INTR,:SNOW]],
                            ref_NLAYER14.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"θ_")], ref_NLAYER21.θ) < 0.04
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"ψ_")], ref_NLAYER21.ψ) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.6,5.0)
    @test RMS_differences(sim.above[:,[:time,:GWAT,:INTS,:INTR,:SNOW]],
                            ref_NLAYER21.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"θ_")], ref_NLAYER70.θ) < 0.04
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"ψ_")], ref_NLAYER70.ψ) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.7,5.0)
    @test RMS_differences(sim.above[:,[:time,:GWAT,:INTS,:INTR,:SNOW]],
                            ref_NLAYER70.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)

    # TODO(bernhard): we could run multiple LWFBrook90.jl simulations and compare with the
    # finest LWFBrook90R simulation only.

    if !is_a_CI_system && plot_flag
        # if some error appears, the following code can be used to plot the solutions
        using Plots, Measures
        # fname_illustrations = "out/$(today())/"
        fname_illustrations = "out/$git_status_string/"
        mkpath(dirname(fname_illustrations))

        pl_θ = plot(sim.θψ.time,
                Matrix(sim.θψδ[:,[:θ_100mm, :θ_500mm, :θ_1000mm]]),
                line = :solid, labels = "LWFBrook90.jl:" .* string.(depth_to_read_out_mm) .* "mm",
                ylabel = "θ (-)", legend_position = :bottomright)
        plot!(Matrix(ref_NLAYER7.θ[:,Not(:time)]), line = :dash, color = :black,
                labels = "LWFBrook90R_NLayer7:" .* string.(depth_to_read_out_mm) .* "mm")
        # plot!(Matrix(ref_NLAYER14.θ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer14")
        # plot!(Matrix(ref_NLAYER21.θ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer21")
        # plot!(Matrix(ref_NLAYER70.θ[:,Not(:time)]), line = :dash, color = :black, labels = "LWFBrook90R_NLayer70")
        pl_ψ = plot(sim.θψ.time,
                Matrix(sim.θψδ[:,[:ψ_100mm, :ψ_500mm, :ψ_1000mm]]),
                line = :solid, labels = "LWFBrook90.jl:" .* string.(depth_to_read_out_mm) .* "mm",
                ylabel = "ψ (kPa)", legend_position = :bottomright)
        plot!(Matrix(ref_NLAYER7.ψ[:,Not(:time)]), line = :dash, color = :black,
                labels = "LWFBrook90R_NLayer7:" .* string.(depth_to_read_out_mm) .* "mm")
        # plot!(Matrix(ref_NLAYER14.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer14")
        # plot!(Matrix(ref_NLAYER21.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer21")
        # plot!(Matrix(ref_NLAYER70.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer70")
        pl_a = plot(sim.above.time,
                Matrix(#sim.above[:,Not(:time)]),
                       # label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"],
                        sim.above[:,[:GWAT,:INTS,:INTR,:SNOW]]),
                        label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"],
                        line = :solid)
        plot!(Matrix(ref_NLAYER7.above[:,[:intr,:ints,:snow,:gwat]]),
                line = :dash, color = :black, labels = "LWFBrook90R_NLayer7")
        # plot!(Matrix(ref_NLAYER14.above[:,[:intr,:ints,:snow,:gwat]]),
        #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
        # plot!(Matrix(ref_NLAYER21.above[:,[:intr,:ints,:snow,:gwat]]),
        #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
        # plot!(Matrix(ref_NLAYER70.above[:,[:intr,:ints,:snow,:gwat]]),
        #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
        plot(pl_θ, pl_ψ, pl_a, layout = (3,1), size = (600,800),
            leftmargin = 5mm)
        savefig(fname_illustrations*"TESTSET_site-θ-ψ-aboveground-states_$(site).png")
    end
end




# prepare plotting functions for Hammel (with/without δ)
function my_scatter!(pl, df; args...)
    scatter!(pl, float.(df[:,:time]), Matrix(df[:,Not(:time)]); args...)
end
function plot_Hammel_Dense(sim, ref, hyd, depth_to_read_out_mm, title; subtitle = "", args...)
    #LWFBrook90.jl:
    pl_θ = plot(sim.θψδ_dense.time, Matrix(sim.θψδ_dense[:,r"θ_"]); line = :solid, color = [1 2 3 4 5], labels = "LWFBrook90.jl: " .* permutedims(names(sim.θψδ_dense[:,r"θ_"])))
    pl_ψ = plot(sim.θψδ_dense.time, Matrix(sim.θψδ_dense[:,r"ψ_"]); line = :solid, color = [1 2 3 4 5], labels = "LWFBrook90.jl: " .* permutedims(names(sim.θψδ_dense[:,r"ψ_"])))
    #LWFBrook90R
    my_scatter!(pl_θ, ref.θ;       line = :dot,   color = :black,      labels = ["LWFBrook90R" "" "" "" ""], markersize=2.5)
    my_scatter!(pl_ψ, ref.ψ;       line = :dot,   color = :black,      labels = ["LWFBrook90R" "" "" "" ""], markersize=2.5)
    #Hydrus
    my_scatter!(pl_θ,hyd.θdense, line = :dash,  color = [1 2 3 4 5], markersize=2, markerstrokecolor=:auto, labels = "Hydrus: " .* string.(permutedims(depth_to_read_out_mm)) .* " mm")
    my_scatter!(pl_ψ,hyd.ψdense, line = :dash,  color = [1 2 3 4 5], markersize=2, markerstrokecolor=:auto, labels = "Hydrus: " .* string.(permutedims(depth_to_read_out_mm)) .* " mm")
    # my_plot!(   pl_θ,hyd.θdense, line = :dash,  color = [1 2 3 4 5],               labels = nothing)
    # my_plot!(   pl_ψ,hyd.ψdense, line = :dash,  color = [1 2 3 4 5],               labels = nothing)

    simulate_isotopes = !ismissing(sim.θψδ_dense.δ18O_100mm[1])
    if (simulate_isotopes)
        #LWFBrook90.jl:
        pl_δ18O = plot(sim.θψδ_dense.time, Matrix(sim.θψδ_dense[:,r"δ18O_"]); line = :solid, color = [1 2 3 4 5], labels = "LWFBrook90.jl: " .* permutedims(names(sim.θψδ_dense[:,r"δ18O_"])))
        pl_δ2H  = plot(sim.θψδ_dense.time, Matrix(sim.θψδ_dense[:,r"δ18O_"]);  line = :solid, color = [1 2 3 4 5], labels = "LWFBrook90.jl: " .* permutedims(names(sim.θψδ_dense[:,r"δ18O_"])))
        #Hydrus
        my_scatter!(pl_δ18O,hyd.δ18Odense, line = :dash,  color = [1 2 3 4 5], markersize=2, markerstrokecolor=:auto, labels = "Hydrus: " .* string.(permutedims(depth_to_read_out_mm)) .* " mm")
        my_scatter!(pl_δ2H,hyd.δ2Hdense, line = :dash,  color = [1 2 3 4 5], markersize=2, markerstrokecolor=:auto, labels = "Hydrus: " .* string.(permutedims(depth_to_read_out_mm)) .* " mm")
    else
        pl_δ18O = plot()
        pl_δ2H = plot()
    end

    ODE_method = replace(String(Symbol(sim.ODEsolution.alg)), r"(.*?)(\(|{).*"=>s"\1")
    plot(plot_title = Printf.@sprintf("%s \n %d time steps with method: %s %s",
            title,
            sim.ODEsolution.destats.naccept,
            ODE_method, ifelse(subtitle=="","","\n"*subtitle)),
        filter(!isnothing,
                [plot!(pl_θ,ylabel = "θ (-)", xlabel = "Time (days)"),
                plot!(pl_ψ,ylabel = "ψ (kPa)", xlabel = "Time (days)", legend = false),
                ifelse(simulate_isotopes, plot!(pl_δ18O,ylabel = "δ18O (‰)", xlabel = "Time (days)", legend = false), nothing),
                ifelse(simulate_isotopes, plot!(pl_δ2H, ylabel = "δ2H (‰)", xlabel = "Time (days)", legend = false),  nothing)])... ,
        size = (1200,1200), layout = ifelse(simulate_isotopes, (4,1), (2,1)), leftmargin = 8mm,
        args...)
end
print_timed_statistics = function(time, gctime, gcstats)
    # Printf.@sprintf("%.2f seconds (%.2f M allocations: %.3f MiB, %.2f%% gc time)",
    #         time, gcstats.poolalloc/10^6, gcstats.allocd/1024^2, gctime*100)
    Printf.@sprintf("%.2f seconds (%.2f M allocations: %.3f MiB)",
            time, gcstats.poolalloc/10^6, gcstats.allocd/1024^2)
end

@testset "Hammel-2001-θ-ψ" begin
    # This test set checks that the RMSE of θ and ψ vs reference solutions are below hard-
    # coded limits. Two simulations defined in Hammel et al. (2001) are run.

    @show pwd() # This is to help get the folder right.

    # 1) Run simulation and load references (LWFBrook90R and Hydrus1D)
    (sim4, ref4, hyd4), time4, bytes4, gctime4, gcstats4 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_sand-NLayer-27-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand");
    (sim5, ref5, hyd5), time5, bytes5, gctime5, gcstats5 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_sand-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand");
    high_resolution_flag && begin (sim6, ref6, hyd6), time6, bytes6, gctime6, gcstats6 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_sand-NLayer-400-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand") end;
    (sim1, ref1, hyd1), time1, bytes1, gctime1, gcstats1 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_loam-NLayer-27-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam");
    (sim2, ref2, hyd2), time2, bytes2, gctime2, gcstats2 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_loam-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam");
    high_resolution_flag && begin (sim3, ref3, hyd3), time3, bytes3, gctime3, gcstats3 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_loam-NLayer-400-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam") end;

    # 2) Plot (optional, not done when testing in CI)
    # Illustrate with a plot what will be compared in the tests below
    if !is_a_CI_system && plot_flag
        # if (true) # Do these manually outside of automatic testing in order not to require Plots pkg
        depth_to_read_out_mm = [100, 500, 1000, 1500, 1900]
        # fname_illustrations = "out/$(today())/"

        fname_illustrations = "out/$git_status_string/"
        mkpath(dirname(fname_illustrations))

        # Plot Loam simulations
        pl1 = plot_Hammel_Dense(sim1, ref1, hyd1, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(700,900), dpi=300,
            subtitle = print_timed_statistics(time1, gctime1, gcstats1)*"\n"*git_status_string, topmargin = 15mm)
        pl2 = plot_Hammel_Dense(sim2, ref2, hyd2, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(700,900), dpi=300,
            subtitle = print_timed_statistics(time2, gctime2, gcstats2)*"\n"*git_status_string, topmargin = 15mm)

        # Plot Sand simulations
        pl4 = plot_Hammel_Dense(sim4, ref4, hyd4, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(700,900), dpi=300,
            subtitle = print_timed_statistics(time4, gctime4, gcstats4)*"\n"*git_status_string, topmargin = 15mm)
        pl5 = plot_Hammel_Dense(sim5, ref5, hyd5, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(700,900), dpi=300,
            subtitle = print_timed_statistics(time5, gctime5, gcstats5)*"\n"*git_status_string, topmargin = 15mm)

        high_resolution_flag && (pl3 = plot_Hammel_Dense(sim3, ref3, hyd3, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(700,900), dpi=300,
            subtitle = print_timed_statistics(time3, gctime3, gcstats3)*"\n"*git_status_string, topmargin = 15mm))
        high_resolution_flag && (pl6 = plot_Hammel_Dense(sim6, ref6, hyd6, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(700,900), dpi=300,
            subtitle = print_timed_statistics(time6, gctime6, gcstats6)*"\n"*git_status_string, topmargin = 15mm))

        savefig(pl1, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Loam_sim1.png")
        savefig(pl2, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Loam_sim2.png")
        high_resolution_flag && (savefig(pl3, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Loam_sim3.png"))
        savefig(pl4, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Sand_sim1.png")
        savefig(pl5, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Sand_sim2.png")
        high_resolution_flag && (savefig(pl6, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Sand_sim3.png"))
    end

    # Use sensible accuracy values to compare the two solutions (e.g. θ of 0.02, and ψ of 1 kPa)
    # (e.g. RMSE(reference, simulated) < [hardcoded_value])
    # Floating point accuracy is used by regression tests.
    # (Or could be used, once I compare with a very high-fidelity solution, of which I known
    #  that my code should converge to. E.g. a high-resolution Hydrus simulation.)

    # 3) Test RMSE
    # Compare θ
    # Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.θψδ[:,Cols(:time, r"θ_")], ref1.θ) < 0.0020
    @test RMS_differences(sim2.θψδ[:,Cols(:time, r"θ_")], ref2.θ) < 0.00070
    high_resolution_flag && (@test RMS_differences(sim3.θψδ[:,Cols(:time, r"θ_")], ref3.θ) < 0.00053)
    @test RMS_differences(sim4.θψδ[:,Cols(:time, r"θ_")], ref4.θ) < 0.00035
    @test RMS_differences(sim5.θψδ[:,Cols(:time, r"θ_")], ref5.θ) < 0.00035
    high_resolution_flag && (@test RMS_differences(sim6.θψδ[:,Cols(:time, r"θ_")], ref6.θ) < 0.00035)

    # Compare with Hydrus1D
    @test RMS_differences(sim1.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd1.θ[Not(1),:]) < 0.005
    @test RMS_differences(sim2.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd2.θ[Not(1),:]) < 0.002
    high_resolution_flag && (@test RMS_differences(sim3.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd3.θ[Not(1),:]) < 0.0015)
    @test RMS_differences(sim4.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd4.θ[Not(1),:]) < 0.008
    @test RMS_differences(sim5.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd5.θ[Not(1),:]) < 0.009
    high_resolution_flag && (@test RMS_differences(sim6.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd6.θ[Not(1),:]) < 0.007)

    # Compare ψ
    # Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.θψδ[:,Cols(:time, r"ψ_")], ref1.ψ) < 2.0 # kPa
    @test RMS_differences(sim2.θψδ[:,Cols(:time, r"ψ_")], ref2.ψ) < 0.60 # kPa
    high_resolution_flag && (@test RMS_differences(sim3.θψδ[:,Cols(:time, r"ψ_")], ref3.ψ) < 0.42) # kPa
    @test RMS_differences(sim4.θψδ[:,Cols(:time, r"ψ_")], ref4.ψ) < 0.004 # kPa
    @test RMS_differences(sim5.θψδ[:,Cols(:time, r"ψ_")], ref5.ψ) < 0.0080 # kPa
    high_resolution_flag && (@test RMS_differences(sim6.θψδ[:,Cols(:time, r"ψ_")], ref6.ψ) < 0.0025) # kPa

    # Compare with Hydrus1D
    @test RMS_differences(sim1.θψδ[Not([end-1, end]),Cols(:time, r"ψ_")], hyd1.ψ[Not(1),:]) < 6.1 # kPa
    @test RMS_differences(sim2.θψδ[Not([end-1, end]),Cols(:time, r"ψ_")], hyd2.ψ[Not(1),:]) < 2.8 # kPa
    high_resolution_flag && (@test RMS_differences(sim3.θψδ[Not([end-1, end]),Cols(:time, r"ψ_")], hyd3.ψ[Not(1),:]) < 2.1) # kPa
    @test RMS_differences(sim4.θψδ[Not([end-1, end]),Cols(:time, r"ψ_")], hyd4.ψ[Not(1),:]) < 1.0 # kPa
    @test RMS_differences(sim5.θψδ[Not([end-1, end]),Cols(:time, r"ψ_")], hyd5.ψ[Not(1),:]) < 1.5 # kPa
    high_resolution_flag && (@test RMS_differences(sim6.θψδ[Not([end-1, end]),Cols(:time, r"ψ_")], hyd6.ψ[Not(1),:]) < 1.0) # kPa
end

# # NOTE: locally, i.e. not on CI system, one might need to do manually cd("test")
@testset "Hammel-2001-θ-ψ-δ" begin
    # This test set checks that the RMSE of θ and ψ vs reference solutions are below hard-
    # coded limits. Two simulations defined in Hammel et al. (2001) are run.

    @show pwd() # This is to help get the folder right.

    # # 1) Run simulation and load references (LWFBrook90R and Hydrus1D)
    (sim4, ref4, hyd4), time4, bytes4, gctime4, gcstats4 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_sand-NLayer-27-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand_ISO2",
        simulate_isotopes   = true);
    (sim5, ref5, hyd5), time5, bytes5, gctime5, gcstats5 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_sand-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand_ISO2",
        simulate_isotopes   = true);
    high_resolution_flag && begin (sim6, ref6, hyd6), time6, bytes6, gctime6, gcstats6 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_sand-NLayer-400-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand_ISO2",
        simulate_isotopes   = true) end;

    (sim1, ref1, hyd1), time1, bytes1, gctime1, gcstats1 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_loam-NLayer-27-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam_ISO2",
        simulate_isotopes   = true);
    (sim2, ref2, hyd2), time2, bytes2, gctime2, gcstats2 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_loam-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam_ISO2",
        simulate_isotopes   = true);
    high_resolution_flag && begin (sim3, ref3, hyd3), time3, bytes3, gctime3, gcstats3 = @timed prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_loam-NLayer-400-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam_ISO2",
        simulate_isotopes   = true) end;

    # 2) Plot (optional, not done when testing in CI)
    # Illustrate with a plot what will be compared in the tests below
    if !is_a_CI_system && plot_flag
        # if (true) # Do these manually outside of automatic testing in order not to require Plots pkg
        using Plots, Measures
        depth_to_read_out_mm = [100, 500, 1000, 1500, 1900]
        # fname_illustrations = "out/$(today())/"

        fname_illustrations = "out/$git_status_string/"
        mkpath(dirname(fname_illustrations))

        # Plot Loam simulations
        pl1 = plot_Hammel_Dense(sim1, ref1, hyd1, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time1, gctime1, gcstats1)*"\n"*git_status_string,
            topmargin = 17mm, layout = (2,2))
        pl2 = plot_Hammel_Dense(sim2, ref2, hyd2, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time2, gctime2, gcstats2)*"\n"*git_status_string,
            topmargin = 17mm, layout = (2,2))
        high_resolution_flag && (pl3 = plot_Hammel_Dense(sim3, ref3, hyd3, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time3, gctime3, gcstats3)*"\n"*git_status_string,
            topmargin = 17mm, layout = (2,2)))
        # Plot Sand simulations
        pl4 = plot_Hammel_Dense(sim4, ref4, hyd4, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time4, gctime4, gcstats4)*"\n"*git_status_string,
            topmargin = 17mm, layout = (2,2))
        pl5 = plot_Hammel_Dense(sim5, ref5, hyd5, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time5, gctime5, gcstats5)*"\n"*git_status_string,
            topmargin = 17mm, layout = (2,2))
        high_resolution_flag && (pl6 = plot_Hammel_Dense(sim6, ref6, hyd6, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time6, gctime6, gcstats6)*"\n"*git_status_string,
            topmargin = 17mm, layout = (2,2)))

        savefig(pl1, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Loam_sim1.png")
        savefig(pl2, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Loam_sim2.png")
        high_resolution_flag && (savefig(pl3, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Loam_sim3.png"))
        savefig(pl4, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Sand_sim1.png")
        savefig(pl5, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Sand_sim2.png")
        high_resolution_flag && (savefig(pl6, fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Sand_sim3.png"))
    end


    # Use sensible accuracy values to compare the two solutions (e.g. θ of 0.02, and ψ of 1 kPa)
    # (e.g. RMSE(reference, simulated) < [hardcoded_value])
    # Floating point accuracy is used by regression tests.
    # (Or could be used, once I compare with a very high-fidelity solution, of which I known
    #  that my code should converge to. E.g. a high-resolution Hydrus simulation.)

    # 3) Test RMSE
    # Compare θ
    # θ: Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.θψδ[:,Cols(:time, r"θ_")], ref1.θ) < 0.00194
    @test RMS_differences(sim2.θψδ[:,Cols(:time, r"θ_")], ref2.θ) < 0.00060
    high_resolution_flag && (@test RMS_differences(sim3.θψδ[:,Cols(:time, r"θ_")], ref3.θ) < 0.00053)
    @test RMS_differences(sim4.θψδ[:,Cols(:time, r"θ_")], ref4.θ) < 0.00040
    @test RMS_differences(sim5.θψδ[:,Cols(:time, r"θ_")], ref5.θ) < 0.00035
    high_resolution_flag && (@test RMS_differences(sim6.θψδ[:,Cols(:time, r"θ_")], ref6.θ) < 0.00035)

    # θ: Compare with Hydrus1D-Iso
    @test RMS_differences(sim1.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd1.θ[Not(1),:]) < 0.005
    @test RMS_differences(sim2.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd2.θ[Not(1),:]) < 0.002
    high_resolution_flag && (@test RMS_differences(sim3.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd3.θ[Not(1),:]) < 0.0013)
    @test RMS_differences(sim4.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd4.θ[Not(1),:]) < 0.007
    @test RMS_differences(sim5.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd5.θ[Not(1),:]) < 0.009
    high_resolution_flag && (@test RMS_differences(sim6.θψδ[Not([end-1, end]),Cols(:time, r"θ_")], hyd6.θ[Not(1),:]) < 0.007)

    # Compare ψ
    # ψ: Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.θψδ[:,Cols(:time, r"ψ_")], ref1.ψ) < 1.72 # kPa
    @test RMS_differences(sim2.θψδ[:,Cols(:time, r"ψ_")], ref2.ψ) < 0.50 # kPa
    high_resolution_flag && (@test RMS_differences(sim3.θψδ[:,Cols(:time, r"ψ_")], ref3.ψ) < 0.42) # kPa
    @test RMS_differences(sim4.θψδ[:,Cols(:time, r"ψ_")], ref4.ψ) < 0.003 # kPa
    @test RMS_differences(sim5.θψδ[:,Cols(:time, r"ψ_")], ref5.ψ) < 0.004 # kPa
    high_resolution_flag && (@test RMS_differences(sim6.θψδ[:,Cols(:time, r"ψ_")], ref6.ψ) < 0.0035) # kPa

    # ψ: Compare with Hydrus1D-Iso
    @test RMS_differences(sim1.θψδ[Not([end-1, end]), Cols(:time, r"ψ_")], hyd1.ψ[Not(1),:]) < 6 # kPa
    @test RMS_differences(sim2.θψδ[Not([end-1, end]), Cols(:time, r"ψ_")], hyd2.ψ[Not(1),:]) < 2.8 # kPa
    high_resolution_flag && (@test RMS_differences(sim3.θψδ[Not([end-1, end]), Cols(:time, r"ψ_")], hyd3.ψ[Not(1),:]) < 2.1) # kPa
    @test RMS_differences(sim4.θψδ[Not([end-1, end]), Cols(:time, r"ψ_")], hyd4.ψ[Not(1),:]) < 0.6 # kPa
    @test RMS_differences(sim5.θψδ[Not([end-1, end]), Cols(:time, r"ψ_")], hyd5.ψ[Not(1),:]) < 1.1 # kPa
    high_resolution_flag && (@test RMS_differences(sim6.θψδ[Not([end-1, end]), Cols(:time, r"ψ_")], hyd6.ψ[Not(1),:]) < 0.5) # kPa

    # Compare δ18O with Hydrus1D-Iso
    @test RMS_differences(sim1.θψδ[Not([end-1, end]), Cols(:time, r"δ18O_")], hyd1.δ18O[Not(1),:]) < 0.5 # unit: ‰
    @test RMS_differences(sim2.θψδ[Not([end-1, end]), Cols(:time, r"δ18O_")], hyd2.δ18O[Not(1),:]) < 0.3 # unit: ‰
    high_resolution_flag && (@test RMS_differences(sim3.θψδ[Not([end-1, end]), Cols(:time, r"δ18O_")], hyd3.δ18O[Not(1),:]) < 0.5) # unit: ‰
    @test RMS_differences(sim4.θψδ[Not([end-1, end]), Cols(:time, r"δ18O_")], hyd4.δ18O[Not(1),:]) < 0.5 # unit: ‰
    @test RMS_differences(sim5.θψδ[Not([end-1, end]), Cols(:time, r"δ18O_")], hyd5.δ18O[Not(1),:]) < 0.4 # unit: ‰
    high_resolution_flag && (@test RMS_differences(sim6.θψδ[Not([end-1, end]), Cols(:time, r"δ18O_")], hyd6.δ18O[Not(1),:]) < 0.5) # unit: ‰
    # Compare δ2H with Hydrus1D-Iso
    @test RMS_differences(sim1.θψδ[Not([end-1, end]), Cols(:time, r"δ2H_")], hyd1.δ2H[Not(1),:]) < 3    # unit: ‰
    @test RMS_differences(sim2.θψδ[Not([end-1, end]), Cols(:time, r"δ2H_")], hyd2.δ2H[Not(1),:]) < 2    # unit: ‰
    high_resolution_flag && (@test RMS_differences(sim3.θψδ[Not([end-1, end]), Cols(:time, r"δ2H_")], hyd3.δ2H[Not(1),:]) < 2.9)  # unit: ‰
    @test RMS_differences(sim4.θψδ[Not([end-1, end]), Cols(:time, r"δ2H_")], hyd4.δ2H[Not(1),:]) < 2.3  # unit: ‰
    @test RMS_differences(sim5.θψδ[Not([end-1, end]), Cols(:time, r"δ2H_")], hyd5.δ2H[Not(1),:]) < 1.65 # unit: ‰
    high_resolution_flag && (@test RMS_differences(sim6.θψδ[Not([end-1, end]), Cols(:time, r"δ2H_")], hyd6.δ2H[Not(1),:]) < 3.7)  # unit: ‰

end
