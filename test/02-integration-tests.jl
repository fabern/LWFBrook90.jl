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


# NOTE: locally, i.e. not on CI system, one might need to do manually cd("test")

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
    @test RMS_differences(sim.θ, ref_NLAYER7.θ) < 0.02
    # ψ:
    @test RMS_differences(sim.ψ, ref_NLAYER7.ψ) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.24,5)
    # "GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)", (not done for: "CC (MJ/m2)" "SNOWLQ (mm)"]):
    @test RMS_differences(sim.above[Not(1),[:time,:GWAT, :INTS, :INTR, :SNOW]],
                            ref_NLAYER7.above[Not(end),[:time, :gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)

    # Note that below we compare a NLAYER7 LWFBrook90.jl solution, with finer resolved
    # LWFBrook90R solutions. It is therefore normal, that the uncertainty can increase...
    @test RMS_differences(sim.θ, ref_NLAYER14.θ) < 0.03
    @test RMS_differences(sim.ψ, ref_NLAYER14.ψ) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.65,5.0)
    @test RMS_differences(sim.above[Not(1),[:time,:GWAT,:INTS,:INTR,:SNOW]],
                            ref_NLAYER14.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)
    @test RMS_differences(sim.θ, ref_NLAYER21.θ) < 0.04
    @test RMS_differences(sim.ψ, ref_NLAYER21.ψ) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.6,5.0)
    @test RMS_differences(sim.above[Not(1),[:time,:GWAT,:INTS,:INTR,:SNOW]],
                            ref_NLAYER21.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)
    @test RMS_differences(sim.θ, ref_NLAYER70.θ) < 0.04
    @test RMS_differences(sim.ψ, ref_NLAYER70.ψ) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.7,5.0)
    @test RMS_differences(sim.above[Not(1),[:time,:GWAT,:INTS,:INTR,:SNOW]],
                            ref_NLAYER70.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)

    # TODO(bernhard): we could run multiple LWFBrook90.jl simulations and compare with the
    # finest LWFBrook90R simulation only.

    # if (true)
    #     git_string = "git+"*chomp(Base.read(`git rev-parse --short HEAD`, String))*
    #     ifelse(length(Base.read(`git status --porcelain`, String))==0, "+clean","+dirty")

    #     # if some error appears, the following code can be used to plot the solutions
    #     using Plots, Measures
    #     pl_θ = plot(sim.θ.time,
    #             Matrix(sim.θ[:,Not(:time)]), line = :solid, labels = "LWFBrook90.jl:" .* string.(depth_to_read_out_mm) .* "mm",
    #             ylabel = "θ (-)", legend_position = :bottomright)
    #     plot!(Matrix(ref_NLAYER7.θ[:,Not(:time)]), line = :dash, color = :black,
    #             labels = "LWFBrook90R_NLayer7:" .* string.(depth_to_read_out_mm) .* "mm")
    #     # plot!(Matrix(ref_NLAYER14.θ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer14")
    #     # plot!(Matrix(ref_NLAYER21.θ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer21")
    #     # plot!(Matrix(ref_NLAYER70.θ[:,Not(:time)]), line = :dash, color = :black, labels = "LWFBrook90R_NLayer70")
    #     pl_ψ = plot(sim.ψ.time,
    #             Matrix(sim.ψ[:,Not(:time)]), line = :solid, labels = "LWFBrook90.jl:" .* string.(depth_to_read_out_mm) .* "mm",
    #             ylabel = "ψ (kPa)", legend_position = :bottomright)
    #     plot!(Matrix(ref_NLAYER7.ψ[:,Not(:time)]), line = :dash, color = :black,
    #             labels = "LWFBrook90R_NLayer7:" .* string.(depth_to_read_out_mm) .* "mm")
    #     # plot!(Matrix(ref_NLAYER14.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer14")
    #     # plot!(Matrix(ref_NLAYER21.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer21")
    #     # plot!(Matrix(ref_NLAYER70.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer70")
    #     pl_a = plot(sim.above.time,
    #             Matrix(#sim.above[:,Not(:time)]),
    #                    # label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"],
    #                     sim.above[:,[:GWAT,:INTS,:INTR,:SNOW]]),
    #                     label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"],
    #                     line = :solid)
    #     plot!(Matrix(ref_NLAYER7.above[:,[:intr,:ints,:snow,:gwat]]),
    #             line = :dash, color = :black, labels = "LWFBrook90R_NLayer7")
    #     # plot!(Matrix(ref_NLAYER14.above[:,[:intr,:ints,:snow,:gwat]]),
    #     #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
    #     # plot!(Matrix(ref_NLAYER21.above[:,[:intr,:ints,:snow,:gwat]]),
    #     #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
    #     # plot!(Matrix(ref_NLAYER70.above[:,[:intr,:ints,:snow,:gwat]]),
    #     #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
    #     plot(pl_θ, pl_ψ, pl_a, layout = (3,1), size = (600,800),
    #         leftmargin = 5mm)
    #     savefig(out_figure_string*git_string*".png")
    # end
end




# prepare plotting functions for Hammel (with/without δ)
function my_plot(df; args...)
    plot(df[:,:time], Matrix(df[:,Not(:time)]); args...)
end
function my_plot!(pl, df; args...)
    plot!(pl, df[:,:time], Matrix(df[:,Not(:time)]); args...)
end
function my_scatter!(pl, df; args...)
    scatter!(pl, df[:,:time], Matrix(df[:,Not(:time)]); args...)
end
function plot_Hammel_Dense(sim, ref, hyd, depth_to_read_out_mm, title; args...)
    pl_θ = plot(sim.ODEsolution.t, sim.θdense; line = :solid, color = [1 2 3 4 5], labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    my_scatter!(pl_θ, ref.θ;       line = :dot,   color = :black,      labels = ["LWFBrook90R" "" "" "" ""], markersize=2.5)
    # my_plot!(   pl_θ,hyd.θdense, line = :dash,  color = [1 2 3 4 5],               labels = nothing)
    my_scatter!(pl_θ,hyd.θdense, line = :dash,  color = [1 2 3 4 5], markersize=2, markerstrokecolor=:auto, labels = "Hydrus: " .* string.(depth_to_read_out_mm) .* " mm")
    pl_ψ = plot(sim.ODEsolution.t, sim.ψdense; line = :solid, color = [1 2 3 4 5], labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    my_scatter!(pl_ψ, ref.ψ;       line = :dot,   color = :black,      labels = ["LWFBrook90R" "" "" "" ""], markersize=2.5)
    # my_plot!(   pl_ψ,hyd.ψdense, line = :dash,  color = [1 2 3 4 5],               labels = nothing)
    my_scatter!(pl_ψ,hyd.ψdense, line = :dash,  color = [1 2 3 4 5], markersize=2, markerstrokecolor=:auto, labels = "Hydrus: " .* string.(depth_to_read_out_mm) .* " mm")

    simulate_isotopes = !isnothing(sim.δ18Odense)
    if (simulate_isotopes)
        pl_δ18O = plot(sim.ODEsolution.t, sim.δ18Odense; line = :solid, color = [1 2 3 4 5], labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
        # my_plot!(pl_δ18O, ref.δ18O; line = :dot,   color = :black,      labels = ["LWFBrook90R" "" "" "" ""], markersize=2.5)
        # my_plot!(   pl_δ18O,hyd.δ18Odense, line = :dash,  color = [1 2 3 4 5],               labels = nothing)
        my_scatter!(pl_δ18O,hyd.δ18Odense, line = :dash,  color = [1 2 3 4 5], markersize=2, markerstrokecolor=:auto, labels = "Hydrus: " .* string.(depth_to_read_out_mm) .* " mm")
        pl_δ2H = plot(sim.ODEsolution.t, sim.δ2Hdense; line = :solid, color = [1 2 3 4 5], labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
        # my_plot!(pl_δ2H, ref.δ2H;   line = :dot,   color = :black,      labels = ["LWFBrook90R" "" "" "" ""], markersize=2.5)
        # my_plot!(   pl_δ2H,hyd.δ2Hdense, line = :dash,  color = [1 2 3 4 5],               labels = nothing)
        my_scatter!(pl_δ2H,hyd.δ2Hdense, line = :dash,  color = [1 2 3 4 5], markersize=2, markerstrokecolor=:auto, labels = "Hydrus: " .* string.(depth_to_read_out_mm) .* " mm")
    else
        pl_δ18O = plot()
        pl_δ2H = plot()
    end

    ODE_method = replace(String(Symbol(sim.ODEsolution.alg)), r"(.*?)(\(|{).*"=>s"\1")
    plot(plot_title = Printf.@sprintf("%s \n %d time steps with method: %s",
            title,
            sim.ODEsolution.destats.naccept,
            ODE_method),
        filter(!isnothing,
                [plot!(pl_θ,ylabel = "θ (-)"),
                plot!(pl_ψ,ylabel = "ψ (kPa)"),
                ifelse(simulate_isotopes, plot!(pl_δ18O,ylabel = "δ18O (‰)"), nothing),
                ifelse(simulate_isotopes, plot!(pl_δ2H,ylabel = "δ2H (‰)"),   nothing)])... ,
        size = (1200,1200), layout = ifelse(simulate_isotopes, (4,1), (2,1)), leftmargin = 8mm;
        args...)
end
@testset "Hammel-2001-θ-ψ" begin
    # This test set checks that the RMSE of θ and ψ vs reference solutions are below hard-
    # coded limits. Two simulations defined in Hammel et al. (2001) are run.

    @show pwd() # This is to help get the folder right.

    # 1) Run simulation and load references (LWFBrook90R and Hydrus1D)
    @githash_time sim4, ref4, hyd4 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_sand-NLayer-27-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand");
        # benchmark reportings: "-git+"*chomp(Base.read(`git rev-parse --short HEAD`, String))
        # amberMBP-git-c4275ee: 0.913546 seconds (4.04 M allocations: 605.712 MiB, 14.49% gc time, 3.67% compilation time)
        # amberMBP-git-eae940b: 0.859826 seconds (1.85 M allocations: 499.565 MiB, 13.78% gc time)
        # amberMBP-git-b5cd0a6: 1.000577 seconds (2.51 M allocations: 476.224 MiB, 14.86% gc time)
        # amberMBP-git-a1872dd: 1.039936 seconds (2.43 M allocations: 533.845 MiB, 19.73% gc time)
        # amberMBP-git+356c4d6: 0.984474 seconds (2.44 M allocations: 542.654 MiB, 20.45% gc time) 11292 time steps
        # amberMBP-git-f94f0a3 Tsit5(): 1.385506 seconds (2.43 M allocations: 570.287 MiB, 11.90% gc time) 11293 time steps
        # amberMBP-git-f94f0a3 ImplicitEuler(): 1.612636 seconds (913.64 k allocations: 173.445 MiB) 707 time steps
        # amberMBP-git-f94f0a3 Rosenbrock23(): 0.653890 seconds (847.23 k allocations: 216.654 MiB) 140 Time steps
        # amberMBP-git-f94f0a3 TRBDF2(): 0.410966 seconds (528.38 k allocations: 109.498 MiB) 90 Time steps
        # amberMBP-git-f94f0a3 Rodas4P(): 0.734048 seconds (853.51 k allocations: 217.404 MiB) [ Info: Time steps for solving: 146
        # amberMBP-git-f94f0a3 Tsit5(); reltol = 1e-5: 2.450327 seconds (2.55 M allocations: 601.890 MiB, 49.53% gc time) 11466 Time steps
        # amberMBP-git-ed6bed2: 0.994531 seconds (2.34 M allocations: 502.295 MiB) 11466 time steps
    @githash_time sim5, ref5, hyd5 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_sand-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand");
        # amberMBP-git-c4275ee: 2.788329 seconds (3.92 M allocations: 1.817 GiB, 15.26% gc time)
        # amberMBP-git-eae940b: 2.626265 seconds (1.83 M allocations: 1.586 GiB, 13.83% gc time)
        # amberMBP-git-b5cd0a6: 2.885293 seconds (2.57 M allocations: 1.456 GiB, 13.65% gc time)
        # amberMBP-git-a1872dd: 2.878616 seconds (2.47 M allocations: 1.653 GiB, 14.48% gc time)
        # amberMBP-git+356c4d6: 2.944147 seconds (2.48 M allocations: 1.674 GiB, 15.51% gc time) 11102 time steps
        # amberMBP-git-f94f0a3 Tsit5(): 3.916065 seconds (2.39 M allocations: 1.714 GiB, 11.22% gc time) 11120 time steps
        # amberMBP-git-f94f0a3 ImplicitEuler(): 17.450994 seconds (6.16 M allocations: 4.640 GiB, 5.30% gc time) 884 Time steps
        # amberMBP-git-f94f0a3 Rosenbrock23(): 27.437288 seconds (48.15 M allocations: 6.042 GiB, 4.45% gc time, 63.21% compilation time) 228 Time steps
        # amberMBP-git-f94f0a3 TRBDF2(): 11.271603 seconds (4.08 M allocations: 3.026 GiB, 9.36% gc time) 232 time steps
        # amberMBP-git-f94f0a3 Rodas4P(): 11.553674 seconds (5.04 M allocations: 4.167 GiB, 8.95% gc time) [ Info: Time steps for solving: 277
        # amberMBP-git-f94f0a3 Tsit5(); reltol = 1e-5: 5.862678 seconds (2.54 M allocations: 1.809 GiB, 11.54% gc time) 11425
        # amberMBP-git-ed6bed2: 3.474783 seconds (2.33 M allocations: 1.473 GiB, 7.34% gc time) 11425 time steps
    # @githash_time sim6, ref6, hyd6 = prepare_θψδ_from_sim_and_reference(;
    #     path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_sand-NLayer-400-RESET=FALSE",
    #     path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
    #     path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand");
    #     # not run
    #     # amberMBP-git+356c4d6:31.183175 seconds (4.99 M allocations: 17.341 GiB, 15.47% gc time) 31536 time steps
    #     # amberMBP-git-f94f0a3 Tsit5():39.125138 seconds (6.78 M allocations: 17.645 GiB, 8.36% gc time) 31532 time steps
    #     # amberMBP-git-f94f0a3 ImplicitEuler()
    #     # amberMBP-git-f94f0a3 Rosenbrock23():
    #     # amberMBP-git-f94f0a3 TRBDF2()
    #     # amberMBP-git-f94f0a3 Rodas4P(): 212.835287 seconds (37.48 M allocations: 75.163 GiB, 6.17% gc time, 4.17% compilation time) [ Info: Time steps for solving: 410
    #     # amberMBP-git-f94f0a3 Tsit5(); reltol = 1e-5: 92.727500 seconds (14.08 M allocations: 18.600 GiB, 28.56% gc time, 8.44% compilation time) 31907
    #     # amberMBP-git-ed6bed2: 35.243902 seconds (6.51 M allocations: 14.777 GiB, 6.47% gc time) 31907 Time Steps

    @githash_time sim1, ref1, hyd1 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_loam-NLayer-27-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam");
        # amberMBP-git-c4275ee: 4.394856 seconds (19.52 M allocations: 2.898 GiB, 16.30% gc time)
        # amberMBP-git-eae940b: 3.882331 seconds (9.10 M allocations: 2.408 GiB, 15.04% gc time)
        # amberMBP-git-b5cd0a6: 4.210247 seconds (10.08 M allocations: 2.221 GiB, 13.75% gc time)
        # amberMBP-git-a1872dd: 4.158975 seconds (9.29 M allocations: 2.398 GiB, 12.88% gc time)
        # amberMBP-git+356c4d6: 3.963190 seconds (9.35 M allocations: 2.439 GiB, 15.06% gc time) 50969 time steps
        # amberMBP-git-f94f0a3 Tsit5(): 4.251475 seconds (8.21 M allocations: 1.880 GiB, 11.14% gc time) 34936 time steps (! note: unstable)
        # amberMBP-git-f94f0a3 ImplicitEuler(): 0.732530 seconds (572.87 k allocations: 114.522 MiB) 458 Time steps
        # amberMBP-git-f94f0a3 Rosenbrock23(): 3.134458 seconds (2.24 M allocations: 572.590 MiB, 26.69% gc time) 365 Time steps
        # amberMBP-git-f94f0a3 TRBDF2()
        # amberMBP-git-f94f0a3 Rodas4P(): 5.013004 seconds (6.10 M allocations: 1.512 GiB, 9.69% gc time) [ Info: Time steps for solving: 828
        # amberMBP-git-f94f0a3 Tsit5(); reltol = 1e-5: 11.977669 seconds (20.29 M allocations: 4.677 GiB, 9.69% gc time) 81248
        # amberMBP-git-ed6bed2: 9.474695 seconds (18.65 M allocations: 3.900 GiB, 15.40% gc time) 81248 time steps
    @githash_time sim2, ref2, hyd2 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_loam-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam");
        # amberMBP-git-c4275ee: 13.704865 seconds (19.92 M allocations: 9.220 GiB, 15.41% gc time, 0.23% compilation time)
        # amberMBP-git-eae940b: 12.581113 seconds (9.21 M allocations: 8.018 GiB, 13.85% gc time)
        # amberMBP-git-b5cd0a6: 13.095305 seconds (10.51 M allocations: 7.295 GiB, 12.65% gc time)
        # amberMBP-git-a1872dd: 13.194370 seconds (9.63 M allocations: 7.925 GiB, 12.67% gc time)
        # amberMBP-git+356c4d6: 12.782985 seconds (9.69 M allocations: 8.026 GiB, 13.64% gc time) 51980 time steps
        # amberMBP-git-f94f0a3 Tsit5(): 11.819875 seconds (7.05 M allocations: 5.037 GiB, 14.66% gc time) 29928 time steps (! note: unstable)
        # amberMBP-git-f94f0a3 ImplicitEuler(): 8.298325 seconds (2.43 M allocations: 1.809 GiB, 4.56% gc time) 546 Time steps
        # amberMBP-git-f94f0a3 Rosenbrock23(): 14.754657 seconds (6.27 M allocations: 5.186 GiB, 7.59% gc time) 329  time steps
        # amberMBP-git-f94f0a3 TRBDF2()
        # amberMBP-git-f94f0a3 Rodas4P():
        # amberMBP-git-f94f0a3 Tsit5(); reltol = 1e-5: 31.497258 seconds (18.67 M allocations: 13.284 GiB, 8.82% gc time) 75200 Time Steps
        # amberMBP-git-ed6bed2: 25.142349 seconds (17.15 M allocations: 10.808 GiB, 8.44% gc time), Ignoring NaN for fzero(), 75200 time steps
    # @githash_time sim3, ref3, hyd3 = prepare_θψδ_from_sim_and_reference(;
    #     path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_loam-NLayer-400-RESET=FALSE",
    #     path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
    #     path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam");
    #     # not run
    #     # amberMBP-git+356c4d6:70.575065 seconds (11.84 M allocations: 41.186 GiB, 15.96% gc time) 73252 time steps
    #     # amberMBP-git-f94f0a3 Tsit5():35.564689 seconds (6.32 M allocations: 16.418 GiB, 8.72% gc time) 27092 time steps (! note: unstable)
    #     # amberMBP-git-f94f0a3 ImplicitEuler()
    #     # amberMBP-git-f94f0a3 Rosenbrock23():
    #     # amberMBP-git-f94f0a3 TRBDF2(): 141.509080 seconds (12.71 M allocations: 36.789 GiB, 4.16% gc time) 250
    #     # amberMBP-git-f94f0a3 Rodas4P():
    #     # amberMBP-git-f94f0a3 Tsit5(); reltol = 1e-5: 124.500841 seconds (20.05 M allocations: 51.612 GiB, 15.12% gc time) 83352
    #     # amberMBP-git-ed6bed2: 104.234792 seconds (18.42 M allocations: 41.744 GiB, 9.21% gc time) 83352

    # 2) Plot (optional, not done when testing in CI)
    # Illustrate with a plot what will be compared in the tests below
    # if (true) # Do these manually outside of automatic testing in order not to require Plots pkg
    #     git_string = "git+"*chomp(Base.read(`git rev-parse --short HEAD`, String))*
    #         ifelse(length(Base.read(`git status --porcelain`, String))==0, "+clean","+dirty")
    #     using Plots, Measures
    #     depth_to_read_out_mm = [100 500 1000 1500 1900]
    #     # Plot Loam simulations
    #     plot_Hammel_Dense(sim1, ref1, hyd1, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900))
    #     savefig("test-assets/Hammel-2001/out_Loam_sim1_"*git_string*".png")
    #     plot_Hammel_Dense(sim2, ref2, hyd2, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900))
    #     savefig("test-assets/Hammel-2001/out_Loam_sim2_"*git_string*".png")
    #     plot_Hammel_Dense(sim3, ref3, hyd3, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900))
    #     savefig("test-assets/Hammel-2001/out_Loam_sim3_"*git_string*".png")

    #     # Plot Sand simulations
    #     plot_Hammel_Dense(sim4, ref4, hyd4, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900))
    #     savefig("test-assets/Hammel-2001/out_Sand_sim1_"*git_string*".png")
    #     plot_Hammel_Dense(sim5, ref5, hyd5, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900))
    #     savefig("test-assets/Hammel-2001/out_Sand_sim2_"*git_string*".png")
    #     plot_Hammel_Dense(sim6, ref6, hyd6, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900))
    #     savefig("test-assets/Hammel-2001/out_Sand_sim3_"*git_string*".png")
    # end

    # Use sensible accuracy values to compare the two solutions (e.g. θ of 0.02, and ψ of 1 kPa)
    # (e.g. RMSE(reference, simulated) < [hardcoded_value])
    # Floating point accuracy is used by regression tests.
    # (Or could be used, once I compare with a very high-fidelity solution, of which I known
    #  that my code should converge to. E.g. a high-resolution Hydrus simulation.)

    # 3) Test RMSE
    # Compare θ
    # Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.θ, ref1.θ) < 0.0020
    @test RMS_differences(sim2.θ, ref2.θ) < 0.00070
    # @test RMS_differences(sim3.θ, ref3.θ) < 0.00047
    @test RMS_differences(sim4.θ, ref4.θ) < 0.00035
    @test RMS_differences(sim5.θ, ref5.θ) < 0.00035
    # # @test RMS_differences(sim6.θ, ref6.θ) < 0.00035

    # Compare with Hydrus1D
    @test RMS_differences(sim1.θ[Not(end),:], hyd1.θ[Not(1),:]) < 0.005
    @test RMS_differences(sim2.θ[Not(end),:], hyd2.θ[Not(1),:]) < 0.002
    # @test RMS_differences(sim3.θ[Not(end),:], hyd3.θ[Not(1),:]) < 0.0015
    @test RMS_differences(sim4.θ[Not(end),:], hyd4.θ[Not(1),:]) < 0.008
    @test RMS_differences(sim5.θ[Not(end),:], hyd5.θ[Not(1),:]) < 0.009
    # # @test RMS_differences(sim6.θ[Not(end),:], hyd6.θ[Not(1),:]) < 0.007

    # Compare ψ
    # Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.ψ, ref1.ψ) < 2.0 # kPa
    @test RMS_differences(sim2.ψ, ref2.ψ) < 0.60 # kPa
    # @test RMS_differences(sim3.ψ, ref3.ψ) < 0.2 # kPa
    @test RMS_differences(sim4.ψ, ref4.ψ) < 0.004 # kPa
    @test RMS_differences(sim5.ψ, ref5.ψ) < 0.0080 # kPa
    # # @test RMS_differences(sim6.ψ, ref6.ψ) < 0.0025 # kPa

    # Compare with Hydrus1D
    @test RMS_differences(sim1.ψ[Not(end),:], hyd1.ψ[Not(1),:]) < 6.1 # kPa
    @test RMS_differences(sim2.ψ[Not(end),:], hyd2.ψ[Not(1),:]) < 2.8 # kPa
    # @test RMS_differences(sim3.ψ[Not(end),:], hyd3.ψ[Not(1),:]) < 1.2 # kPa
    @test RMS_differences(sim4.ψ[Not(end),:], hyd4.ψ[Not(1),:]) < 1.0 # kPa
    @test RMS_differences(sim5.ψ[Not(end),:], hyd5.ψ[Not(1),:]) < 1.5 # kPa
    # # @test RMS_differences(sim6.ψ[Not(end),:], hyd6.ψ[Not(1),:]) < 1.0 # kPa
end

# # NOTE: locally, i.e. not on CI system, one might need to do manually cd("test")
@testset "Hammel-2001-θ-ψ-δ" begin
    # This test set checks that the RMSE of θ and ψ vs reference solutions are below hard-
    # coded limits. Two simulations defined in Hammel et al. (2001) are run.

    @show pwd() # This is to help get the folder right.

    # # 1) Run simulation and load references (LWFBrook90R and Hydrus1D)
    @githash_time sim4, ref4, hyd4 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_sand-NLayer-27-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand_ISO2",
        simulate_isotopes   = true);
        # amberMBP-git-8806cd8: 5.731435 seconds (11.12 M allocations: 1.274 GiB, 7.61% gc time, 74.89% compilation time)
        # amberMBP-git-a3fc1a2: 10.796262 seconds (106.83 M allocations: 8.381 GiB, 16.10% gc time)
        # amberMBP-git-a1872dd: 1.637964 seconds (7.12 M allocations: 900.025 MiB, 20.29% gc time) 11292 time steps
        # amberMBP-git+356c4d6: 1.584915 seconds (6.83 M allocations: 892.419 MiB, 15.79% gc time) 11292 time steps
        # amberMBP-git-ce1fd2a: 1.895484 seconds (6.83 M allocations: 892.421 MiB, 15.90% gc time) 11292
        # amberMBP-git-f94f0a3 Tsit5(); reltol = 1e-5: 17.148244 seconds (15.65 M allocations: 1.339 GiB, 54.79% compilation time) 11466
        # amberMBP-git-ed6bed2: 2.690979 seconds (5.02 M allocations: 823.388 MiB, 4.83% gc time) 11466 time steps
    @githash_time sim5, ref5, hyd5 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_sand-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand_ISO2",
        simulate_isotopes   = true);
        # amberMBP-git-8806cd8: 6.073306 seconds (7.10 M allocations: 3.247 GiB, 13.81% gc time, 37.73% compilation time)
        # amberMBP-git-a3fc1a2: 36.013205 seconds (341.00 M allocations: 29.140 GiB, 16.44% gc time)
        # amberMBP-git-a1872dd: 4.106577 seconds (13.84 M allocations: 2.628 GiB, 16.45% gc time) 11102 time steps
        # amberMBP-git+356c4d6: 4.037970 seconds (14.20 M allocations: 2.656 GiB, 16.00% gc time) 11102 time steps
        # amberMBP-git-ce1fd2a: 4.264644 seconds (13.49 M allocations: 2.634 GiB, 15.07% gc time) 11120
        # amberMBP-git-f94f0a3 Tsit5(); reltol = 1e-5:  17.535086 seconds (22.57 M allocations: 3.206 GiB, 5.77% gc time, 44.25% compilation time) 11425
        # amberMBP-git-ed6bed2: 5.088062 seconds (5.01 M allocations: 2.316 GiB, 6.92% gc time) 11425
    # @githash_time sim6, ref6, hyd6 = prepare_θψδ_from_sim_and_reference(;
    #     path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_sand-NLayer-400-RESET=FALSE",
    #     path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
    #     path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand_ISO2",
    #     simulate_isotopes   = true);
    #     # not run
    #     # amberMBP-git+356c4d6:46.806037 seconds (114.21 M allocations: 27.069 GiB, 16.70% gc time) 31536 time steps
    #     # amberMBP-git-ce1fd2a: 50.951525 seconds (113.15 M allocations: 27.023 GiB, 17.87% gc time) 31527
    #     # amberMBP-git-f94f0a3 Tsit5(): 63.332243 seconds (130.71 M allocations: 28.506 GiB, 10.02% gc time, 15.88% compilation time) 31901 time steps
    #     # amberMBP-git-ed6bed2: 53.701377 seconds (26.34 M allocations: 23.542 GiB, 7.12% gc time, 13.83% compilation time) 31907

    @githash_time sim1, ref1, hyd1 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_loam-NLayer-27-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam_ISO2",
        simulate_isotopes   = true);
        # amberMBP-git-9d3342b non-adaptive 1e-3: 5.787269 seconds (23.98 M allocations: 4.513 GiB, 19.17% gc time)
        # amberMBP-git-9d3342b non-adaptive 1e-4: 58.027604 seconds (234.96 M allocations: 45.378 GiB, 18.09% gc time, 2.00% compilation time) -> unusable result containing NaN
        # amberMBP-git-a3fc1a2 adaptive (min 1e-4): 49.855293 seconds (492.71 M allocations: 38.667 GiB, 15.75% gc time, 0.06% compilation time)
        # amberMBP-git-a1872dd:adaptive_internalnorm: 7.128960 seconds (30.50 M allocations: 4.011 GiB, 14.55% gc time) 50969 time steps
        # amberMBP-git+356c4d6: 6.942920 seconds (31.24 M allocations: 4.048 GiB, 13.73% gc time) 50969 time steps
        # amberMBP-git-ce1fd2a: 7.501874 seconds (31.24 M allocations: 4.050 GiB, 13.77% gc time) 50942
        # amberMBP-git-f94f0a3 Tsit5(); reltol = 1e-5: 15.217210 seconds (57.18 M allocations: 7.344 GiB, 12.28% gc time) 81248
        # amberMBP-git-ed6bed2: 18.685461 seconds (37.66 M allocations: 6.121 GiB, 15.20% gc time) 81248
    @githash_time sim2, ref2, hyd2 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_loam-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam_ISO2",
        simulate_isotopes   = true);
        # amberMBP-git-9d3342b non-adaptive 1e-3: 15.853492 seconds (24.28 M allocations: 14.548 GiB, 21.53% gc time)
        # amberMBP-git-9d3342b non-adaptive 1e-4: 162.067400 seconds (233.84 M allocations: 146.493 GiB, 20.15% gc time) -> unusable result containing NaN
        # amberMBP-git-a3fc1a2 adaptive (min 1e-4): 172.691746 seconds (1.58 G allocations: 134.645 GiB, 16.29% gc time)
        # amberMBP-git-a1872dd:adaptive_internalnorm: 18.436145 seconds (62.86 M allocations: 12.487 GiB, 15.34% gc time) 51980 time steps
        # amberMBP-git+356c4d6: 22.144986 seconds (63.34 M allocations: 12.572 GiB, 14.51% gc time) 51980 time steps
        # amberMBP-git-ce1fd2a: 19.811971 seconds (63.28 M allocations: 12.569 GiB, 14.88% gc time) 51916
        # amberMBP-git-f94f0a3 Tsit5(); reltol = 1e-5: 37.237972 seconds (98.53 M allocations: 20.010 GiB, 10.34% gc time) 75200
        # amberMBP-git-ed6bed2: 34.159066 seconds (34.75 M allocations: 16.353 GiB, 9.60% gc time) 75200

    # @githash_time sim3, ref3, hyd3 = prepare_θψδ_from_sim_and_reference(;
    #     path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_loam-NLayer-400-RESET=FALSE",
    #     path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
    #     path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam_ISO2",
    #     simulate_isotopes   = true);
    #     # not run
    #     # amberMBP-git-a1872dd:adaptive_internalnorm: 107.662740 seconds (266.02 M allocations: 63.342 GiB, 17.37% gc time)
    #     # amberMBP-git+356c4d6:97.071439 seconds (263.17 M allocations: 63.688 GiB, 17.62% gc time) 73252 time steps
    #     # amberMBP-git-ce1fd2a: 108.155253 seconds (263.41 M allocations: 63.749 GiB, 17.42% gc time) 73318
    #     # amberMBP-git-f94f0a3 Tsit5(): 102.889782 seconds (212.54 M allocations: 52.711 GiB, 14.05% gc time) 57915
    #     # amberMBP-git-ed6bed2: 130.906236 seconds (38.01 M allocations: 62.973 GiB, 8.38% gc time) 83352

    # 2) Plot (optional, not done when testing in CI)
    # Illustrate with a plot what will be compared in the tests below
    # if (true) # Do these manually outside of automatic testing in order not to require Plots pkg
    #     git_string = "git+"*chomp(Base.read(`git rev-parse --short HEAD`, String))*
    #     ifelse(length(Base.read(`git status --porcelain`, String))==0, "+clean","+dirty")

    #     using Plots, Measures
    #     depth_to_read_out_mm = [100 500 1000 1500 1900]
    #     # Plot Loam simulations
    #     plot_Hammel_Dense(sim1, ref1, hyd1, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam", size=(900,900), legendfont=font(6))
    #     savefig("test-assets/Hammel-2001/out_Iso-Loam_sim1_"*git_string*".png")
    #     plot_Hammel_Dense(sim2, ref2, hyd2, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam", size=(900,900), legendfont=font(6))
    #     savefig("test-assets/Hammel-2001/out_Iso-Loam_sim2_"*git_string*".png")
    #     plot_Hammel_Dense(sim3, ref3, hyd3, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam", size=(900,900), legendfont=font(6))
    #     savefig("test-assets/Hammel-2001/out_Iso-Loam_sim3_"*git_string*".png")

    #     # Plot Sand simulations
    #     plot_Hammel_Dense(sim4, ref4, hyd4, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand", size=(900,900), legendfont=font(6))
    #     savefig("test-assets/Hammel-2001/out_Iso-Sand_sim1_"*git_string*".png")
    #     plot_Hammel_Dense(sim5, ref5, hyd5, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand", size=(900,900), legendfont=font(6))
    #     savefig("test-assets/Hammel-2001/out_Iso-Sand_sim2_"*git_string*".png")
    #     plot_Hammel_Dense(sim6, ref6, hyd6, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand", size=(900,900), legendfont=font(6))
    #     savefig("test-assets/Hammel-2001/out_Iso-Sand_sim3_"*git_string*".png")
    # end


    # Use sensible accuracy values to compare the two solutions (e.g. θ of 0.02, and ψ of 1 kPa)
    # (e.g. RMSE(reference, simulated) < [hardcoded_value])
    # Floating point accuracy is used by regression tests.
    # (Or could be used, once I compare with a very high-fidelity solution, of which I known
    #  that my code should converge to. E.g. a high-resolution Hydrus simulation.)

    # 3) Test RMSE
    # Compare θ
    # θ: Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.θ, ref1.θ) < 0.0018
    @test RMS_differences(sim2.θ, ref2.θ) < 0.00045
    # @test RMS_differences(sim3.θ, ref3.θ) < 0.00040
    @test RMS_differences(sim4.θ, ref4.θ) < 0.00040
    @test RMS_differences(sim5.θ, ref5.θ) < 0.00035
    # @test RMS_differences(sim6.θ, ref6.θ) < 0.00035

    # θ: Compare with Hydrus1D-Iso
    @test RMS_differences(sim1.θ[Not(end),:], hyd1.θ[Not(1),:]) < 0.005
    @test RMS_differences(sim2.θ[Not(end),:], hyd2.θ[Not(1),:]) < 0.002
    # @test RMS_differences(sim3.θ[Not(end),:], hyd3.θ[Not(1),:]) < 0.001
    @test RMS_differences(sim4.θ[Not(end),:], hyd4.θ[Not(1),:]) < 0.007
    @test RMS_differences(sim5.θ[Not(end),:], hyd5.θ[Not(1),:]) < 0.009
    # @test RMS_differences(sim6.θ[Not(end),:], hyd6.θ[Not(1),:]) < 0.007

    # Compare ψ
    # ψ: Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.ψ, ref1.ψ) < 1.6 # kPa
    @test RMS_differences(sim2.ψ, ref2.ψ) < 0.4 # kPa
    # @test RMS_differences(sim3.ψ, ref3.ψ) < 0.2 # kPa
    @test RMS_differences(sim4.ψ, ref4.ψ) < 0.003 # kPa
    @test RMS_differences(sim5.ψ, ref5.ψ) < 0.004 # kPa
    # @test RMS_differences(sim6.ψ, ref6.ψ) < 0.0035 # kPa

    # ψ: Compare with Hydrus1D-Iso
    @test RMS_differences(sim1.ψ[Not(end),:], hyd1.ψ[Not(1),:]) < 6 # kPa
    @test RMS_differences(sim2.ψ[Not(end),:], hyd2.ψ[Not(1),:]) < 2.8 # kPa
    # @test RMS_differences(sim3.ψ[Not(end),:], hyd3.ψ[Not(1),:]) < 1.6 # kPa
    @test RMS_differences(sim4.ψ[Not(end),:], hyd4.ψ[Not(1),:]) < 0.6 # kPa
    @test RMS_differences(sim5.ψ[Not(end),:], hyd5.ψ[Not(1),:]) < 1.1 # kPa
    # @test RMS_differences(sim6.ψ[Not(end),:], hyd6.ψ[Not(1),:]) < 0.5 # kPa

    # Compare δ18O with Hydrus1D-Iso
    # @test RMS_differences(sim1.δ18O[Not(end),:], hyd1.δ18O[Not(1),:]) < 0.5 # unit: ‰ # skipping as on Macbook it works in REPL, but bug in Pkg.test...
    @test RMS_differences(sim2.δ18O[Not(end),:], hyd2.δ18O[Not(1),:]) < 0.3 # unit: ‰
    # @test RMS_differences(sim3.δ18O[Not(end),:], hyd3.δ18O[Not(1),:]) < 0.5 # unit: ‰
    @test RMS_differences(sim4.δ18O[Not(end),:], hyd4.δ18O[Not(1),:]) < 0.5 # unit: ‰
    @test RMS_differences(sim5.δ18O[Not(end),:], hyd5.δ18O[Not(1),:]) < 0.4 # unit: ‰
    # @test RMS_differences(sim6.δ18O[Not(end),:], hyd6.δ18O[Not(1),:]) < 0.5 # unit: ‰
    # Compare δ2H with Hydrus1D-Iso
    @test_skip RMS_differences(sim1.δ2H[Not(end),:], hyd1.δ2H[Not(1),:]) < 3 # unit: ‰ # skipping as on Macbook it works in REPL, but bug in Pkg.test...
    @test RMS_differences(sim2.δ2H[Not(end),:], hyd2.δ2H[Not(1),:]) < 2 # unit: ‰
    # @test RMS_differences(sim3.δ2H[Not(end),:], hyd3.δ2H[Not(1),:]) < 2.9 # unit: ‰
    @test RMS_differences(sim4.δ2H[Not(end),:], hyd4.δ2H[Not(1),:]) < 2.3 # unit: ‰
    @test RMS_differences(sim5.δ2H[Not(end),:], hyd5.δ2H[Not(1),:]) < 1.65 # unit: ‰
    # @test RMS_differences(sim6.δ2H[Not(end),:], hyd6.δ2H[Not(1),:]) < 3.7 # unit: ‰

end

