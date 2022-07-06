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
@testset "Hammel-2001-θ-ψ" begin
    # This test set checks that the RMSE of θ and ψ vs reference solutions are below hard-
    # coded limits. Two simulations defined in Hammel et al. (2001) are run.

    @show pwd() # This is to help get the folder right.

    # 1) Run simulation and load references (LWFBrook90R and Hydrus1D)
    @githash_time sim4, ref4, hyd4 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_sand-NLayer-27-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand");
        # amberMBP-git-c4275ee: 0.913546 seconds (4.04 M allocations: 605.712 MiB, 14.49% gc time, 3.67% compilation time)
        # amberMBP-git-eae940b: 0.859826 seconds (1.85 M allocations: 499.565 MiB, 13.78% gc time)
        # amberMBP-git-b5cd0a6: 1.000577 seconds (2.51 M allocations: 476.224 MiB, 14.86% gc time)
    @githash_time sim5, ref5, hyd5 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_sand-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand");
        # amberMBP-git-c4275ee: 2.788329 seconds (3.92 M allocations: 1.817 GiB, 15.26% gc time)
        # amberMBP-git-eae940b: 2.626265 seconds (1.83 M allocations: 1.586 GiB, 13.83% gc time)
        # amberMBP-git-b5cd0a6: 2.885293 seconds (2.57 M allocations: 1.456 GiB, 13.65% gc time)
    # @githash_time sim6, ref6, hyd6 = prepare_θψδ_from_sim_and_reference(;
    #     path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_sand-NLayer-400-RESET=FALSE",
    #     path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
    #     path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand")
        # not run

    @githash_time sim1, ref1, hyd1 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_loam-NLayer-27-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam");
        # amberMBP-git-c4275ee: 4.394856 seconds (19.52 M allocations: 2.898 GiB, 16.30% gc time)
        # amberMBP-git-eae940b: 3.882331 seconds (9.10 M allocations: 2.408 GiB, 15.04% gc time)
        # amberMBP-git-b5cd0a6: 4.210247 seconds (10.08 M allocations: 2.221 GiB, 13.75% gc time)
    @githash_time sim2, ref2, hyd2 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_loam-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam");
        # amberMBP-git-c4275ee: 13.704865 seconds (19.92 M allocations: 9.220 GiB, 15.41% gc time, 0.23% compilation time)
        # amberMBP-git-eae940b: 12.581113 seconds (9.21 M allocations: 8.018 GiB, 13.85% gc time)
        # amberMBP-git-b5cd0a6: 13.095305 seconds (10.51 M allocations: 7.295 GiB, 12.65% gc time)
    # @githash_time sim3, ref3, hyd3 = prepare_θψδ_from_sim_and_reference(;
    #     path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_loam-NLayer-400-RESET=FALSE",
    #     path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
    #     path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam")
        # not run

    # 2) Plot (optional, not done when testing in CI)
    # # Illustrate with a plot what will be compared in the tests below
    # using Plots, Measures
    # function my_plot(df; args...)
    #     plot(df[:,:time], Matrix(df[:,Not(:time)]); args...)
    # end
    # function my_plot!(pl, df; args...)
    #     plot!(pl, df[:,:time], Matrix(df[:,Not(:time)]); args...)
    # end
    # depth_to_read_out_mm = [100 500 1000 1500 1900]
    # ## Go for simulation 2
    # pl2_θ = my_plot(sim2.θ;      line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    # my_plot!(pl2_θ, ref2.θ;      line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""])
    # my_plot!(pl2_θ, hyd2.θ;      line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""])
    # # my_plot!(pl2_θ, hyd2.θdense,line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""])
    # pl2_ψ = my_plot(sim2.ψ; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    # my_plot!(pl2_ψ, ref2.ψ; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""])
    # my_plot!(pl2_ψ, hyd2.ψ; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""])
    # # pl2_δ18O = my_plot(sim2.δ18O; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    # # my_plot!(pl2_δ18O, ref2.δ18O; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""])
    # # my_plot!(pl2_δ18O, hyd2.δ18O; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""])
    # # pl2_δ2H = my_plot(sim2.δ2H; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    # # my_plot!(pl2_δ2H, ref2.δ2H; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""])
    # # my_plot!(pl2_δ2H, hyd2.δ2H; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""])
    # plot(title = "Simulation from Hammel et al. (2001) - Loam",
    #      plot!(pl2_θ,ylabel = "θ (-)"),
    #      plot!(pl2_ψ,ylabel = "ψ (kPa)"),
    #      #pl2_δ18O, pl2_δ2H,
    #      size = (1200,1200), layout = (2,1), leftmargin = 8mm)
    # savefig("test-assets/Hammel-2001/out_Loam.png")
    # ## Go for simulation 5
    # pl5_θ = my_plot(sim5.θ;      line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    # my_plot!(pl5_θ, ref5.θ;      line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""])
    # my_plot!(pl5_θ, hyd5.θ;      line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""])
    # # my_plot!(pl5_θ, hyd5.θdense,line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""])
    # pl5_ψ = my_plot(sim5.ψ; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    # my_plot!(pl5_ψ, ref5.ψ; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""])
    # my_plot!(pl5_ψ, hyd5.ψ; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""])
    # # pl5_δ18O = my_plot(sim5.δ18O; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    # # my_plot!(pl5_δ18O, ref5.δ18O; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""])
    # # my_plot!(pl5_δ18O, hyd5.δ18O; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""])
    # # pl5_δ2H = my_plot(sim5.δ2H; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    # # my_plot!(pl5_δ2H, ref5.δ2H; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""])
    # # my_plot!(pl5_δ2H, hyd5.δ2H; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""])
    # plot(title = "Simulation from Hammel et al. (2001) - Sand",
    #      plot!(pl5_θ,ylabel = "θ (-)"),
    #      plot!(pl5_ψ,ylabel = "ψ (kPa)"),
    #      #pl5_δ18O, pl5_δ2H,
    #      size = (1200,1200), layout = (2,1), leftmargin = 8mm)
    # savefig("test-assets/Hammel-2001/out_Sand.png")

    # Use sensible accuracy values to compare the two solutions (e.g. θ of 0.02, and ψ of 1 kPa)
    # (e.g. RMSE(reference, simulated) < [hardcoded_value])
    # Floating point accuracy is used by regression tests.
    # (Or could be used, once I compare with a very high-fidelity solution, of which I known
    #  that my code should converge to. E.g. a high-resolution Hydrus simulation.)

    # 3) Test RMSE
    # Compare θ
    # Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.θ, ref1.θ) < 0.0015
    @test RMS_differences(sim2.θ, ref2.θ) < 0.00040
    # @test RMS_differences(sim3.θ, ref3.θ) < 0.00045
    @test RMS_differences(sim4.θ, ref4.θ) < 0.00035
    @test RMS_differences(sim5.θ, ref5.θ) < 0.00035
    # # @test RMS_differences(sim6.θ, ref6.θ) < 0.00035

    # Compare with Hydrus1D
    @test RMS_differences(sim1.θ[Not(end),:], hyd1.θ[Not(1),:]) < 0.005
    @test RMS_differences(sim2.θ[Not(end),:], hyd2.θ[Not(1),:]) < 0.002
    # @test RMS_differences(sim3.θ[Not(end),:], hyd3.θ[Not(1),:]) < 0.001
    @test RMS_differences(sim4.θ[Not(end),:], hyd4.θ[Not(1),:]) < 0.006
    @test RMS_differences(sim5.θ[Not(end),:], hyd5.θ[Not(1),:]) < 0.007
    # # @test RMS_differences(sim6.θ[Not(end),:], hyd6.θ[Not(1),:]) < 0.005

    # Compare ψ
    # Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.ψ, ref1.ψ) < 2.0 # kPa
    @test RMS_differences(sim2.ψ, ref2.ψ) < 0.26 # kPa
    # @test RMS_differences(sim3.ψ, ref3.ψ) < 0.2 # kPa
    @test RMS_differences(sim4.ψ, ref4.ψ) < 0.003 # kPa
    @test RMS_differences(sim5.ψ, ref5.ψ) < 0.004 # kPa
    # # @test RMS_differences(sim6.ψ, ref6.ψ) < 0.0025 # kPa

    # Compare with Hydrus1D
    @test RMS_differences(sim1.ψ[Not(end),:], hyd1.ψ[Not(1),:]) < 6 # kPa
    @test RMS_differences(sim2.ψ[Not(end),:], hyd2.ψ[Not(1),:]) < 2.6 # kPa
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
        # amberMBP-git-TODO:    1.956504 seconds (7.55 M allocations: 909.316 MiB, 15.30% gc time)
    @githash_time sim5, ref5, hyd5 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_sand-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand_ISO2",
        simulate_isotopes   = true);
        # amberMBP-git-8806cd8: 6.073306 seconds (7.10 M allocations: 3.247 GiB, 13.81% gc time, 37.73% compilation time)
        # amberMBP-git-a3fc1a2: 36.013205 seconds (341.00 M allocations: 29.140 GiB, 16.44% gc time)
        # amberMBP-git-TODO:    4.430388 seconds (14.26 M allocations: 2.638 GiB, 15.92% gc time)
    # @githash_time sim6, ref6, hyd6 = prepare_θψδ_from_sim_and_reference(;
    #     path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_sand-NLayer-400-RESET=FALSE",
    #     path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_sand-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
    #     path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Sand_ISO2",
    #     simulate_isotopes   = true);
        # not run

    @githash_time sim1, ref1, hyd1 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_loam-NLayer-27-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam_ISO2",
        simulate_isotopes   = true);
        # amberMBP-git-9d3342b non-adaptive 1e-3: 5.787269 seconds (23.98 M allocations: 4.513 GiB, 19.17% gc time)
        # amberMBP-git-9d3342b non-adaptive 1e-4: 58.027604 seconds (234.96 M allocations: 45.378 GiB, 18.09% gc time, 2.00% compilation time) -> unusable result containing NaN
        # amberMBP-git-a3fc1a2 adaptive (min 1e-4): 49.855293 seconds (492.71 M allocations: 38.667 GiB, 15.75% gc time, 0.06% compilation time)
        # amberMBP-git-TODO:   adaptive_internalnorm: 7.167678 seconds (32.53 M allocations: 4.054 GiB, 13.37% gc time)
    @githash_time sim2, ref2, hyd2 = prepare_θψδ_from_sim_and_reference(;
        path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_loam-NLayer-103-RESET=FALSE",
        path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-103-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
        path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam_ISO2",
        simulate_isotopes   = true);
        # amberMBP-git-9d3342b non-adaptive 1e-3: 15.853492 seconds (24.28 M allocations: 14.548 GiB, 21.53% gc time)
        # amberMBP-git-9d3342b non-adaptive 1e-4: 162.067400 seconds (233.84 M allocations: 146.493 GiB, 20.15% gc time) -> unusable result containing NaN
        # amberMBP-git-a3fc1a2 adaptive (min 1e-4): 172.691746 seconds (1.58 G allocations: 134.645 GiB, 16.29% gc time)
        # amberMBP-git-TODO:   adaptive_internalnorm: 19.682619 seconds (64.91 M allocations: 12.537 GiB, 15.37% gc time)
    # @githash_time sim3, ref3, hyd3 = prepare_θψδ_from_sim_and_reference(;
    #     path_jl_prefix      = "test-assets/Hammel-2001/input-files-ISO/Hammel_loam-NLayer-400-RESET=FALSE",
    #     path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-400-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv",
    #     path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam_ISO2",
    #     simulate_isotopes   = true);
    #     # not run
    #     # amberMBP-git-TODO:   adaptive_internalnorm: 107.662740 seconds (266.02 M allocations: 63.342 GiB, 17.37% gc time)

    # 2) Plot (optional, not done when testing in CI)
    # Illustrate with a plot what will be compared in the tests below
    if (false) # Do these manually outside of automatic testing in order not to require Plots pkg
        git_string = "git+"*chomp(read(`git rev-parse --short HEAD`, String))*
        ifelse(length(read(`git status --porcelain`, String))==0, "+clean","+dirty")

        using Plots, Measures
        function my_plot(df; args...)
            plot(df[:,:time], Matrix(df[:,Not(:time)]); args...)
        end
        function my_plot!(pl, df; args...)
            plot!(pl, df[:,:time], Matrix(df[:,Not(:time)]); args...)
        end
        depth_to_read_out_mm = [100 500 1000 1500 1900]
        ## Go for simulation 2
        pl2_θ = my_plot(sim2.θ;      line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm");
        my_plot!(pl2_θ, ref2.θ;      line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""]);
        my_plot!(pl2_θ, hyd2.θ;      line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""]);
        #my_plot!(pl2_θ,hyd2.θdense,line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""]);
        pl2_ψ = my_plot(sim2.ψ; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm");
        my_plot!(pl2_ψ, ref2.ψ; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""]);
        my_plot!(pl2_ψ, hyd2.ψ; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""]);
        pl2_δ18O = my_plot(sim2.δ18O; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm");
        #my_plot!(pl2_δ18O,ref2.δ18O; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""]);
        my_plot!(pl2_δ18O, hyd2.δ18O; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""]);
        pl2_δ2H = my_plot( sim2.δ2H; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm");
        #my_plot!(pl2_δ2H, ref2.δ2H; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""])
        my_plot!(pl2_δ2H,  hyd2.δ2H; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""])
        plot(title = "Simulation from Hammel et al. (2001) - Loam",
            plot!(pl2_θ,ylabel = "θ (-)"),
            plot!(pl2_ψ,ylabel = "ψ (kPa)"),
            plot!(pl2_δ18O,ylabel = "δ18O (‰)"),
            plot!(pl2_δ2H,ylabel = "δ2H (‰)"),
            size = (1200,1200), layout = (4,1), leftmargin = 8mm)
        savefig("test-assets/Hammel-2001/out_Iso-Loam_"*git_string*".png")
        ## Go for simulation 5
        pl5_θ = my_plot(sim5.θ;      line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm");
        my_plot!(pl5_θ, ref5.θ;      line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""]);
        my_plot!(pl5_θ, hyd5.θ;      line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""]);
        # my_plot!(pl5_θ, hyd5.θdense,line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""]);
        pl5_ψ = my_plot(sim5.ψ; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm");
        my_plot!(pl5_ψ, ref5.ψ; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""]);
        my_plot!(pl5_ψ, hyd5.ψ; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""]);
        pl5_δ18O = my_plot(sim5.δ18O; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm");
        #my_plot!(pl5_δ18O,ref5.δ18O; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""]);
        my_plot!(pl5_δ18O, hyd5.δ18O; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""]);
        pl5_δ2H = my_plot(sim5.δ2H; line = :solid,                labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm");
        #my_plot!(pl5_δ2H,ref5.δ2H; line = :dot, color = :black, labels = ["LWFBrook90R" "" "" "" ""]);
        my_plot!(pl5_δ2H, hyd5.δ2H; line = :dash, color = :green, labels = ["Hydrus" "" "" "" ""]);
        plot(title = "Simulation from Hammel et al. (2001) - Sand",
            plot!(pl5_θ,ylabel = "θ (-)"),
            plot!(pl5_ψ,ylabel = "ψ (kPa)"),
            plot!(pl5_δ18O,ylabel = "δ18O (‰)"),
            plot!(pl5_δ2H,ylabel = "δ2H (‰)"),
            size = (1200,1200), layout = (4,1), leftmargin = 8mm)
        savefig("test-assets/Hammel-2001/out_Iso-Sand_"*git_string*".png")
    end

    # Use sensible accuracy values to compare the two solutions (e.g. θ of 0.02, and ψ of 1 kPa)
    # (e.g. RMSE(reference, simulated) < [hardcoded_value])
    # Floating point accuracy is used by regression tests.
    # (Or could be used, once I compare with a very high-fidelity solution, of which I known
    #  that my code should converge to. E.g. a high-resolution Hydrus simulation.)

    # 3) Test RMSE
    # Compare θ
    # θ: Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.θ, ref1.θ) < 0.0015
    @test RMS_differences(sim2.θ, ref2.θ) < 0.00035
    # @test RMS_differences(sim3.θ, ref3.θ) < 0.00045
    @test RMS_differences(sim4.θ, ref4.θ) < 0.00035
    @test RMS_differences(sim5.θ, ref5.θ) < 0.00035
    # @test RMS_differences(sim6.θ, ref6.θ) < 0.00035

    # θ: Compare with Hydrus1D-Iso
    @test RMS_differences(sim1.θ[Not(end),:], hyd1.θ[Not(1),:]) < 0.005
    @test RMS_differences(sim2.θ[Not(end),:], hyd2.θ[Not(1),:]) < 0.002
    # @test RMS_differences(sim3.θ[Not(end),:], hyd3.θ[Not(1),:]) < 0.001
    @test RMS_differences(sim4.θ[Not(end),:], hyd4.θ[Not(1),:]) < 0.006
    @test RMS_differences(sim5.θ[Not(end),:], hyd5.θ[Not(1),:]) < 0.007
    # @test RMS_differences(sim6.θ[Not(end),:], hyd6.θ[Not(1),:]) < 0.005

    # Compare ψ
    # ψ: Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1.ψ, ref1.ψ) < 2.0 # kPa
    @test RMS_differences(sim2.ψ, ref2.ψ) < 0.2 # kPa
    # @test RMS_differences(sim3.ψ, ref3.ψ) < 0.2 # kPa
    @test RMS_differences(sim4.ψ, ref4.ψ) < 0.003 # kPa
    @test RMS_differences(sim5.ψ, ref5.ψ) < 0.004 # kPa
    # @test RMS_differences(sim6.ψ, ref6.ψ) < 0.0025 # kPa

    # ψ: Compare with Hydrus1D-Iso
    @test RMS_differences(sim1.ψ[Not(end),:], hyd1.ψ[Not(1),:]) < 6 # kPa
    @test RMS_differences(sim2.ψ[Not(end),:], hyd2.ψ[Not(1),:]) < 2.2 # kPa
    # @test RMS_differences(sim3.ψ[Not(end),:], hyd3.ψ[Not(1),:]) < 1.2 # kPa
    @test RMS_differences(sim4.ψ[Not(end),:], hyd4.ψ[Not(1),:]) < 1.0 # kPa
    @test RMS_differences(sim5.ψ[Not(end),:], hyd5.ψ[Not(1),:]) < 1.5 # kPa
    # @test RMS_differences(sim6.ψ[Not(end),:], hyd6.ψ[Not(1),:]) < 1.0 # kPa

    # TODO(bernhard): fix LWFBrook90.jl simulation 1 and 2
    # Compare δ18O with Hydrus1D-Iso
    @test RMS_differences(sim1.δ18O[Not(end),:], hyd1.δ18O[Not(1),:]) < 5 # unit: ‰
    @test RMS_differences(sim2.δ18O[Not(end),:], hyd2.δ18O[Not(1),:]) < 2.2 # unit: ‰
    # @test RMS_differences(sim3.δ18O[Not(end),:], hyd3.δ18O[Not(1),:]) < 1.2 # unit: ‰
    @test RMS_differences(sim4.δ18O[Not(end),:], hyd4.δ18O[Not(1),:]) < 0.5 # unit: ‰
    @test RMS_differences(sim5.δ18O[Not(end),:], hyd5.δ18O[Not(1),:]) < 1.0 # unit: ‰
    # @test RMS_differences(sim6.δ18O[Not(end),:], hyd6.δ18O[Not(1),:]) < 1.0 # unit: ‰
    # Compare δ2H with Hydrus1D-Iso
    @test RMS_differences(sim1.δ2H[Not(end),:], hyd1.δ2H[Not(1),:]) < 5 # unit: ‰
    @test RMS_differences(sim2.δ2H[Not(end),:], hyd2.δ2H[Not(1),:]) < 7 # unit: ‰
    # @test RMS_differences(sim3.δ2H[Not(end),:], hyd3.δ2H[Not(1),:]) < 1.2 # unit: ‰
    @test RMS_differences(sim4.δ2H[Not(end),:], hyd4.δ2H[Not(1),:]) < 2.0 # unit: ‰
    @test RMS_differences(sim5.δ2H[Not(end),:], hyd5.δ2H[Not(1),:]) < 10. # unit: ‰
    # @test RMS_differences(sim6.δ2H[Not(end),:], hyd6.δ2H[Not(1),:]) < 1.0 # unit: ‰

end


# NOTE: locally, i.e. not on CI system, one might need to do manually cd("test")
@testset "BEA-2016-θ-ψ-aboveground-states" begin

    # @show pwd() # This is to help get the folder right.

    # # Check the RMSE of θ in simulations is below a limit
    sim, ref_NLAYER7, ref_NLAYER14, ref_NLAYER21, ref_NLAYER70 =
        prepare_sim_and_ref_for_BEA_2016("test-assets/BEA-2016","BEA2016-reset-FALSE")

    # # Use sensible accuracy values to compare the two solutions (e.g. θ of 0.02, and ψ of 1 kPa)
    # # (e.g. RMSE(reference, simulated) < [hardcoded_value])
    # # Floating point accuracy is used by regression tests.
    # # (Or could be used, once I compare with a very high-fidelity solution, of which I known
    # #  that my code should converge to. E.g. a high-resolution Hydrus simulation.)

    # Compare with LWFBrook90R as reference solution
    # θ:
    @test RMS_differences(sim.θ, ref_NLAYER7.θ) < 0.007
    # ψ:
    @test RMS_differences(sim.ψ, ref_NLAYER7.ψ) < 0.07
    # "GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)", (not done for: "CC (MJ/m2)" "SNOWLQ (mm)"]):
    @test RMS_differences(sim.above[Not(1),[:time,:GWAT, :INTS, :INTR, :SNOW]],
                            ref_NLAYER7.above[Not(end),[:time, :gwat,:ints,:intr,:snow]]) < 0.0005

    # Note that below we compare a NLAYER7 LWFBrook90.jl solution, with finer resolved
    # LWFBrook90R solutions. It is therefore normal, that the uncertainty can increase...
    @test RMS_differences(sim.θ, ref_NLAYER14.θ) < 0.03
    @test RMS_differences(sim.ψ, ref_NLAYER14.ψ) < 0.6
    @test RMS_differences(sim.above[Not(1),[:time,:GWAT,:INTS,:INTR,:SNOW]],
                            ref_NLAYER14.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < 0.0005
    @test RMS_differences(sim.θ, ref_NLAYER21.θ) < 0.04
    @test RMS_differences(sim.ψ, ref_NLAYER21.ψ) < 0.6
    @test RMS_differences(sim.above[Not(1),[:time,:GWAT,:INTS,:INTR,:SNOW]],
                            ref_NLAYER21.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < 0.0005
    @test RMS_differences(sim.θ, ref_NLAYER70.θ) < 0.04
    @test RMS_differences(sim.ψ, ref_NLAYER70.ψ) < 0.7
    @test RMS_differences(sim.above[Not(1),[:time,:GWAT,:INTS,:INTR,:SNOW]],
                            ref_NLAYER70.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < 0.0005

    # TODO(bernhard): we could run multiple LWFBrook90.jl simulations and compare with the
    # finest LWFBrook90R simulation only.

    # # if some error appears, the following code can be used to plot the solutions
    # using Plots
    #     # Compare with LWFBrook90R as reference solution
    # # θ:
    # @test RMS_differences(sim.θ, ref_NLAYER7.θ) < 0.007
    # # ψ:
    # @test RMS_differences(sim.ψ, ref_NLAYER7.ψ) < 0.07
    # # "GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)", (not done for: "CC (MJ/m2)" "SNOWLQ (mm)"]):
    # @test RMS_differences(sim.above[Not(1),[:time,:GWAT, :INTS, :INTR, :SNOW]],
    #                         ref_NLAYER7.above[Not(end),[:time, :gwat,:ints,:intr,:snow]]) < 0.0005

    # pl_θ = plot(sim.θ.time,
    #         Matrix(sim.θ[:,Not(:time)]), line = :solid, labels = "LWFBrook90.jl",
    #         ylabel = "θ (-)")
    # plot!(Matrix(ref_NLAYER7.θ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
    # # plot!(Matrix(ref_NLAYER14.θ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer14")
    # # plot!(Matrix(ref_NLAYER21.θ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer21")
    # # plot!(Matrix(ref_NLAYER70.θ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer70")
    # pl_ψ = plot(sim.ψ.time,
    #         Matrix(sim.ψ[:,Not(:time)]), line = :solid, labels = "LWFBrook90.jl",
    #         ylabel = "ψ (kPa)")
    # plot!(Matrix(ref_NLAYER7.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
    # # plot!(Matrix(ref_NLAYER14.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer14")
    # # plot!(Matrix(ref_NLAYER21.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer21")
    # # plot!(Matrix(ref_NLAYER70.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer70")
    # pl_a = plot(sim.above.time,
    #         Matrix(sim.above[:,Not(:time)]), line = :solid, label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])
    # plot!(Matrix(ref_NLAYER7.above[:,[:intr,:ints,:snow,:gwat]]),
    #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
    # # plot!(Matrix(ref_NLAYER14.above[:,[:intr,:ints,:snow,:gwat]]),
    # #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
    # # plot!(Matrix(ref_NLAYER21.above[:,[:intr,:ints,:snow,:gwat]]),
    # #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
    # # plot!(Matrix(ref_NLAYER70.above[:,[:intr,:ints,:snow,:gwat]]),
    # #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
    # plot(pl_θ, pl_ψ, pl_a, layout = (3,1), size = (600,800))
    # savefig("test-assets/BEA-2016/out.png")

end
