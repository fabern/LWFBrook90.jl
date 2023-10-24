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
    @test RMS_differences(sim.above[:,[:time, :GWAT_mm, :INTS_mm, :INTR_mm, :SNOW_mm]],
                            ref_NLAYER7.above[Not(end),[:time, :gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)

    # Note that below we compare a NLAYER7 LWFBrook90.jl solution, with finer resolved
    # LWFBrook90R solutions. It is therefore normal, that the uncertainty can increase...
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"θ_")], ref_NLAYER14.θ) < 0.03
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"ψ_")], ref_NLAYER14.ψ) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.65,5.0)
    @test RMS_differences(sim.above[:,[:time, :GWAT_mm, :INTS_mm, :INTR_mm, :SNOW_mm]],
                            ref_NLAYER14.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"θ_")], ref_NLAYER21.θ) < 0.04
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"ψ_")], ref_NLAYER21.ψ) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.6,5.0)
    @test RMS_differences(sim.above[:,[:time, :GWAT_mm, :INTS_mm, :INTR_mm, :SNOW_mm]],
                            ref_NLAYER21.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"θ_")], ref_NLAYER70.θ) < 0.04
    @test RMS_differences(sim.θψδ[:,Cols(:time, r"ψ_")], ref_NLAYER70.ψ) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.7,5.0)
    @test RMS_differences(sim.above[:,[:time, :GWAT_mm, :INTS_mm, :INTR_mm, :SNOW_mm]],
                            ref_NLAYER70.above[Not(end),[:time,:gwat,:ints,:intr,:snow]]) < ifelse(input_prefix=="BEA2016-reset-FALSE",0.51,1.8)

    # TODO(bernhard): we could run multiple LWFBrook90.jl simulations and compare with the
    # finest LWFBrook90R simulation only.

    if !is_a_CI_system && plot_flag
        # fname_illustrations = "out/$(today())/"
        fname_illustrations = "out/$git_status_string/"
        mkpath(dirname(fname_illustrations))
        # if some error appears, the following code can be used to plot the solutions
        # using Plots, Measures

        # pl_θ = Plots.plot(sim.θψδ.time, Matrix(sim.θψδ[:, r"^θ"]),
        #         line = :solid, labels = reshape("LWFBrook90.jl:" .* names(sim.θψδ[:, r"^θ"]), 1,:),
        #         ylabel = "θ (-)", legend_position = :bottomright)
        # Plots.plot!(Matrix(ref_NLAYER7.θ[:,Not(:time)]), line = :dash, color = :black,
        #         labels = reshape("LWFBrook90R_NLayer7:" .* names(ref_NLAYER7.θ[:,Not(:time)]), 1,:))
        # # Plots.plot!(Matrix(ref_NLAYER14.θ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer14")
        # # Plots.plot!(Matrix(ref_NLAYER21.θ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer21")
        # # Plots.plot!(Matrix(ref_NLAYER70.θ[:,Not(:time)]), line = :dash, color = :black, labels = "LWFBrook90R_NLayer70")

        # sim
        # pl_ψ = Plots.plot(sim.θψδ.time,
        #         Matrix(sim.θψδ[:,r"^ψ"]),
        #         line = :solid, labels = reshape("LWFBrook90.jl:" .* names(sim.θψδ[:,r"^ψ"]), 1,:),
        #         ylabel = "ψ (kPa)", legend_position = :bottomright)
        # Plots.plot!(Matrix(ref_NLAYER7.ψ[:,Not(:time)]), line = :dash, color = :black,
        #         labels = reshape("LWFBrook90R_NLayer7:" .* names(ref_NLAYER7.ψ[:,Not(:time)]), 1,:))
        # # Plots.plot!(Matrix(ref_NLAYER14.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer14")
        # # Plots.plot!(Matrix(ref_NLAYER21.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer21")
        # # Plots.plot!(Matrix(ref_NLAYER70.ψ[:,Not(:time)]), line = :dot, color = :black, labels = "LWFBrook90R_NLayer70")
        # pl_a = Plots.plot(sim.above.time,
        #         Matrix(#sim.above[:,Not(:time)]),
        #                # label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"],
        #                 sim.above[:,[ :GWAT_mm, :INTS_mm, :INTR_mm, :SNOW_mm]]),
        #                 label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"],
        #                 line = :solid)
        # Plots.plot!(Matrix(ref_NLAYER7.above[:,[:intr,:ints,:snow,:gwat]]),
        #         line = :dash, color = :black, labels = "LWFBrook90R_NLayer7")
        # # Plots.plot!(Matrix(ref_NLAYER14.above[:,[:intr,:ints,:snow,:gwat]]),
        # #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
        # # Plots.plot!(Matrix(ref_NLAYER21.above[:,[:intr,:ints,:snow,:gwat]]),
        # #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
        # # Plots.plot!(Matrix(ref_NLAYER70.above[:,[:intr,:ints,:snow,:gwat]]),
        # #         line = :dot, color = :black, labels = "LWFBrook90R_NLayer7")
        # Plots.plot(pl_θ, pl_ψ, pl_a, layout = (3,1), size = (600,800),
        #     leftmargin = 5mm)
        # savefig(fname_illustrations*"TESTSET_site-θ-ψ-aboveground-states_$(site).png")


        # Make publication figure with Makie # https://docs.makie.org/stable/explanations/figure_size/#vector_graphics
        # size_inches = (3.60, 4*3) # often either 7.25 inches or 3.60 inches wide
        size_inches = (7.25, 12) # often either 7.25 inches or 3.60 inches wide
        size_pt = 72 .* size_inches
        fig = Makie.Figure(resolution = size_pt, fontsize = 12)
        ax1 = Makie.Axis(fig[1,1], xlabel = "Time (days)", ylabel = "θ (-)")#,   title = "Title")
        ax2 = Makie.Axis(fig[2,1], xlabel = "Time (days)", ylabel = "ψ (kPa)",#, title = "Title")
            yscale = Makie.pseudolog10, yticks = [-100, -30, -10, -3, -1, 0])
        # ax3 = Makie.Axis(fig[3,1], xlabel = "Time (days)", ylabel = " ")#, title = "Title")
        # Makie.linkxaxes!(ax1, ax2, ax3); [Makie.hidexdecorations!(ax, grid = false) for ax in [ax1, ax2] ]
        scalars_jl = [:INTS_mm,:INTR_mm,:SNOW_mm]#, :GWAT]
        scalars_R  = [:ints,:intr,:snow]#, :gwat]
        scalars_label = ["INTS (mm)", "INTR (mm)", "SNOW (mm)"]#, "GWAT (mm)", "CC (MJ/m2)", "SNOWLQ (mm)"]
        gl = fig[3,1] = GridLayout()
        ax3_list = [Axis(gl[i, 1], ylabel = lbl) for (i, lbl) in enumerate(scalars_label)]
        Makie.linkxaxes!(ax1, ax2, ax3_list...);
        # [Makie.hidexdecorations!(ax, grid = false) for ax in [ax1, ax2, contents(gl)[Not(end)]...] ]
        [Makie.hidexdecorations!(ax, grid = false) for ax in [contents(gl)[Not(end)]...] ]

        Makie.series!(ax1, sim.θψδ.time, transpose(Matrix(sim.θψδ[:,r"^θ"])), color = Makie.wong_colors(), labels = "LWFBrook90.jl:" .* names(sim.θψδ[:,r"^θ"]))
        Makie.series!(ax2, sim.θψδ.time, transpose(Matrix(sim.θψδ[:,r"^ψ"])), color = Makie.wong_colors(), labels = "LWFBrook90.jl:" .* names(sim.θψδ[:,r"^ψ"]))
        # Makie.series!(ax3, sim.above.time, transpose(Matrix(sim.above[:,scalars_jl])),
        #             color = Makie.wong_colors(), labels = scalars_label)
        [Makie.lines!(ax3_list[i], sim.above.time, sim.above[:,col],
            color = Makie.wong_colors()[4]) for (i, col) in enumerate(scalars_jl)]

        # add LWFBrook90R solution
        Makie.series!(ax1, ref_NLAYER7.θ.time, transpose(Matrix(ref_NLAYER7.θ[:,Not(:time)])),
            linestyle = :dash, solid_color = :black, labels = "LWFBrook90R_NLayer7:" .* names(ref_NLAYER7.θ[:,Not(:time)]))
        Makie.series!(ax2, ref_NLAYER7.ψ.time, transpose(Matrix(ref_NLAYER7.ψ[:,Not(:time)])),
                    linestyle = :dash, solid_color = :black, labels = "LWFBrook90R_NLayer7:" .* names(ref_NLAYER7.ψ[:,Not(:time)]))
        # Makie.series!(ax3, ref_NLAYER7.above.time, transpose(Matrix(ref_NLAYER7.above[:,scalars_R])),
        #             linestyle = :dash, solid_color = :black, labels = ["LWFBrook90R_NLayer7_$i" for i in 1:4])
        [Makie.lines!(ax3_list[i], ref_NLAYER7.above.time, ref_NLAYER7.above[:,col], linestyle = :dash, color = :black) for (i, col) in enumerate(scalars_R)]
        # add subplot labels
        for (label, layout) in zip(["a", "b", "c"], [fig[1,1], fig[2,1], gl])
            Label(layout[1, 1, TopLeft()], label,
                fontsize = 26, font = :bold,
                padding = (0, 5, 5, 0),
                halign = :right)
        end

        # add automatic legends
        # fig[1, 2] = Legend(fig[1, 1], ax1, tellwidth = false)
        # fig[2, 2] = Legend(fig[2, 1], ax2, tellwidth = false)
        # # fig[3, 2] = Legend(fig[3, 1], ax3, tellwidth = false)
        # add manual legend

        # Legend(fig[1, 1, BottomLeft()], ## tellwidth = false, tellheight = true)
        axislegend(ax2,
            [LineElement(color    = Makie.wong_colors()[1], linestyle = nothing),
                LineElement(color = Makie.wong_colors()[2], linestyle = nothing),
                LineElement(color = Makie.wong_colors()[3], linestyle = nothing),
                LineElement(color = Makie.wong_colors()[4], linestyle = nothing),
                LineElement(color = :black, linestyle = :dash)],
            ["LWFBrook90.jl: ".*replace.(replace.(names(sim.θψδ[:,r"^θ"]), "θ_" => ""), "mm" => " mm")...,
                "LWFBrook90.jl: aboveground",
                "LWFBrook90R"],
            patchsize = (15, 5), rowgap = 0, labelsize = 10,
            position = :lb)

        save(fname_illustrations*"TESTSET_site-θ-ψ-aboveground-states_$(site).png", fig, pt_per_unit = 1)
        save(fname_illustrations*"TESTSET_site-θ-ψ-aboveground-states_$(site).pdf", fig, pt_per_unit = 1)
    end
end




# prepare plotting functions for Hammel (with/without δ)
function my_scatter!(pl, df; args...)
    Plots.scatter!(pl, float.(df[:,:time]), Matrix(df[:,Not(:time)]); args...)
end

function plot_Hammel_Dense(sim, ref, hyd, depth_to_read_out_mm, title; subtitle = "", fig = nothing, args...)
    # create title
    ODE_method = replace(String(Symbol(sim.ODEsolution.alg)), r"(.*?)(\(|{).*"=>s"\1")
    final_title = title
    final_subtitle = Printf.@sprintf("%d time steps with method: %s %s",
            sim.ODEsolution.destats.naccept,
            ODE_method,
            ifelse(subtitle=="","","\n"*subtitle))

    # Make publication figure with Makie # https://docs.makie.org/stable/explanations/figure_size/#vector_graphics
    # size_inches = (3.60, 4*3) # often either 7.25 inches or 3.60 inches wide
    if isnothing(fig)
        size_inches = (7.25, 12); # often either 7.25 inches or 3.60 inches wide
        size_pt = 72 .* size_inches;
        fig = Makie.Figure(resolution = size_pt, fontsize = 12);
    else
        # fig comes from outside
    end
    is_loam_simulation = maximum(sim.θψδ_dense.time)>20.
    @show is_loam_simulation
    if is_loam_simulation
        x_ticks = range(extrema(sim.θψδ_dense.time)..., step = 5);
        y_ticks_ψ = [-100, -30, -10, -3, -1, 0]
    else
        x_ticks = range(extrema(sim.θψδ_dense.time)..., step = 2);
        y_ticks_ψ = Makie.automatic; #[-100, -30, -10, -3, -1, 0]
        # y_ticks_ψ = round.(collect(10 .^ [-1:0.5:5...]))#[-100, -30, -10, -3, -1, 0]
        # y_ticks_ψ = round.(collect(10 .^ [-1:0.2:5...]))#
    end
    ax1 = Makie.Axis(fig[1,1], xlabel = "Time (days)", ylabel = "θ (-)",   xticks = x_ticks,
        title = final_title, subtitle = final_subtitle);
    ax2 = Makie.Axis(fig[2,1], xlabel = "Time (days)", ylabel = "ψ (kPa)", xticks = x_ticks,
        yscale = Makie.pseudolog10, yticks = y_ticks_ψ);
    simulate_isotopes = !ismissing(sim.θψδ_dense.δ18O_100mm[1])
    if (simulate_isotopes)
        ax3 = Makie.Axis(fig[3,1], xlabel = "Time (days)", ylabel = "δ¹⁸O (‰)", xticks = x_ticks)
        ax4 = Makie.Axis(fig[4,1], xlabel = "Time (days)", ylabel = "δ²H (‰)",  xticks = x_ticks)
        Makie.linkxaxes!(ax1, ax2, ax3, ax4); [Makie.hidexdecorations!(ax, grid = false) for ax in [ax1, ax2, ax3] ] # Not ax4
    else
        Makie.linkxaxes!(ax1, ax2); [Makie.hidexdecorations!(ax, grid = false) for ax in [ax1] ] # Not ax2
    end

    #LWFBrook90.jl:
    Makie.series!(ax1, sim.θψδ_dense.time, transpose(Matrix(sim.θψδ_dense[:,r"^θ"])), color = Makie.wong_colors())
    Makie.series!(ax2, sim.θψδ_dense.time, transpose(Matrix(sim.θψδ_dense[:,r"^ψ"])), color = Makie.wong_colors())
    if (simulate_isotopes)
        Makie.series!(ax3, sim.θψδ_dense.time, transpose(Matrix(sim.θψδ_dense[:,r"δ18O_"])), color = Makie.wong_colors())
        Makie.series!(ax4, sim.θψδ_dense.time, transpose(Matrix(sim.θψδ_dense[:,r"δ2H_"])),  color = Makie.wong_colors())
    end

    #Hydrus-1D
    Makie.series!(ax1, hyd.θdense.time, transpose(Matrix(disallowmissing(hyd.θdense[:,Not(:time)]))), color = Makie.wong_colors(), linewidth=0, marker = :circle, markersize = 8)
    # [Makie.scatter!(ax1, hyd.θdense.time, c, color = Makie.wong_colors()) for c in eachcol(Matrix(hyd.θdense[:,Not(:time)]))]
    Makie.series!(ax2, hyd.ψdense.time, transpose(Matrix(disallowmissing(hyd.ψdense[:,Not(:time)]))), color = Makie.wong_colors(), linewidth=0, marker = :circle, markersize = 8)
    if (simulate_isotopes)
        Makie.series!(ax3, hyd.δ18Odense.time, transpose(Matrix(disallowmissing(hyd.δ18Odense[:,Not(:time)]))), color = Makie.wong_colors(), linewidth=0, marker = :circle, markersize = 8)
        Makie.series!(ax4, hyd.δ2Hdense.time, transpose(Matrix(disallowmissing(hyd.δ2Hdense[:,Not(:time)]))), color = Makie.wong_colors(), linewidth=0, marker = :circle, markersize = 8)

    end

    #LWFBrook90R
    Makie.series!(ax1, ref.θ.time, transpose(Matrix(disallowmissing(ref.θ[:,Not(:time)]))), solid_color = :black, linewidth=0, marker = 'x', markersize = 8)
    Makie.series!(ax2, ref.ψ.time, transpose(Matrix(disallowmissing(ref.ψ[:,Not(:time)]))), solid_color = :black, linewidth=0, marker = 'x', markersize = 8)

    # add manual legend
    # axislegend(ax2,
    #     [LineElement(color    = Makie.wong_colors()[1], linestyle = nothing),
    #         LineElement(color = Makie.wong_colors()[2], linestyle = nothing),
    #         LineElement(color = Makie.wong_colors()[3], linestyle = nothing),
    #         LineElement(color = Makie.wong_colors()[4], linestyle = nothing),
    #         LineElement(color = Makie.wong_colors()[5], linestyle = nothing),
    #         MarkerElement(color = :black, marker = :circle), # markersize = 15, strokecolor = :black,
    #         MarkerElement(color = :black, marker = 'x'), # markersize = 15, strokecolor = :black, LineElement(color = :black, linestyle = :dash)
    #         ],
    #     ["LWFBrook90.jl: ".*replace.(replace.(names(sim.θψδ_dense[:,r"^θ"]), "θ_" => ""), "mm" => " mm")...,
    #         "Hydrus-1D", "LWFBrook90R"],
    #     patchsize = (15, 5), rowgap = 0, labelsize = 10,
    #     position = :rt)
    markers = [:hline, :circle, 'x']
    markers_string = ["LWFBrook90.jl","Hydrus-1D","LWFBrook90R"]
    # group_depth = [MarkerElement(marker = mar, color = :black, strokecolor = :transparent) for mar in markers]
    group_depth = [
        LineElement(color = :black, linestyle = nothing),
        MarkerElement(marker = markers[2], color = :black, strokecolor = :transparent),
        MarkerElement(marker = markers[3], color = :black, strokecolor = :transparent)]
    depth_string = replace.(replace.(names(sim.θψδ_dense[:,r"^θ"]), "θ_" => ""), "mm" => " mm")
    group_color = [PolyElement(color = color, strokecolor = :transparent) for color in Makie.wong_colors()[eachindex(depth_string)]]
    axislegend(ax1, [group_depth..., group_color...], [markers_string..., depth_string...],
        patchsize = (15, 5), rowgap = 0, labelsize = 10,
        position = :rt)

    return fig
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
    if plot_flag # if we plot out, compute again as first 3 runs were affected by compilation
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
    end
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
        using Makie
        depth_to_read_out_mm = [100, 500, 1000, 1500, 1900]
        # fname_illustrations = "out/$(today())/"

        fname_illustrations = "out/$git_status_string/"
        mkpath(dirname(fname_illustrations))

        # Plot Loam simulations
        pl1 = plot_Hammel_Dense(sim1, ref1, hyd1, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time1, gctime1, gcstats1)*"\n"*git_status_string)
        save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Loam_sim1.png", pl1, pt_per_unit = 1)
        # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Loam_sim1.pdf", pl1, pt_per_unit = 1)
        pl2 = plot_Hammel_Dense(sim2, ref2, hyd2, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time2, gctime2, gcstats2)*"\n"*git_status_string)
        save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Loam_sim2.png", pl2, pt_per_unit = 1)
        # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Loam_sim2.pdf", pl2, pt_per_unit = 1)
        if high_resolution_flag
            pl3 = plot_Hammel_Dense(sim3, ref3, hyd3, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900), dpi=300,
                subtitle = print_timed_statistics(time3, gctime3, gcstats3)*"\n"*git_status_string)
            save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Loam_sim3.png", pl3, pt_per_unit = 1)
            # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Loam_sim3.pdf", pl3, pt_per_unit = 1)
        end
        # Plot Sand simulations
        pl4 = plot_Hammel_Dense(sim4, ref4, hyd4, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time4, gctime4, gcstats4)*"\n"*git_status_string)
        save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Sand_sim1.png", pl4, pt_per_unit = 1)
        # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Sand_sim1.pdf", pl4, pt_per_unit = 1)
        pl5 = plot_Hammel_Dense(sim5, ref5, hyd5, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time5, gctime5, gcstats5)*"\n"*git_status_string)
        save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Sand_sim2.png", pl5, pt_per_unit = 1)
        # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Sand_sim2.pdf", pl5, pt_per_unit = 1)
        if high_resolution_flag
            pl6 = plot_Hammel_Dense(sim6, ref6, hyd6, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900), dpi=300,
                subtitle = print_timed_statistics(time6, gctime6, gcstats6)*"\n"*git_status_string)
            save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Sand_sim3.png", pl6, pt_per_unit = 1)
            # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Sand_sim3.pdf", pl6, pt_per_unit = 1)
        end

        # Also make a combined plot (for article):
        if high_resolution_flag
            size_inches = (7.25, 12); # often either 7.25 inches or 3.60 inches wide
            size_pt = 72 .* size_inches;
            fig_combined = Makie.Figure(resolution = size_pt, fontsize = 12);
            fig_a = fig_combined[1,1]; fig_b = fig_combined[1,2]
            plot_Hammel_Dense(sim3, ref3, hyd3, depth_to_read_out_mm, "Loam"; fig = fig_a)
            plot_Hammel_Dense(sim6, ref6, hyd6, depth_to_read_out_mm, "Sand"; fig = fig_b)
            # add maintitle
                # remove subtitles
                content(fig_a[1,1]).subtitle = ""
                content(fig_b[1,1]).subtitle = ""
                # add maintitle (super title)
                Label(fig_combined[0,:], "Simulation from Hammel et al. (2001)";
                    tellheight = true, tellwidth = false,
                    valign = :bottom, padding = (0, 0, 0, 0),
                    font = "TeX Gyre Heros Bold")
            # add subplot labels
            for (label, layout) in zip(["a", "b"], [fig_a, fig_b])
                Label(layout[1, 1, TopLeft()], label,
                    fontsize = 26, font = :bold,
                    padding = (0, 5, 5, 0),
                    halign = :right)
            end
            save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Sand3-Loam3.png", fig_combined, pt_per_unit = 1)
            save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ_Sand3-Loam3.pdf", fig_combined, pt_per_unit = 1)
        end
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
    if plot_flag # if we plot out, compute again as first 3 runs were affected by compilation
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
    end

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
        using Makie
        depth_to_read_out_mm = [100, 500, 1000, 1500, 1900]
        # fname_illustrations = "out/$(today())/"

        fname_illustrations = "out/$git_status_string/"
        mkpath(dirname(fname_illustrations))

        # Plot Loam simulations
        pl1 = plot_Hammel_Dense(sim1, ref1, hyd1, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time1, gctime1, gcstats1)*"\n"*git_status_string)
        save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Loam_sim1.png", pl1, pt_per_unit = 1)
        # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Loam_sim1.pdf", pl1, pt_per_unit = 1)

        pl2 = plot_Hammel_Dense(sim2, ref2, hyd2, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time2, gctime2, gcstats2)*"\n"*git_status_string)
        save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Loam_sim2.png", pl2, pt_per_unit = 1)
        # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Loam_sim2.pdf", pl2, pt_per_unit = 1)
        if high_resolution_flag
            pl3 = plot_Hammel_Dense(sim3, ref3, hyd3, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Loam"; size=(900,900), dpi=300,
                subtitle = print_timed_statistics(time3, gctime3, gcstats3)*"\n"*git_status_string)
            save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Loam_sim3.png", pl3, pt_per_unit = 1)
            # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Loam_sim3.pdf", pl3, pt_per_unit = 1)
        end
        # Plot Sand simulations
        pl4 = plot_Hammel_Dense(sim4, ref4, hyd4, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time4, gctime4, gcstats4)*"\n"*git_status_string)
        save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Sand_sim1.png", pl4, pt_per_unit = 1)
        # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Sand_sim1.pdf", pl4, pt_per_unit = 1)
        pl5 = plot_Hammel_Dense(sim5, ref5, hyd5, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900), dpi=300,
            subtitle = print_timed_statistics(time5, gctime5, gcstats5)*"\n"*git_status_string)
        save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Sand_sim2.png", pl5, pt_per_unit = 1)
        # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Sand_sim2.pdf", pl5, pt_per_unit = 1)
        if high_resolution_flag
            pl6 = plot_Hammel_Dense(sim6, ref6, hyd6, depth_to_read_out_mm, "Simulation from Hammel et al. (2001) - Sand"; size=(900,900), dpi=300,
                subtitle = print_timed_statistics(time6, gctime6, gcstats6)*"\n"*git_status_string)
            save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Sand_sim3.png", pl6, pt_per_unit = 1)
            # save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Sand_sim3.pdf", pl6, pt_per_unit = 1)
        end

        # Also make a combined plot (for article):
        if high_resolution_flag
            size_inches = (7.25, 12); # often either 7.25 inches or 3.60 inches wide
            size_pt = 72 .* size_inches;
            fig_combined = Makie.Figure(resolution = size_pt, fontsize = 12);
            fig_a = fig_combined[1,1]; fig_b = fig_combined[1,2]
            plot_Hammel_Dense(sim3, ref3, hyd3, depth_to_read_out_mm, "Loam"; fig = fig_a)
            plot_Hammel_Dense(sim6, ref6, hyd6, depth_to_read_out_mm, "Sand"; fig = fig_b)
            # add maintitle
                # remove subtitles
                content(fig_a[1,1]).subtitle = ""
                content(fig_b[1,1]).subtitle = ""
                # add maintitle (super title)
                Label(fig_combined[0,:], "Simulation from Hammel et al. (2001)";
                    tellheight = true, tellwidth = false,
                    valign = :bottom, padding = (0, 0, 0, 0),
                    font = "TeX Gyre Heros Bold")
            # add subplot labels
            for (label, layout) in zip(["a", "b"], [fig_a, fig_b])
                Label(layout[1, 1, TopLeft()], label,
                    fontsize = 26, font = :bold,
                    padding = (0, 5, 5, 0),
                    halign = :right)
            end
            save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Sand3-Loam3.png", fig_combined, pt_per_unit = 1)
            save(fname_illustrations*"TESTSET_Hammel-2001-θ-ψ-δ_Sand3-Loam3.pdf", fig_combined, pt_per_unit = 1)
        end
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
