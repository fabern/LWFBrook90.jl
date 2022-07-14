using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5, init#, step!
# example = LWFBrook90.run_example()
using Plots, Measures

# Read in input data
# input_prefix = "isoBEA2016-reset-FALSE"
# input_path = "examples/isoBEA2016-reset-FALSE-input/"
# input_prefix = "isoBEA2010-18-reset-FALSE";
# input_path = "examples/isoBEA2010-18-reset-FALSE-input/";
input_prefix = "isoBEAdense2010-18-reset-FALSE";
input_path = "examples/isoBEAdense2010-18-reset-FALSE-input/";
# input_prefix = "BEA2010-2021";
# input_path = "../../../LWF-Brook90.jl-calibration/Meteo-Data/BEA2010-2021/";
# input_prefix = "LAU2010-2021";
# input_path = "../../../LWF-Brook90.jl-calibration/Meteo-Data/LAU2010-2021/";
# input_prefix = "LAE2010-2021";
# input_path = "../../../LWF-Brook90.jl-calibration/Meteo-Data/LAE2010-2021/";
# input_prefix = "DAV2010-2021";
# input_path = "../../../LWF-Brook90.jl-calibration/Meteo-Data/DAV2010-2021/";
# input_prefix = "WaldLab";
# input_path = "../../../LWF-Brook90.jl-calibration/Meteo-Data/WaldLab/";
# input_prefix = "DAV_LW1_def";
# input_path = "../../../LWF-Brook90.jl-calibration/Meteo-Data/DAV_LW1_def/";

####################
simulate_isotopes = true
(input_meteoveg,
    input_meteoiso,
    input_meteoveg_reference_date,
    input_param,
    input_storm_durations,
    input_initial_conditions,
    input_soil_horizons,
    simOption_FLAG_MualVanGen) = read_inputData(input_path, input_prefix;
                                                simulate_isotopes = simulate_isotopes);
# Define grid for spatial discretization as well as initial conditions and root densities
# a) either read the discretization from a file `soil_discretization.csv`
unused = discretize_soil(input_path, input_prefix)

# b) or define them manually
ψ_init = unused.uAux_PSIM_init_kPa[1]
# Δz_m = round.(-diff([unused.Upper_m[1]; unused.Lower_m]), digits=5)
# Δz_m = [fill(0.005, 4); fill(0.02, 99)]                          # grid spacing (heterogenous), meter (N=103)
# Δz_m = round.(diff(0.0:0.05:-minimum(unused.Lower_m)), digits=5) # grid spacing (heterogenous), meter

# As a test: subsequently increase resolution of the top layers.
Δz_m = [0.04, 0.04, 0.12, 0.25, 0.3, 0.35, 0.1];                           # grid spacing (heterogenous), meter (N=7)
# Δz_m = [0.04, 0.04, 0.04, 0.08, 0.25, 0.3, 0.35, 0.1]                    # grid spacing (heterogenous), meter (N=8)
# Δz_m = [0.04, 0.04, 0.04, 0.04, 0.04, 0.25, 0.3, 0.35, 0.1]              # grid spacing (heterogenous), meter (N=9)
# Δz_m = [fill(0.04, 5); fill(0.05, 5); 0.3; 0.35; 0.1]                    # grid spacing (heterogenous), meter (N=13)
# Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); 0.35; 0.1]          # grid spacing (heterogenous), meter (N=17)
Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 5); 0.1]; # grid spacing (heterogenous), meter (N=21)
# Δz_m = [fill(0.02, 10); fill(0.025, 10); fill(0.03, 10); fill(0.035, 10); 0.1]; # grid spacing (heterogenous), meter (N=41)
# Δz_m = [fill(0.02, 60)... ]; # grid spacing (homgenous), meter (N=60)
# Δz_m = [fill(0.01, 20); fill(0.0125, 20); fill(0.015, 20); fill(0.0175, 20); 0.1]; # grid spacing (heterogenous), meter (N=81)
# Δz_m = [0.04, 0.04, 0.12, 0.25, 0.3, 0.35, 0.1];                           # grid spacing (heterogenous), meter (N=7)
# Δz_m = [1.20];# grid spacing (heterogenous), meter (N=1)


for Δz_m in (
    [0.04, 0.04, 0.12, 0.25, 0.3, 0.35, 0.1],
    [fill(0.04, 5); fill(0.05, 5); 0.3;            0.35;         0.1], # N=13
    [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 5); 0.1], # N=21
    [fill(0.02, 10); fill(0.025, 10); fill(0.03, 10); fill(0.035, 10); 0.1], # N=41
    [fill(0.02, 60)... ], # N=60, 2cm similar to Pollacco et al. 2022 suggestions
    )

    f1 = (Δz_m) -> LWFBrook90.Rootden_beta_(0.97, Δz_m = Δz_m)  # function for root density as f(Δz)
    f2 = (Δz_m) -> fill(-6.3, length(Δz_m))          # function for initial conditions as f(Δz)
    input_soil_discretization = discretize_soil(;
        Δz_m = Δz_m,
        Rootden_ = f1,
        uAux_PSIM_init_kPa = f2,
        u_delta18O_init_permil = ifelse.(cumsum(Δz_m) .<= 0.2, -13., -10.),
        u_delta2H_init_permil  = ifelse.(cumsum(Δz_m) .<= 0.2, -95., -70.))
    # input_soil_discretization = unused
    ####################

    ####################
    # Define solver options
    Reset = false;                          # currently only Reset = 0 implemented
    compute_intermediate_quantities = true; # Flag whether ODE containes additional quantities than only states

    # Override input file settings
    # Here possibility to check and override dataframes input_[...] manually
        # # E.g:
        # # Soil hydraulic model
        # input_param[1,"NOOUTF"] = true # `true` if outflow from roots prevented, `false` if allowed
        # input_soil_horizons.gravel_volFrac .= 0.00
    ####################

    ####################
    # Define parameters for differential equation
    (ψM_initial,δ18O_initial,δ2H_initial), p = define_LWFB90_p(
        input_meteoveg,
        input_meteoiso,
        input_meteoveg_reference_date,
        input_param,
        input_storm_durations,
        input_soil_horizons,
        input_soil_discretization,
        simOption_FLAG_MualVanGen;
        Reset = Reset,
        compute_intermediate_quantities = compute_intermediate_quantities,
        simulate_isotopes = simulate_isotopes,
        # soil_output_depths = [-0.35, -0.42, -0.48, -1.05]
        # soil_output_depths = collect(-0.05:-0.05:-1.1)
        );

    # using Plots
    # hline([0; cumsum(p[1][1].p_THICK)], yflip = true, xticks = false,
    #     title = "N_layer = "*string(p[1][1].NLAYER))
    ####################

    ####################
    # Define initial states of differential equation
    # state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI
    # Create u0 for DiffEq.jl
    u0, p = define_LWFB90_u0(p, input_initial_conditions,
        ψM_initial, δ18O_initial, δ2H_initial,
        compute_intermediate_quantities;
        simulate_isotopes = simulate_isotopes);
    ####################

    ####################
    # Define ODE problem which consists of
    #   - definition of right-hand-side (RHS) function f
    #   - definition of callback function cb
    #   - u0:     initial condition of states
    #   - tspan:  definition of simulation time span
    #   - p:      parameters

    # Define simulation time span:
    # tspan = (0.,  5.) # simulate 5 days
    # tspan = (0.,  100.) # simulate 100 days # NOTE: KAU bugs in "branch 005-" when at least 3*365
    tspan = (minimum(input_meteoveg[:,"days"]),
            maximum(input_meteoveg[:,"days"])) # simulate all available days
    # tspan = (LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1980,1,1), reference_date),
    #          LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1985,1,1), reference_date)) # simulates selected period
    # tspan = (0., 700.)
    # tspan = (0.,  60)
    # tspan = (0.,  3)
    ####################

    ####################
    ## Solve ODE:
    # sol_LWFBrook90 = solve_LWFB90(u0, tspan, p);
    ####################

    ####################
    ## Benchmarking of input_path = "examples/isoBEAdense2010-18-reset-FALSE-input/"
    @time sol_LWFBrook90 = solve_LWFB90(u0, tspan, p);
        # 700 days in: 40 seconds dt=1e-3, adaptive = false (isoBea default spacing: NLAYER = 7)
        # 700 days in: 90 seconds dt=1e-3, adaptive = false (isoBea dense spacing (0.05m): NLAYER = 26)
        #  git+c4d37eb+gitdirty: NLAYER=7:  25.944645 seconds (196.27 M allocations: 16.947 GiB, 15.09% gc time)
        #  git+c4d37eb+gitdirty: NLAYER=13: 34.060465 seconds (213.98 M allocations: 24.347 GiB, 17.25% gc time)
        #  git+c4d37eb+gitdirty: NLAYER=21: 47.000202 seconds (237.66 M allocations: 34.254 GiB, 16.13% gc time)
        #  git+61a19ed+gitclean: NLAYER=7:  22.871329 seconds (175.84 M allocations: 16.795 GiB, 16.08% gc time)
        #  git+61a19ed+gitclean: NLAYER=13: 31.875148 seconds (176.05 M allocations: 24.064 GiB, 17.32% gc time)
        #  git+61a19ed+gitclean: NLAYER=21: 45.809290 seconds (176.38 M allocations: 33.797 GiB, 16.53% gc time)
        #  git+a3fc1a2+gitclean: NLAYER=7:  31.444282 seconds (289.68 M allocations: 18.872 GiB, 14.39% gc time, 27.43% compilation time)
        #  git+a3fc1a2+gitclean: NLAYER=13: 37.432710 seconds (379.27 M allocations: 27.501 GiB, 17.63% gc time)
        #  git+a3fc1a2+gitclean: NLAYER=21: 53.301924 seconds (527.72 M allocations: 40.511 GiB, 16.98% gc time)
        #  git+0fc6e1b+gitclean: NLAYER=7:  20.076741 seconds (164.78 M allocations: 15.548 GiB, 16.56% gc time)
        #  git+0fc6e1b+gitclean: NLAYER=13: 27.925882 seconds (164.96 M allocations: 22.139 GiB, 16.92% gc time)
        #  git+0fc6e1b+gitclean: NLAYER=21: 40.381837 seconds (165.29 M allocations: 30.971 GiB, 16.54% gc time)
        #  git+d7a34b5+gitclean: NLAYER=7: 206.039105 seconds (1.66 G allocations: 155.787 GiB, 16.17% gc time, 4.00% compilation time)
        #  git+d7a34b5+gitclean: NLAYER=13: 285.350438 seconds (1.65 G allocations: 221.057 GiB, 17.95% gc time)
        #  git+d7a34b5+gitclean: NLAYER=21: 401.913824 seconds (1.65 G allocations: 309.421 GiB, 17.12% gc time)
        #  git+25e3bd0+gitclean: NLAYER=7 : 34.753962 seconds (289.67 M allocations: 18.871 GiB, 14.32% gc time, 27.87% compilation time)
        #  git+25e3bd0+gitclean: NLAYER=13: 35.102646 seconds (379.27 M allocations: 27.501 GiB, 18.29% gc time)
        #  git+25e3bd0+gitclean: NLAYER=21: 52.302226 seconds (527.72 M allocations: 40.511 GiB, 17.46% gc time)
        #  git+9fe5a08+gitclean: NLAYER= 7: 0.583609 seconds (5.93 M allocations: 453.802 MiB, 14.74% gc time)
        #  git+9fe5a08+gitclean: NLAYER=13: 0.947032 seconds (8.42 M allocations: 730.844 MiB, 18.08% gc time)
        #  git+9fe5a08+gitclean: NLAYER=21: 1.828127 seconds (7.66 M allocations: 726.068 MiB, 50.04% gc time)
        #  git+d1eb2f8+gitclean: NLAYER=21: 1.126634 seconds (4.22 M allocations: 683.778 MiB, 14.38% gc time)
        #  git+eea6abd+gitdirty: NLAYER=21: 4.136514 seconds (16.66 M allocations: 2.509 GiB, 13.65% gc time)
        #  git+40de77b+gitdirty: NLAYER=21: 4.030769 seconds (14.84 M allocations: 2.162 GiB, 12.15% gc time)
        #  git+3851514+gitdirty-iso+false NLAYER= 7: 0.540389 seconds (2.63 M allocations: 292.582 MiB, 13.41% gc time, 0.42% compilation time)
        #  git+3851514+gitdirty-iso+false NLAYER=13: 0.677501 seconds (2.78 M allocations: 438.114 MiB, 22.44% gc time)
        #  git+3851514+gitdirty-iso+false NLAYER=21: 4.009734 seconds (10.14 M allocations: 2.174 GiB, 14.54% gc time)
        #  git+3851514+gitdirty-iso+true  NLAYER= 7:  5.625118 seconds (26.38 M allocations: 2.188 GiB, 10.43% gc time)
        #  git+3851514+gitdirty-iso+true  NLAYER=13:  5.301599 seconds (26.90 M allocations: 2.967 GiB, 13.61% gc time)
        #  git+3851514+gitdirty-iso+true  NLAYER=21: 15.839506 seconds (56.62 M allocations: 8.301 GiB, 12.66% gc time)
        #  git+3851514+gitdirty-iso+true  NLAYER=41: 35.529361 seconds (84.29 M allocations: 20.755 GiB, 14.10% gc time)
        #  git+3851514+gitdirty-iso+true  NLAYER=60: XXXXX seconds
        #  git+e6ad7ad+gitdirty-iso+true  NLAYER= 7: 0.231580 seconds (1.07 M allocations: 104.647 MiB, 29.50% gc time) 2203 Δt's
        #  git+e6ad7ad+gitdirty-iso+true  NLAYER=13: 0.312734 seconds (1.26 M allocations: 165.862 MiB) 2555 Δt's
        #  git+e6ad7ad+gitdirty-iso+true  NLAYER=21: 0.857902 seconds (3.19 M allocations: 548.413 MiB, 10.90% gc time) 7020 Δt's
        #  git+e6ad7ad+gitdirty-iso+true  NLAYER=41: 2.700105 seconds (6.28 M allocations: 1.763 GiB, 15.09% gc time) 14469 Δt's
        #  git+e6ad7ad+gitdirty-iso+true  NLAYER=60: 3.271614 seconds (4.72 M allocations: 1.875 GiB, 13.70% gc time) 10619 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+7-git+544f265+gitclean-iso+true:  1.100431 seconds (6.02 M allocations: 582.128 MiB, 12.61% gc time) 12417 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+13-git+544f265+gitclean-iso+true: 1.689267 seconds (6.84 M allocations: 893.392 MiB, 15.44% gc time) 14124 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+21-git+544f265+gitclean-iso+true: 1.261623 seconds (4.01 M allocations: 703.352 MiB, 14.22% gc time) 8411 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+14-git+544f265+gitclean-iso+true: 3.212237 seconds (5.96 M allocations: 1.726 GiB, 33.60% gc time) 12811 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+62-git+544f265+gitclean-iso+true: 1.224711 seconds (1.96 M allocations: 823.858 MiB, 16.57% gc time) 4089 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+7-git+544f265+gitclean-iso+true gravelFrac=0:  0.249753 seconds (1.07 M allocations: 104.647 MiB, 23.68% gc time) 2203 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+13-git+544f265+gitclean-iso+true gravelFrac=0: 0.307845 seconds (1.26 M allocations: 165.862 MiB, 22.77% gc time) 2555 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+21-git+544f265+gitclean-iso+true gravelFrac=0: 1.090218 seconds (3.19 M allocations: 548.413 MiB, 12.95% gc time) 7020 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+14-git+544f265+gitclean-iso+true gravelFrac=0: 2.708308 seconds (6.28 M allocations: 1.763 GiB, 16.18% gc time) 14469 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+62-git+544f265+gitclean-iso+true gravelFrac=0: 2.620070 seconds (4.72 M allocations: 1.875 GiB, 12.91% gc time) 10619 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+7-git+4cdb5e8+gitclean-iso+true gravelFrac=0: 0.156029 seconds (1.11 M allocations: 105.185 MiB) 2203 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+13-git+4cdb5e8+gitclean-iso+true gravelFrac=0: 0.536214 seconds (1.30 M allocations: 166.486 MiB, 50.42% gc time) 2555 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+21-git+4cdb5e8+gitclean-iso+true gravelFrac=0: 0.994925 seconds (3.30 M allocations: 550.127 MiB, 16.16% gc time) 7020 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+14-git+4cdb5e8+gitclean-iso+true gravelFrac=0: 3.067565 seconds (6.52 M allocations: 1.767 GiB, 14.51% gc time) 14469 Δt's
        # isoBEAdense2010-18-reset-FALSE_NLAYER+62-git+4cdb5e8+gitclean-iso+true gravelFrac=0: 2.926338 seconds (4.89 M allocations: 1.878 GiB, 13.04% gc time) 10619 Δt's
                # TODO(bernhard): there is a bug in NLAYER=60...

                    # integrator =  init(ode_LWFBrook90, Tsit5(); progress = true,
                    #             unstable_check = unstable_check_function) # = (dt,u,p,t) -> false, #any(isnan,u))
                    # du = copy(integrator.u)
                    # du .= 0
                    # LWFBrook90.f_LWFBrook90R(du, integrator.u, p, 0);
                    # du
                    # step!(integrator)

    # using BenchmarkTools # for benchmarking
    # sol_LWFBrook90 = @btime solve(ode_LWFBrook90; dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity
    # sol_LWFBrook90 = @btime solve(ode_LWFBrook90; saveat = tspan[1]:tspan[2], dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity
    ####################

    ####################
    ## Plotting

    θ_limits = (minimum(sol_LWFBrook90.prob.p[1][1].p_θr), maximum(sol_LWFBrook90.prob.p[1][1].p_THSAT))
    ψ_pF_limits = (-2, 7)

    # 0) defined in module_ISO.jl
    # using Plots, Measures
    mkpath("out")
    using DataFrames
    using CSV: File
    using Plots, Measures
    using Dates: DateTime, format
    fname = joinpath(
        "out",
        input_prefix*"_NLAYER+" * string(sol_LWFBrook90.prob.p[1][1].NLAYER)*
        "-git+"*chomp(Base.read(`git rev-parse --short HEAD`, String))*
        ifelse(length(Base.read(`git status --porcelain`, String))==0, "+gitclean","+gitdirty")*
        "-iso+"*string(simulate_isotopes))

    if simulate_isotopes
        optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
        # pl1 = LWFBrook90.ISO.plotisotopes(sol_LWFBrook90);
        pl2 = LWFBrook90.ISO.plotisotopes(
            sol_LWFBrook90, optim_ticks;
            layout = grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4]),
            size=(1000,1400), dpi = 300, leftmargin = 15mm);
        plot!(pl2, link = :x);
        savefig(pl2, fname*"_plotRecipe.png")
    end
    if input_prefix == "WaldLab"
        pl2b = LWFBrook90.ISO.plotisotopes(
                sol_LWFBrook90, optim_ticks;
                layout = grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4]),
                size=(1000,1400), dpi = 300, leftmargin = 15mm,
                xlim = (DateTime("2020-01-01"), DateTime("2020-12-31")))
        savefig(pl2b, fname*"_plotRecipe-2020.png")
    end
    # plot RWU time series
    # idx_u_scalar_isotopes_d18O = sol_LWFBrook90.prob.p[1][4][8     ]
    #     idx_u_vector_isotopes_d18O = sol_LWFBrook90.prob.p[1][4][9     ]
    #     idx_u_scalar_isotopes_d2H  = sol_LWFBrook90.prob.p[1][4][10     ]
    #     idx_u_vector_isotopes_d2H  = sol_LWFBrook90.prob.p[1][4][11     ]
    if simulate_isotopes
        row_RWU_d18O  = reshape(sol_LWFBrook90[p[1][4].row_idx_scalars.totalRWU, sol_LWFBrook90.prob.p[1][4].col_idx_d18O, :, 1], 1, :)
        row_XYL_d18O  = reshape(sol_LWFBrook90[p[1][4].row_idx_scalars.XylemV,   sol_LWFBrook90.prob.p[1][4].col_idx_d18O, :, 1], 1, :)
        plot([transpose(row_RWU_d18O) transpose(row_XYL_d18O)], labels = ["δ_RWU" "δ_XylemV"])
    end


    if input_prefix == "BEA2016-reset-FALSE" || input_prefix == "isoBEA2016-reset-FALSE" || input_prefix == "isoBEAdense2016-reset-FALSE"
        ref_aboveground =
            DataFrame(File(
            # "test/test-assets-external/BEA-2016/BEA-IntegrationTests-LWFBrook90/output_LWFBrook90R/BEA2016-reset-FALSE_NLAYER14_LWFBrook90R-0.4.5daily_output.csv"))[:,
            "out/BEA2016-reset-FALSE_NLAYER14_LWFBrook90R-0.4.5daily_output.csv"))[:,
                [:yr, :mo, :da, :doy, :intr, :ints, :snow, :gwat]]
        pl_ab_3 = plot(sol_LWFBrook90; vars = [2, 3, 4],
            label=["INTS (mm)" "INTR (mm)" "SNOW (mm)"])
        plot!(pl_ab_3,
                [ref_aboveground.intr,
                ref_aboveground.ints,
                ref_aboveground.snow], label = "LWFBrook90R", line = :dash, color = :black)
        savefig(fname*"_plot-INTS_INTR_SNOW.png")
        pl_ab_4 = plot(sol_LWFBrook90; vars = [2, 3],
            label=["INTS (mm)" "INTR (mm)"])
        plot!(pl_ab_4,
                [ref_aboveground.intr,
                ref_aboveground.ints], label = "LWFBrook90R", line = :dash, color = :black)
        savefig(fname*"_plot-INTS_INTR.png")
    end


    if (true)
        # theme(:default) # theme(:sand)
        PREC_color = :black
        depth_to_read_out_mm = [150 500 800 1500]
        if simulate_isotopes
            δ_resultsSoil = LWFBrook90.get_δsoil(depth_to_read_out_mm, sol_LWFBrook90)
            δ_results = get_δ(sol_LWFBrook90)
        end

        pl_θ = plot(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
            LWFBrook90.get_θ(depth_to_read_out_mm, sol_LWFBrook90),
            labels = string.(depth_to_read_out_mm) .* "mm",
            xlabel = "Date",
            ylabel = "θ\n[-]", ylims = θ_limits,
            legend = :outerright)
        pl_ψ = plot(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
            # -LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90) .+ 1, yaxis = :log, yflip = true,
            # LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90),  ylabel = "ψ\n[kPa]",
            log10.(-10 .* LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90)),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
            ylims = ψ_pF_limits,
            labels = string.(depth_to_read_out_mm) .* "mm",
            xlabel = "Date",
            legend = :outerright);
        if simulate_isotopes
            pl_δ18O = plot(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
                δ_resultsSoil[1],
                labels = string.(depth_to_read_out_mm) .* "mm",
                xlabel = "Date",
                ylabel = "δ¹⁸O soil\n[‰]",
                legend = :outerright);
            pl_δ2H = plot(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
                δ_resultsSoil[2],
                labels = string.(depth_to_read_out_mm) .* "mm",
                xlabel = "Date",
                ylabel = "δ²H soil\n[‰]",
                legend = :outerright);
            # add precipitation to soil δ
            plot!(pl_δ2H,
                LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
                δ_results.PREC.d2H', labels = "PREC", color = PREC_color, linestyle = :dot);
            plot!(pl_δ18O,
                LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
                δ_results.PREC.d18O', labels = "PREC", color = PREC_color, linestyle = :dot);
                # # add raw data values to check interpolation
                # data_days_to_plot = (input_meteoiso.days .> tspan[1] .&& input_meteoiso.days .< tspan[2])
                # scatter!(LWFBrook90.RelativeDaysFloat2DateTime.(input_meteoiso.days[data_days_to_plot], input_meteoveg_reference_date),
                #         input_meteoiso.delta18O_permil[data_days_to_plot])

            row_RWU_d18O  = reshape(sol_LWFBrook90[p[1][4].row_idx_scalars.totalRWU, sol_LWFBrook90.prob.p[1][4].col_idx_d18O, :, 1], 1, :)
            row_XYL_d18O  = reshape(sol_LWFBrook90[p[1][4].row_idx_scalars.XylemV,   sol_LWFBrook90.prob.p[1][4].col_idx_d18O, :, 1], 1, :)
            plot!(pl_δ18O, linestyle = :dash,
                 LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
                 [transpose(row_RWU_d18O) transpose(row_XYL_d18O)], labels = ["δ_RWU" "δ_XylemV"])
        else
            pl_δ18O = plot();
            pl_δ2H = plot();
        end
        pl_PREC = plot(
            LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
            sol_LWFBrook90.prob.p[2][8].(sol_LWFBrook90.t),
            t = :bar, color=PREC_color,
            legend = :outerright, labels = "PREC    ", # whitespace for hardcoded alignment of legend
            ylabel = "PREC\n[mm]");
        plot(plot(pl_PREC, xlab = "", xticks = :none, topmargin = 5mm, bottommargin = 0mm),
            plot(pl_θ;     xlab = "", xticks = :none, topmargin = 0mm, bottommargin = 0mm),
            plot(pl_ψ;     xlab = "", xticks = :none, topmargin = 0mm, bottommargin = 0mm),
            plot(pl_δ18O;  xlab = "", xticks = :none, topmargin = 0mm, bottommargin = 0mm),
            plot(pl_δ2H;   xtick_direction=:out     , topmargin = 0mm, bottommargin = 5mm),
            link = :x,
            layout = grid(5, 1, heights=[0.1, 0.25 ,0.25, 0.2, 0.2]),
            size=(600,500), dpi = 300, margin = 5mm);

        savefig(fname*"_plot-θ-ψ-δ.png")
    end

    aux_indices = sol_LWFBrook90.prob.p[1][4].row_idx_accum
    aux_names = sol_LWFBrook90.prob.p[1][4].names_accum
    plot(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
                [sol_LWFBrook90[aux_indices[30],:] sol_LWFBrook90[aux_indices[31],:]],
                legend = :outerright, labels = aux_names[:, 30:31],
                ylabel = "Water balance error [mm]")
    savefig(fname*"_plot-water-balance-error.png")

    if (true)
        # Depth-vs-time heatmap plots of θ, ψ, RWU
        # For all heatmaps
        t_ref = sol_LWFBrook90.prob.p[2][17]
        x = RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, t_ref);
        y = cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK) - sol_LWFBrook90.prob.p[1][1].p_THICK/2
        optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
        y_soil_ticks = optim_ticks(0., round(maximum(cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
        y_ticks= y_soil_ticks
        y_labels=round.(y_soil_ticks; digits=0)
        # a) for θ and ψ:
        (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
                        LWFBrook90.get_auxiliary_variables(sol_LWFBrook90)
        t_to_plot = 1:length(x)
        pl_θ = heatmap(x, y,
                    u_aux_θ', colorbar_title = "θ [-]", clims = θ_limits, c=cgrad(:inferno, rev=true),
                    yflip = true,
                    xlabel = "",#"Date",
                    ylabel = "Depth [mm]");
        pl_ψ = heatmap(x, y,
                    # u_aux_PSIM', colorbar_title = "ψ_m [kPa]",
                    log10.(-10 .* u_aux_PSIM'),  colorbar_title = "pF = log₁₀(-ψ hPa)", clims = ψ_pF_limits, c=cgrad(:inferno, rev=false),
                    yflip = true,
                    xlabel = "",#"Date",
                    ylabel = "Depth [mm]");
        # -LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90) .+ 1, yaxis = :log, yflip = true,
            # LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90),  ylabel = "ψ\n[kPa]",



        # plot(pl_θ,pl_ψ,layout = (2,1))
        # b) for RWU
        # b1) RWU heatmap (depths vs time)
        rows_RWU_mmDay = sol_LWFBrook90[sol_LWFBrook90.prob.p[1][4].row_idx_RWU, 1, :]
        rows_RWU = rows_RWU_mmDay ./ sol_LWFBrook90.prob.p[1][1].p_THICK
        pl_RWU = heatmap(x,y,rows_RWU; yticks = (y_ticks, y_labels),
                yflip = true,
                c = :diverging_bwr_20_95_c54_n256, clim = maximum(abs.(rows_RWU)) .* (-1, 1),
                ylabel = "Depth [mm]",
                colorbar_title = "RWU [mm water/day per mm soil depth]",
                size=(1400,800), dpi = 300, leftmargin = 15mm);
        # b2) distribution of RWU and roots
        plot_avgRWU = plot(sum(rows_RWU,dims=2), y,
            yflip = true, ylabel = "Depth [mm]", #xlabel = "Total RWU over\nsimulation period [mm]")
            #xlabel = "RWU (mm)"
            legend = :bottomright
        )
        t_toPlot = range(extrema(sol_LWFBrook90.prob.tspan)..., length=5)
        root_data_start_end = [sol_LWFBrook90.prob.p[2].p_RELDEN.(t, 1:sol_LWFBrook90.prob.p[1][1].NLAYER) for t in t_toPlot]
        pl_roots = plot(root_data_start_end, y,
            labels = (Dates.format.(permutedims(RelativeDaysFloat2DateTime.(t_toPlot, t_ref)), "yyyy-mm")),
            linestyle = [:solid :solid :solid :solid :dash],
            yflip = true, ylabel = "Depth [mm]", xlabel = "Relative root\ndensity [-]", legend = :bottomright);
        empty_plot = plot(legend=false,grid=false,foreground_color_subplot=:white)
        # pl_RWU2 = plot(pl_θ, pl_ψ, pl_RWU,
        #         empty_plot, plot(pl_roots, ylabel = ""), plot(plot_avgRWU, ylabel = ""),
        #         # layout = @layout [[a{0.333h};b;c{0.333h}] grid(3,1){0.2w}]
        #         title = permutedims(["a) θ (m3/m3):", "b) ψ (kPa):", "c) RWU (mm/day):",
        #                              " ",
        #                              "d) Root density distribution:",
        #                              "e) Total RWU (mm) averaged over time:"]),
        #         titleloc = :left, titlefont = font(8),
        #         size=(1400,800), dpi = 300, leftmargin = 15mm, bottommargin = 10mm,
        #         # layout = @layout [[a{0.333h};b;c{0.333h}] [f;d{0.333h};e{0.333h}] # does not work with {0.2w}
        #         layout = @layout [[a{0.333h};b;c{0.333h}] grid(3,1){0.2w}]
        # )
        pl_RWU2 = plot(
                pl_θ, empty_plot,
                pl_ψ, plot(pl_roots,xlabel=""),
                pl_RWU, plot(plot_avgRWU,xlabel=""),

                title = permutedims(["a) θ (m3/m3):",                              " ",
                                    "b) ψ (kPa):",                                 "d) Root density (-) distribution:",
                                    "c) RWU (mm water/day per mm soil depth):",    "e) Total RWU (mm water / mm soil)\naveraged over time:"]),
                titleloc = :left, titlefont = font(8),
                size=(1400,800), dpi = 300, leftmargin = 15mm, bottommargin = 10mm,
                link = :y,
                layout = @layout [a{0.333h} d{0.2w}; b e; c{0.333h} f]
        );
        savefig(pl_RWU2, fname*"_RWU2.png")
        # savefig(plot(pl_θ, pl_ψ, pl_RWU, layout = (3,1)),
        #         fname*"_RWU.png")
        # investigate what happens on 2010-05-22 in BEA (artifact in θ simulation)?
        # ode_LWFBrook90_2, unstable_check_function_2 = define_LWFB90_ODE(sol_LWFBrook90.u[130], (130,140), p);
        # @time sol_LWFBrook90_2 = solve(ode_LWFBrook90_2, Tsit5(); progress = true,
        #     unstable_check = unstable_check_function_2, # = (dt,u,p,t) -> false, #any(isnan,u),
        #     saveat = 130:140,
        #     # When adding transport of isotopes adaptive time stepping has difficulties reducing dt below dtmin.
        #     # We solve it either by using a constant time step
        #     # or by setting force_dtmin to true, which has the drawback that tolerances are ignored.
        #     # # Variant1: less sophisticated but robust:
        #     # adaptive = false, dt=1e-3 #dt=5/60/24 # fixed 5 minutes time steps
        #     # Variant2: more sophisticated but giving errors:
        #     adaptive = true, dtmin = 1e-5, dt=1e-4,  # if adaptive dt is just the starting value of the time steps
        #     force_dtmin = true,# without this callbacks generate an abort "dt <= dtmin",
        #                     # e.g. see: https://github.com/SciML/DifferentialEquations.jl/issues/648
        #     maxiters = (tspan[2]-tspan[1])/1e-4 # when using force_dtmin also use maxiters
        #                     # TODO(bernhard): regarding maxiters it seems to be the callback for δ of SNOW, INTR, INTS that causes the dtmin to be so small to reach large maxiters
        #     );
        # (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        #                 LWFBrook90.get_auxiliary_variables(sol_LWFBrook90_2)
        # heatmap(RelativeDaysFloat2DateTime.(sol_LWFBrook90_2.t, t_ref),
        #         cumsum(sol_LWFBrook90_2.prob.p[1][1].p_THICK) - sol_LWFBrook90_2.prob.p[1][1].p_THICK/2,
        #         u_aux_θ, yflip = true)
    end
end

# # # 1) very basic
# # using Plots # Install plot package at first use with `]` and then `add Plots`
# # # plot(p[2][15](1:300))
# # # plot(p[2][16](1:300))

# # # # # Plot 1
# # plot(sol_LWFBrook90; vars = [1, 2, 3, 4, 5, 6],
# #      label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

# # idx_u_vector_amounts       = sol_LWFBrook90.prob.p[1][4].row_idx_SWATI
# # idx_u_scalar_isotopes_d18O = sol_LWFBrook90.prob.p[1][4].row_idx_accum
# # idx_u_vector_isotopes_d18O = sol_LWFBrook90.prob.p[1][4][6  ]
# # idx_u_scalar_isotopes_d2H  = sol_LWFBrook90.prob.p[1][4].names_accum
# # idx_u_vector_isotopes_d2H  = sol_LWFBrook90.prob.p[1][4][8  ]
# # idx_u_vector_accumulators  = sol_LWFBrook90.prob.p[1][4][9  ]


# # a = plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d18O,
# #     ylabel = "δ18O [‰]",
# #     label=["GWAT (‰)" "INTS (‰)" "INTR (‰)" "SNOW (‰)"])
# # plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), label = "PREC (‰)")
# # # plot!(ylims = (-50,50))

# # plot(plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d18O[2], ylabel = "δ18O [‰]", label=["INTS (‰)"]),
# #     plot(sol_LWFBrook90; vars = [2],label=["INTS (mm)"]),
# #     layout=(2,1))
# # # plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), label = "PREC (‰)")
# # plot(plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d18O[3], ylabel = "δ18O [‰]", label=["INTR (‰)"]),
# #     plot(sol_LWFBrook90; vars = [3],label=["INTR (mm)"]),
# #     layout=(2,1))
# # # plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), label = "PREC (‰)")
# # plot(plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d18O[4], ylabel = "δ18O [‰]", label=["SNOW (‰)"]),
# #     plot(sol_LWFBrook90; vars = [4],label=["SNOW (mm)"]),#, ylims = (0,0.01)),
# #     layout=(2,1))
# # # plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), label = "PREC (‰)")
# # plot(plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d18O[1], ylabel = "δ18O [‰]", label=["GWAT (‰)"]),
# #     plot(sol_LWFBrook90; vars = [1],label=["GWAT (mm)"]),
# #     layout=(2,1))
# # # plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), label = "PREC (‰)")


# # # sol_LWFBrook90[26,1,1:10]
# # plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d2H,
# #     ylabel = "δ2H [‰]",
# #     label=["GWAT (‰)" "INTS (‰)" "INTR (‰)" "SNOW (‰)"])
# # # plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][16].(sol_LWFBrook90.t), label = "PREC (‰)")
# # scatter([
# #         (sol_LWFBrook90[13 + 1,:], sol_LWFBrook90[24 + 1,:]),
# #         (sol_LWFBrook90[13 + 2,:], sol_LWFBrook90[24 + 2,:]),
# #         (sol_LWFBrook90[13 + 3,:], sol_LWFBrook90[24 + 3,:]),
# #         (sol_LWFBrook90[13 + 4,:], sol_LWFBrook90[24 + 4,:])
# #     ];
# #     xlabel = "δ18O [‰]", ylabel = "δ2H [‰]",
# #     label = ["GWAT" "INTS" "INTR" "SNOW"])
# # # plot!(ylims = (-100,-60), xlims = (-15, -8))
# # # plot!(sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), sol_LWFBrook90.prob.p[2][16].(sol_LWFBrook90.t), label = "PREC")

# # # # Plot 2
# # # # http://docs.juliaplots.org/latest/generated/gr/#gr-ref43
# # x = LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date)
# # y = cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK)
# # n = sol_LWFBrook90.prob.p[1][1].NLAYER
# # y_centers = [ 0; cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK)[1:(n-1)] ] +
# #     sol_LWFBrook90.prob.p[1][1].p_THICK / 2

# # z = sol_LWFBrook90[7 .+ (0:sol_LWFBrook90.prob.p[1][1].NLAYER-1),
# #                     1,
# #                     :]./sol_LWFBrook90.prob.p[1][1].p_THICK;
# # z2 = sol_LWFBrook90[idx_u_vector_isotopes_d18O,1,:];
# # z3 = sol_LWFBrook90[idx_u_vector_isotopes_d2H,1,:];

# # heatmap(x, y_centers, z, yflip = true,
# #         xlabel = "Date",
# #         ylabel = "Depth [mm]",
# #         colorbar_title = "θ [-]")
# # hline!([0; cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK)], yflip = true, xticks = false,
# #     color = :black, linestyle = :dot
# #     #title = "N_layer = "*string(sol_LWFBrook90.prob.[1][1].NLAYER)
# #     )
# # hline!(y_centers, yflip = true, xticks = false,
# #     color = :blue, linestyle = :dot)
# # heatmap(x, y, z2, yflip = true,
# #         xlabel = "Date",
# #         ylabel = "Depth [mm]",
# #         colorbar_title = "δ18O [‰]")
# # heatmap(x, y, z3, yflip = true,
# #         xlabel = "Date",
# #         ylabel = "Depth [mm]",
# #         colorbar_title = "δ2H [‰]")

# # TODO: edges of cells in heatmap are not entirely correct. Find a way to override heatmap()
# #       where we provide cell edges (n+1) instead of cell centers (n)
# # TODO: e.g. plots_heatmap_edges: @recipe function f(::Type{Val{:plots_heatmap_edges}}, xe, ye, z)
# # TODO: e.g. plots_heatmap_edges:     m, n = size(z.surf)
# # TODO: e.g. plots_heatmap_edges:     x_pts, y_pts = fill(NaN, 6 * m * n), fill(NaN, 6 * m * n)
# # TODO: e.g. plots_heatmap_edges:     fz = zeros(m * n)
# # TODO: e.g. plots_heatmap_edges:     for i in 1:m # y
# # TODO: e.g. plots_heatmap_edges:         for j in 1:n # x
# # TODO: e.g. plots_heatmap_edges:             k = (j - 1) * m + i
# # TODO: e.g. plots_heatmap_edges:             inds = (6 * (k - 1) + 1):(6 * k - 1)
# # TODO: e.g. plots_heatmap_edges:             x_pts[inds] .= [xe[j], xe[j + 1], xe[j + 1], xe[j], xe[j]]
# # TODO: e.g. plots_heatmap_edges:             y_pts[inds] .= [ye[i], ye[i], ye[i + 1], ye[i + 1], ye[i]]
# # TODO: e.g. plots_heatmap_edges:             fz[k] = z.surf[i, j]
# # TODO: e.g. plots_heatmap_edges:         end
# # TODO: e.g. plots_heatmap_edges:     end
# # TODO: e.g. plots_heatmap_edges:     ensure_gradient!(plotattributes, :fillcolor, :fillalpha)
# # TODO: e.g. plots_heatmap_edges:     fill_z := fz
# # TODO: e.g. plots_heatmap_edges:     line_z := fz
# # TODO: e.g. plots_heatmap_edges:     x := x_pts
# # TODO: e.g. plots_heatmap_edges:     y := y_pts
# # TODO: e.g. plots_heatmap_edges:     z := nothing
# # TODO: e.g. plots_heatmap_edges:     seriestype := :shape
# # TODO: e.g. plots_heatmap_edges:     label := ""
# # TODO: e.g. plots_heatmap_edges:     widen --> false
# # TODO: e.g. plots_heatmap_edges:     ()
# # TODO: e.g. plots_heatmap_edges: end
# # TODO: e.g. plots_heatmap_edges: @deps plots_heatmap_edges shape
# # TODO: e.g. plots_heatmap_edges: @shorthands plots_heatmap_edges
# # TODO: e.g. plots_heatmap_edges:
# # TODO: e.g. plots_heatmap_edges: Plots.heatmap(x[1:100], y_centers, z[:,1:100])
# # TODO: e.g. plots_heatmap_edges: Plots.heatmap(x[1:100], y_centers, z[:,1:100])
# # TODO: e.g. plots_heatmap_edges: plot(t = :heatmap, x[1:50], y_centers, z[:,1:50]) # works
# # TODO: e.g. plots_heatmap_edges: plot(t = :plots_heatmap, x[1:50], y_centers, z[:,1:50]) # doesn't work
# # TODO: e.g. plots_heatmap_edges: plot(t = :plots_heatmap_edges, x[1:50], y_centers, z[:,1:50]) # doesn't work either
# # # ###################

# # plot(x,
# #      z[LWFBrook90.find_indices(depth_to_read_out_mm, sol_LWFBrook90), :]',
# #      labels = string.(depth_to_read_out_mm) .* "mm",
# #      xlabel = "Date",
# #      ylabel = "θ [-]",
# #      legend = :bottomright)
# # savefig(input_prefix*"_θ-feinUndSteinErde_depths_NLAYER"*string(sol_LWFBrook90.prob.p[1][1].NLAYER)*".png")





# ###################################
# ## More complex plots (comparing with LWFBrook90R results if available):
# using Statistics
# using Plots, Dates, StatsPlots, DataFramesMeta
# using Gadfly, Cairo #,Fontconfig
# using CSV

using CSV: read, File
using Dates: DateTime, Millisecond, Second, Day, Month, month, value, dayofyear, format
using Plots
using DataFrames: DataFrame, rename, sort!# ,select
using DataFramesMeta#: @linq, transform, DataFramesMeta
# input_meteoveg = @linq DataFrame(File(path_meteoveg;
#     skipto=3, delim=',', ignorerepeated=false,
#     # Be strict about loading NA's -> error if NA present
#     types=parsing_types, missingstring = nothing, strict=true))  |>
#     transform(:dates = DateTime.(:dates))

if (input_prefix == "DAV2010-2021" || input_prefix == "DAV_LW1_def");
    # File("../../../LWF-Brook90.jl-calibration/Isotope-Data/2-combined-deposition/outputs/DAV2010-2021_meteoiso.csv")
    # File("../../../LWF-Brook90.jl-calibration/Meteo-Data/DAV2010-2021/DAV2010-2021_meteoiso.csv")
    dat_δ_soilSol   = DataFrame(File("../../../LWF-Brook90.jl-calibration/Isotope-Data/3-LWF-soil-solution/outputs/DAV-2022-07_soilsolutioniso.csv";header=1,skipto=3))
    dat_δ_soilbulk  = DataFrame(File("../../../LWF-Brook90.jl-calibration/Isotope-Data/4-LWF-xylem-bulksoil/outputs/DAV-2022-07_bulksoiliso.csv";header=1,skipto=3))
    dat_δ_xylem     = DataFrame(File("../../../LWF-Brook90.jl-calibration/Isotope-Data/4-LWF-xylem-bulksoil/outputs/DAV-2022-07_xylemiso.csv";header=1,skipto=3))
    dat_ψ_tensiometer= DataFrame(File("../../../LWF-Brook90.jl-calibration/Soil-Data/1-LWF-tensiometer/DAV-2022-07_tensiometer.csv";header=1,skipto=3))
    dat_θ_EC5Waldner = DataFrame(File("../../../LWF-Brook90.jl-calibration/Soil-Data/2-LWF-EC5-Waldner/outputs/DAV-2022-07_ec5-Waldner.csv";header=1,skipto=3))
    dat_ψ_tensiomark = DataFrame(File("../../../LWF-Brook90.jl-calibration/Soil-Data/3-LWF-Tensiomark-EC5-Meusburger/outputs/DAV-2022-07_tensiomark-Meusburger.csv";header=1,skipto=3))
    dat_ψ_tensiomarkAVG= DataFrame(File("../../../LWF-Brook90.jl-calibration/Soil-Data/3-LWF-Tensiomark-EC5-Meusburger/outputs/DAV-2022-07_tensiomarkAVG-Meusburger.csv";header=1,skipto=3))

    # # a) Using Wide data sets:
    # df = unstack(dat_ψ_tensiomarkAVG, :depth_cm, :psi_kPa)
    # plot(df.dates,
    #      Matrix(df[:,Not(:dates)]),
    #      labels = permutedims(names(df[:,Not(:dates)])).*"cm",
    #      legend = :bottomright, marker = :dot)

    # df = unstack(dat_ψ_tensiomark, :depth_cm, :psi_kPa; allowduplicates=true)
    # plot(df.dates,
    #      Matrix(df[:,Not(:dates)]),
    #      labels = permutedims(names(df[:,Not(:dates)])).*"cm",
    #      legend = :bottomright, marker = :dot)

    # # b) Keeping Long data sets:
    # scatter(dat_ψ_tensiomark.dates, dat_ψ_tensiomark.psi_kPa, group = dat_ψ_tensiomark.depth_cm)
    # scatter(dat_ψ_tensiomarkAVG.dates, dat_ψ_tensiomarkAVG.psi_kPa, group = dat_ψ_tensiomarkAVG.depth_cm)
    # scatter(dat_θ_EC5Waldner.dates, dat_θ_EC5Waldner.theta_m3m3, group = dat_θ_EC5Waldner.depth_cm)
    # scatter(dat_ψ_tensiometer.dates, dat_ψ_tensiometer.psim_median_kPa, group = dat_ψ_tensiometer.depth_cm)

    @chain dat_ψ_tensiomark    scatter(_.dates, _.psi_kPa,         group = _.depth_cm)
    @chain dat_ψ_tensiomarkAVG scatter(_.dates, _.psi_kPa,         group = _.depth_cm)
    @chain dat_θ_EC5Waldner    scatter(_.dates, _.theta_m3m3,      group = _.depth_cm)
    @chain dat_ψ_tensiometer   scatter(_.dates, _.psim_median_kPa, group = _.depth_cm)
    @chain dat_δ_soilSol       scatter(_.dates, _.delta18O_permil, group = _.depth_cm)
    @chain dat_δ_soilbulk      scatter(_.dates, _.delta18O_permil, group = _.depth_cm)
    @chain dat_δ_xylem         scatter(_.dates, _.delta18O_permil, group = _.treeID)

    depth_to_read_out_mm = [150 500 800 1500]

    @chain dat_ψ_tensiometer[1:100,:]   scatter(convert.(DateTime, _.dates), _.psim_median_kPa, group = _.depth_cm .* 10)
    pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
        # -LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90) .+ 1, yaxis = :log, yflip = true,
        LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90),  ylabel = "ψ\n[kPa]",
        # log10.(-10 .* LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90)),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
        labels = string.(depth_to_read_out_mm) .* "mm",
        xlabel = "Date",
        legend = :outerright)

    @chain dat_θ_EC5Waldner[1:900,:]    scatter(convert.(DateTime, _.dates), _.theta_m3m3,      group = _.depth_cm)
    pl_θ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
        LWFBrook90.get_θ(depth_to_read_out_mm, sol_LWFBrook90),
        labels = string.(depth_to_read_out_mm) .* "mm",
        xlabel = "Date",
        ylabel = "θ\n[-]",
        legend = :outerright)
    xlims!(extrema(convert.(DateTime, dat_θ_EC5Waldner[1:900,:].dates)))

    @chain dat_ψ_tensiometer[1:100,:]   scatter!(_.dates, _.psim_median_kPa, group = _.depth_cm)
    @chain dat_ψ_tensiomarkAVG scatter!(
        _.dates,
        log10.(-10 .* _.psi_kPa),
        group = _.depth_cm)


    # Averages with nominal depth (or rounded depth actually)
    depth_to_read_out_mm = sort(unique(dat_ψ_tensiomarkAVG.depth_nominal_cm)).*10
    col_palette = palette(:inferno, length(depth_to_read_out_mm)+1; rev = true)[Not(1)]
    @chain dat_ψ_tensiomarkAVG scatter(
            convert.(DateTime, _.dates),
            _.psi_kPa,
            group = _.depth_nominal_cm,
            # color_palette = [:blue, :red, :green, :yellow]
            # color_palette = palette(:inferno, 4; rev = true),
            # color_palette = palette(:inferno, 4+2; rev = true)[Not(1,end)],
            # color_palette = palette(:inferno, 4+1; rev = true)[Not(1)],
            color_palette = col_palette,
            marker = (0.3, :o), markerstrokewidth = 0)
    pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
        # -LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90) .+ 1, yaxis = :log, yflip = true,
        LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90),  ylabel = "ψ\n[kPa]",
        # log10.(-10 .* LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90)),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
        labels = reshape(string.(depth_to_read_out_mm) .* "mm",1,:),
        xlabel = "Date", legend = :outerright, color_palette = col_palette,
        linewidth = 3)
    xlims!(extrema(convert.(DateTime, dat_ψ_tensiomarkAVG.dates)));
    ylims!(-15,0);
    plot!(size=(800,500), dpi = 200, leftmargin = 5mm)
    savefig(fname*"_vs-ψ-average-data.png")

    @chain dat_ψ_tensiomarkAVG scatter(
            convert.(DateTime, _.dates),
            log10.(-10 .* _.psi_kPa), yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", # _.psi_kPa,
            group = _.depth_nominal_cm,
            color_palette = col_palette,
            marker = (0.3, :o), markerstrokewidth = 0)
    pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
        log10.(-10 .* LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90)),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
        labels = reshape(string.(depth_to_read_out_mm) .* "mm",1,:),
        xlabel = "Date", legend = :outerright, color_palette = col_palette,
        linewidth = 3)
    xlims!(extrema(convert.(DateTime, dat_ψ_tensiomarkAVG.dates)))
    ylims!(0,3) #ylims!(log10(-10 * -0), log10(-10 * -15)); # ylims!(-15,0);
    plot!(size=(800,500), dpi = 200, leftmargin = 5mm)
    savefig(fname*"_vs-ψ-average-data-log.png")

    # Each sensore with actual depth
    depth_to_read_out_mm = sort(unique(dat_ψ_tensiomark.depth_cm)).*10
    col_palette = palette(:inferno, length(depth_to_read_out_mm)+1; rev = true)[Not(1)]
    @chain dat_ψ_tensiomark scatter(
            convert.(DateTime, _.dates),
            _.psi_kPa,
            group = _.depth_cm,
            # color_palette = [:blue, :red, :green, :yellow]
            # color_palette = palette(:inferno, 4; rev = true),
            # color_palette = palette(:inferno, 4+2; rev = true)[Not(1,end)],
            # color_palette = palette(:inferno, 4+1; rev = true)[Not(1)],
            color_palette = col_palette,
            marker = (0.3, :o), markerstrokewidth = 0)

    pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
        # -LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90) .+ 1, yaxis = :log, yflip = true,
        LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90),  ylabel = "ψ\n[kPa]",
        # log10.(-10 .* LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90)),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
        labels = reshape(string.(depth_to_read_out_mm) .* "mm",1,:),
        xlabel = "Date", legend = :outerright, color_palette = col_palette,
        linewidth = 3)
    xlims!(extrema(convert.(DateTime, dat_ψ_tensiomarkAVG.dates)));
    ylims!(-15,0);
    plot!(size=(800,500), dpi = 200, leftmargin = 5mm)
    savefig(fname*"_vs-ψ-data.png")



    @chain dat_ψ_tensiomark scatter(
            convert.(DateTime, _.dates),
            log10.(-10 .* _.psi_kPa), yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", # _.psi_kPa,
            group = _.depth_cm,
            color_palette = col_palette,
            marker = (0.3, :o), markerstrokewidth = 0)
    pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
        log10.(-10 .* LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90)),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
        labels = reshape(string.(depth_to_read_out_mm) .* "mm",1,:),
        xlabel = "Date", legend = :outerright, color_palette = col_palette,
        linewidth = 3)
    xlims!(extrema(convert.(DateTime, dat_ψ_tensiomarkAVG.dates)));
    ylims!(0,3) #ylims!(log10(-10 * -0), log10(-10 * -15)); # ylims!(-15,0);

    plot!(size=(800,500), dpi = 200, leftmargin = 5mm)
    savefig(fname*"_vs-ψ-data_log.png")




    ### Tensiometer:

    depth_to_read_out_mm = sort(unique(dat_ψ_tensiometer.depth_cm)).*10
    col_palette = palette(:inferno, length(depth_to_read_out_mm)+1; rev = true)[Not(1)]

    @chain dat_ψ_tensiometer   scatter(
            convert.(DateTime, _.dates),
            # log10.(-10 .* _.psim_median_kPa), yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", # _.psi_kPa,
            _.psim_median_kPa,
            group = _.depth_cm,
            color_palette = col_palette,
            marker = (0.3, :o), markerstrokewidth = 0)
    pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
        # log10.(-10 .* LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90)),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
        LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90),  ylabel = "ψ [kPa]",
        labels = reshape(string.(depth_to_read_out_mm) .* "mm",1,:),
        xlabel = "Date", legend = :outerright, color_palette = col_palette,
        linewidth = 3)
    ylims!(-15,0)

    plot!(size=(800,500), dpi = 200, leftmargin = 5mm)
    savefig(fname*"_vs-ψ-tensiometer.png")




    ### Soilsolution isotopes:
    depth_to_read_out_mm = sort(unique(dat_δ_soilSol.depth_cm)).*10
    col_palette = palette(:inferno, length(depth_to_read_out_mm)+1; rev = true)[Not(1)]

    @chain dat_δ_soilSol   scatter(
            convert.(DateTime, _.dates),
            _.delta18O_permil,
            group = _.depth_cm,
            color_palette = col_palette,
            marker = (0.7, :o), markerstrokewidth = 0)
    pl_δsoil = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
        LWFBrook90.get_δsoil(depth_to_read_out_mm  .+ 0.01, sol_LWFBrook90).d18O,
        ylabel = "δ18O [permil]",
        labels = reshape(string.(depth_to_read_out_mm) .* "mm",1,:),
        xlabel = "Date", legend = :outerright, color_palette = col_palette,
        linewidth = 3, alpha = 0.5)
    # add precipitation
    pl_δsoil = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
                    sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t),
                    color = "grey", alpha = 0.5, labels = "Precipiation")

    # ylims!(-15,0)
    xlims!(extrema(convert.(DateTime, dat_δ_soilSol.dates)))
    plot!(size=(800,500), dpi = 200, leftmargin = 5mm)
    savefig(fname*"_vs-δ-soilsolution.png")


    # ### Bulksoil isotopes:
    # # TODO
    # dat_for_heatmap = @chain begin
    #     dat_δ_soilbulk
    #     @subset (:treeID .== "240")
    #     @select($(Not(:delta2H_permil)))
    #     @orderby :depth_cm
    #     unstack(:depth_cm, :delta18O_permil)
    #     @orderby :dates
    # end


    # Plots.heatmap(
    #     dat_for_heatmap.dates, # 1:5,
    #     parse.(Int, names(dat_for_heatmap)[Not(1,2)]),
    #     permutedims(Matrix(dat_for_heatmap[:,Not(1,2)])),
    #     yflip = true
    #     )

    ### Xylem isotopes:
    @chain dat_δ_xylem   scatter(
            convert.(DateTime, _.dates),
            _.delta18O_permil,
            group = _.treeID,
            marker = (0.7, :o), markerstrokewidth = 0)
    pl_δsoil = plot!(
        LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
        permutedims([LWFBrook90.get_δ(sol_LWFBrook90).RWU.d18O;
                    LWFBrook90.get_δ(sol_LWFBrook90).XYL.d18O]),
        ylabel = "δ18O [permil]",
        labels = ["RWU" "Xylem"],
        xlabel = "Date", legend = :outerright,
        linewidth = 3, alpha = 0.5)
    # add precipitation
    pl_δsoil = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
                    sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t),
                    color = "grey", alpha = 0.5, labels = "Precipiation")

    # ylims!(-15,0)
    xlims!(extrema(convert.(DateTime, dat_δ_xylem.dates)) .+ (-Day(100), Day(50)))
    plot!(size=(800,500), dpi = 200, leftmargin = 5mm)
    savefig(fname*"_vs-δ-xylem.png")
end

