using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5, init#, step!
# example = LWFBrook90.run_example()

# Read in input data
# input_prefix = "isoBEA2016-reset-FALSE"
# input_path = "examples/isoBEA2016-reset-FALSE-input/"
# input_prefix = "isoBEA2010-18-reset-FALSE";
# input_path = "examples/isoBEA2010-18-reset-FALSE-input/";
input_prefix = "isoBEAdense2010-18-reset-FALSE";
input_path = "examples/isoBEAdense2010-18-reset-FALSE-input/";
# input_prefix = "BEA2010-2021";
# input_path = "../../../LWF-Brook90.jl-calibration/Meteo-Data/BEA2010-2021/";

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
Δz_m = [fill(0.02, 10); fill(0.025, 10); fill(0.03, 10); fill(0.035, 10); 0.1]; # grid spacing (heterogenous), meter (N=41)
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
    input_param[1,"NOOUTF"]
    input_soil_horizons.gravel_volFrac .= 0.00

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

    # Define ODE:
    ode_LWFBrook90, unstable_check_function = define_LWFB90_ODE(u0, tspan, p);
    ####################

    ####################
    ## Solve ODE:
    # @run sol_LWFBrook90 = solve(ode_LWFBrook90, Tsit5();
    #     unstable_check = (dt,u,p,t) -> false,#any(isnan,u),
    #     progress = true,
    #     saveat = tspan[1]:tspan[2], dt=1e-6, adaptive = true); # dt is initial dt, but adaptive
    ####################

    ####################
    ## Benchmarking
    # @time sol_LWFBrook90 = solve(ode_LWFBrook90, Tsit5(); progress = true,
    #     unstable_check = unstable_check_function, # = (dt,u,p,t) -> false, #any(isnan,u),
    #     saveat = tspan[1]:tspan[2], dt=1e-6, dtmax=1e-3, adaptive = true);
                    # integrator =  init(ode_LWFBrook90, Tsit5(); progress = true,
                    #             unstable_check = unstable_check_function) # = (dt,u,p,t) -> false, #any(isnan,u))

                    # du = copy(integrator.u)
                    # du .= 0
                    # LWFBrook90.f_LWFBrook90R(du, integrator.u, p, 0);
                    # du
                    # step!(integrator)
    @time sol_LWFBrook90 = solve(ode_LWFBrook90, Tsit5(); progress = true,
        unstable_check = unstable_check_function, # = (dt,u,p,t) -> false, #any(isnan,u),
        saveat = tspan[1]:tspan[2],
        # When adding transport of isotopes adaptive time stepping has difficulties reducing dt below dtmin.
        # We solve it either by using a constant time step
        # or by setting force_dtmin to true, which has the drawback that tolerances are ignored.
        # # Variant1: less sophisticated but robust:
        # adaptive = false, dt=1e-3 #dt=5/60/24 # fixed 5 minutes time steps
        # Variant2: more sophisticated but giving errors:
        adaptive = true, dtmin = 1e-4, dt=1e-4,  # if adaptive dt is just the starting value of the time steps
        force_dtmin = true,# without this callbacks generate an abort "dt <= dtmin",
                           # e.g. see: https://github.com/SciML/DifferentialEquations.jl/issues/648
        maxiters = (tspan[2]-tspan[1])/1e-4 # when using force_dtmin also use maxiters
                        # TODO(bernhard): regarding maxiters it seems to be the callback for δ of SNOW, INTR, INTS that causes the dtmin to be so small to reach large maxiters
        );
                using Plots, Measures
                t_ref = sol_LWFBrook90.prob.p[2][17]
                x = RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, t_ref);
                y = cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK) - sol_LWFBrook90.prob.p[1][1].p_THICK/2
                optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
                y_soil_ticks = optim_ticks(0., round(maximum(cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
                y_ticks= y_soil_ticks
                y_labels=round.(y_soil_ticks; digits=0)
                (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
                                LWFBrook90.get_auxiliary_variables(sol_LWFBrook90)
                # heatmap(#x[140:150], y, u_aux_θ'[:,140:150], yflip = true,
                #         x[142:145], y, u_aux_θ'[:,142:145], yflip = true,
                #         xlabel = "Date",
                #         ylabel = "Depth [mm]",
                #         colorbar_title = "θ [-]")
                t_to_plot = 1:length(x)
                plot(
                    heatmap(x[t_to_plot], y, u_aux_θ'[:, t_to_plot], yflip = true,
                            xlabel = "Date",
                            ylabel = "Depth [mm]",
                            colorbar_title = "θ [-]"),
                    heatmap(x[t_to_plot], y, u_aux_PSIM'[:, t_to_plot], yflip = true,
                            xlabel = "Date",
                            ylabel = "Depth [mm]",
                            colorbar_title = "ψ_m [kPa]"),
                    layout = (2,1))
                # du = similar(sol_LWFBrook90.u[142])
                # LWFBrook90.f_LWFBrook90R(du, sol_LWFBrook90.u[142],p,142)
                # integrator = init(ode_LWFBrook90, Tsit5(), unstable_check = unstable_check_function,
                #                 saveat = 0:145, adaptive = true, dtmin = 1e-5, dt=1e-4)
                # step!(integrator)
                # function plot_θ!(pl, integrator)
                #     p_soil = integrator.p[1][1]
                #     idx_u_vector_amounts = integrator.p[1][4].row_idx_SWATI
                #     t_ref = integrator.p[2][17]
                #     # u_SWATI = [solution.u[i][idx_u_vector_amounts] for i = 1:size(integrator.u,2)]
                #     u_SWATI = [integrator.u[idx_u_vector_amounts]]
                #     (_, _, _, u_aux_θ, _) = LWFBrook90.KPT.derive_auxiliary_SOILVAR(u_SWATI[1], p_soil)
                #     heatmap!(pl, RelativeDaysFloat2DateTime.(integrator.t .+ [0; 1], t_ref),
                #             cumsum(integrator.p[1][1].p_THICK) - integrator.p[1][1].p_THICK/2,
                #             [u_aux_θ u_aux_θ], yflip = true)
                #     return pl
                # end
                # for i in 1:140 step!(integrator) end
                # pl = plot()
                # plot_θ!(pl, integrator)
                # begin
                #     step!(integrator)
                #     plot_θ!(pl, integrator)
                # end
                # du

                # # investigate what happens on 2010-05-22 in BEA (artifact in θ simulation)?
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
    # heatmap(x, y, u_aux_PSIM', yflip = true,
    #         xlabel = "Date",
    #         ylabel = "Depth [mm]",
    #         colorbar_title = "ψ [kPa]")


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

    # @time sol_LWFBrook90 = solve(ode_LWFBrook90, progress = true, Euler(); # Note: Euler sometimes hangs
    #     saveat = tspan[1]:tspan[2], dt=1e-1, adaptive = false);
    # using BenchmarkTools # for benchmarking
    # sol_LWFBrook90 = @btime solve(ode_LWFBrook90; dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity
    # sol_LWFBrook90 = @btime solve(ode_LWFBrook90; saveat = tspan[1]:tspan[2], dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity
    ####################

    ####################
    ## Plotting
    # 0) defined in module_ISO.jl
    # using Plots, Measures
    # include("module_ISOplots.jl")
    # pl_final_δ18O, pl_final_δ2H =
    # plot_LWFBrook90_isotopes(sol_LWFBrook90; clims_d18O = (-16, -6), clims_d2H  = (-125, -40));
    # plot(pl_final_δ18O, pl_final_δ2H,
    #     layout = (2,1), size=(1000,1400), leftmargin = 15mm);
    # savefig(input_prefix*".png")
    mkpath("out")

    using DataFrames
    using CSV: File
    using Plots, Measures
    fname = joinpath(
        "out",
        input_prefix*"_NLAYER+" * string(sol_LWFBrook90.prob.p[1][1].NLAYER)*
        "-git+"*chomp(read(`git rev-parse --short HEAD`, String))*
        ifelse(length(read(`git status --porcelain`, String))==0, "+gitclean","+gitdirty")*
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

                # 3b) Heatmap of RWU
                NLAYER = sol_LWFBrook90.prob.p[1][1].NLAYER
                rows_RWU      = sol_LWFBrook90[sol_LWFBrook90.prob.p[1][4].row_idx_RWU, 1, :]
                # NOTE: sol_LWFBrook90[:, :, 1] == u0# is true
                # # y_extended = [-500; -350; -300; -250; -200; -150; -100;   -50;         y;          (maximum(y) .+ 50 .+ [50; 100; 150; 250])]
                # # y_labels   = ["INTS"; ""; "INTR"; ""; "SNOW"; ""; round.(y); "";             "GWAT"]
                # # y_soil_ticks = optimize_ticks(extrema(y)...; k_min = 4)[1]
                # # y_soil_ticks = optimize_ticks(0., round(maximum(cumsum(sol.prob.p[1][1].p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
                # y_soil_ticks = tick_function(0., round(maximum(cumsum(sol.prob.p[1][1].p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
                # y_ticks    = [-500;       -300;       -200;       -100;     y_soil_ticks;          (maximum(y) .+ 50 .+ [    100;    250])]
                # y_labels   = ["PREC";   "INTS";     "INTR";     "SNOW";     round.(y_soil_ticks; digits=0); "GWAT"  ; "RWU"]
                # z2_extended = [row_PREC_d18O; row_NaN; row_INTS_d18O; row_NaN; row_INTR_d18O; row_NaN; row_SNOW_d18O; row_NaN; rows_SWAT_d18O; row_NaN; row_GWAT_d18O; row_NaN; row_RWU_d18O]
                # z3_extended = [row_PREC_d2H;  row_NaN; row_INTS_d2H;  row_NaN; row_INTR_d2H;  row_NaN; row_SNOW_d2H;  row_NaN; rows_SWAT_d2H;  row_NaN; row_GWAT_d2H;  row_NaN; row_RWU_d2H]
                (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
                    LWFBrook90.get_auxiliary_variables(sol_LWFBrook90)
                # z = sol_LWFBrook90[7 .+ (0:sol_LWFBrook90.prob.p[1][1].NLAYER-1), 1, :]./sol_LWFBrook90.prob.p[1][1].p_THICK

                t_ref = sol_LWFBrook90.prob.p[2][17]
                x = RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, t_ref);
                y = cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK) - sol_LWFBrook90.prob.p[1][1].p_THICK/2
                y_soil_ticks = optim_ticks(0., round(maximum(cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
                y_ticks= y_soil_ticks
                y_labels=round.(y_soil_ticks; digits=0)
                heatmap(x, y, u_aux_θ', yflip = true,
                        xlabel = "Date",
                        ylabel = "Depth [mm]",
                        colorbar_title = "θ [-]")
                heatmap(x, y, u_aux_PSIM', yflip = true,
                        xlabel = "Date",
                        ylabel = "Depth [mm]",
                        colorbar_title = "ψ [kPa]")
                z=rows_RWU
                pl_RWU = heatmap(x,y,z; yticks = (y_ticks, y_labels),
                        yflip = true,
                        c = :diverging_bwr_20_95_c54_n256, clim = maximum(abs.(z)) .* (-1, 1),
                        ylabel = "Depth [mm]",
                        colorbar_title = "RWU [mm/day]",
                        size=(1400,800), dpi = 300, leftmargin = 15mm);
                # savefig(pl_RWU, fname*"_RWU.png")
                plot_avgRWU = plot(sum(z,dims=2), y,
                    yflip = true, ylabel = "Depth [mm]", xlabel = "Total RWU over\nsimulation period [mm]")
                # Plot root distribution
                p_RELDEN = sol_LWFBrook90.prob.p[2].p_RELDEN
                root_data_start_end = [p_RELDEN.(t, 1:NLAYER) for t in extrema(sol_LWFBrook90.prob.tspan)]
                pl_roots = plot(root_data_start_end, y,
                    labels = ("t=" .* string.(permutedims(collect(extrema(sol_LWFBrook90.prob.tspan))))),  linestyle = [:solid :dash],
                    yflip = true, ylabel = "Depth [mm]", xlabel = "Relative root\ndensity [-]");
                pl_RWU2 = plot(
                    plot(pl_roots;                 yticks = (y_ticks, y_labels)),
                    plot(plot_avgRWU, ylabel = ""; yticks = []),
                    plot(pl_RWU;      ylabel = "", yticks = [], title = "RWU [mm/day]"),
                    layout=grid(1,3, widths = [0.1,0.1,0.8]), link = :y,
                    size=(1400,800), dpi = 300, leftmargin = 15mm, bottommargin = 10mm)
                savefig(pl_RWU2, fname*"_RWU2.png")


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


    if true
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
            ylabel = "θ\n[-]",
            legend = :outerright);
        pl_ψ = plot(LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date),
            # -LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90) .+ 1, yaxis = :log, yflip = true,
            LWFBrook90.get_ψ(depth_to_read_out_mm, sol_LWFBrook90),
            labels = string.(depth_to_read_out_mm) .* "mm",
            xlabel = "Date",
            ylabel = "ψ\n[kPa]",
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
