using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5, init#, step!
# example = LWFBrook90.run_example()
# using Plots
using Dates
using DataFrames
using Measures
using CSV: read, File
using Dates: DateTime, Millisecond, Second, Day, Month, month, value, dayofyear, format, today
using DataFrames: DataFrame, rename, sort!# ,select
using DataFramesMeta
using Plots
        # # TODO: make a test out of this:
        # parametrizedSPAC = loadSPAC("../../../LWF-Brook90.jl-calibration/Input-Meteo-Data/DAV_LW1_def/", "DAV_LW1_def"; simulate_isotopes = true);
        # Δz = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 5)]; # grid spacing (heterogenous), meter (N=20)
        # soil_output_depths_m = zeros(Float64, 0)
        # parametrizedSPAC.tspan
        # DS = LWFBrook90.setup(parametrizedSPAC::SPAC; Δz = Δz,  tspan = (0,10));
        # # Test the f and cb() functions (LWFBrook90.f_LWFBrook90!, )
        # t0 = 0.0
        # u0 = DS.ODEProblem.u0
        # p = DS.ODEProblem.p
        # du = copy(u0)
        # @btime states = NamedTuple{p[1][4].u0_field_names}(u0.x);

        # LWFBrook90.f_LWFBrook90!(du, u0, p, 1.0)
        # # integrator = init(DS.ODEProblem, Tsit5();)
        # # LWFBrook90.LWFBrook90R_updateAmounts_INTS_INTR_SNOW_CC_SNOWLQ!(integrator)
        # # LWFBrook90.LWFBrook90R_updateIsotopes_INTS_INTR_SNOW!(integrator)
        # # LWFBrook90.LWFBrook90R_updateIsotopes_GWAT_SWAT!(u0, t0, integrator)
        # # LWFBrook90.LWFBrook90R_updateIsotopes_GWAT_SWAT_AdvecDiff!(u0, t0, integrator)
        # # LWFBrook90.LWFBrook90R_check_balance_errors!(integrator)
        # # # @run simulate!(DS) # with ArrayParition this simulating 10 days takes about 0.373707 seconds (1.45 M allocations: 105.396 MiB, 14.36% gc time)
        # DS.ODEProblem.tspan
        # simulate!(DS); # with ArrayParition this simulating 10 days takes about 0.373707 seconds (1.45 M allocations: 105.396 MiB, 14.36% gc time) Time steps for solving: 580 (580 accepted out of 612 total)
        # # simulate!(DS); # with ArrayParition this simulating 10 days takes about 0.373707 seconds (1.45 M allocations: 105.396 MiB, 14.36% gc time) Time steps for solving: 580 (580 accepted out of 612 total)
        # DS = LWFBrook90.setup(parametrizedSPAC::SPAC; Δz = Δz);
        # DS.ODEProblem.tspan
        # # simulate!(DS); # with ArrayParition this simulating all the days takes about 167.115807 seconds (762.30 M allocations: 53.365 GiB, 6.18% gc time, 12.43% compilation time) Time steps for solving: 295894 (295894 accepted out of 314371 total)

function run_main_with_isotopes(;input_prefix, input_path)
    ####################
    # Define simulation model by reading in system definition and input data
    simulate_isotopes = true
    model = loadSPAC(input_path, input_prefix; simulate_isotopes = simulate_isotopes);
    model = loadSPAC(input_path, input_prefix; simulate_isotopes = simulate_isotopes,
                    IC_soil = (PSIM_init_kPa = -6.3, delta18O_init_permil = -13.0, delta2H_init_permil = -95.0),
                    root_distribution = (beta = 0.90, z_rootMax_m = nothing));
    model2 = loadSPAC(input_path, input_prefix;
                simulate_isotopes = false,
                IC_soil = (PSIM_init_kPa = -6.3, delta18O_init_permil = -13.0, delta2H_init_permil = -95.0),
                root_distribution = (beta = 0.97, z_rootMax_m = nothing),
                canopy_evolution = (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel    = 100,
                                    LAI_rel = (DOY_Bstart = 120,
                                           Bduration  = 21,
                                           DOY_Cstart = 270,
                                           Cduration  = 60,
                                           LAI_perc_BtoC = 100,
                                           LAI_perc_CtoB = 20)));
    simulation  = setup(model; soil_output_depths_m = [-0.03, -0.11]);

    # Prepare simulation by discretizing spatial domain
    Δz_m = nothing
    if (input_prefix == "isoBEAdense2010-18-reset-FALSE" || contains(input_prefix, r"infiltrationSaturation"))
        Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 5); 0.1]; # grid spacing (heterogenous), meter (N=21)
        Δz_m = [fill(0.02, 10); fill(0.025, 10); fill(0.03, 10); fill(0.035, 10); 0.1]; # grid spacing (heterogenous), meter (N=21)
    elseif (input_prefix == "DAV_LW1_def")
        # Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 5)]; # grid spacing (heterogenous), meter (N=20)
        Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 11)]; # grid spacing (heterogenous), meter (N=20)
    else #"SCH-2023-03"
        Δz_m = [fill(0.05, 42);]; # grid spacing (heterogenous), meter (N=20)
    end
    model3 = loadSPAC(input_path, input_prefix;
        Δz_thickness_m = Δz_m,
        root_distribution = (beta = 0.98, ),
        IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -9.0, delta2H_init_permil = -11.0))
    simulation3 = setup(model3, soil_output_depths_m = [-0.03, -0.11]);
    model = loadSPAC(input_path, input_prefix; simulate_isotopes = simulate_isotopes,
                    Δz_thickness_m = Δz_m,
                    IC_soil = (PSIM_init_kPa = -6.3, delta18O_init_permil = -13.0, delta2H_init_permil = -95.0),
                    root_distribution = (beta = 0.90, z_rootMax_m = nothing));
    simulation  = setup(model; soil_output_depths_m = [-0.03, -0.11], requested_tspan = (0,300));
    ####################

    ####################
    # using Plots; display(heatmap(simulation.ODEProblem.p.p_fT_RELDEN', ylabel = "SOIL LAYER", xlabel = "Time (days)", yflip=true, colorbar_title = "Root density"))
        # # TODO
        # remake = function(parameterisedSPAC::SPAC; Δz_thickness_m::AbstractVector)
        #     # Redefine resolution, without providing parametrized functions for IC and root density distribution
        #     # a) Check that provided SPAC has parametrized functions for IC and root density distribution
        #     # b) redefine resolution
        # end
        # remake = function(parameterisedSPAC::SPAC; Δz::NamedTuple)
        #     @assert keys(Δz) == (:thickness_m, :functions)
        # end
        # Δz = (thickness_m = [12 12 12],
        #       functions = (rootden = ((Δz_m)->LWFBrook90.Rootden_beta_(0.97, Δz_m = Δz_m)),
        #                    PSIM_init = ((Δz_m)->fill(-6.3, length(Δz_m))),
        #                    δ18Ο_init = ((Δz_m)->ifelse.(cumsum(Δz_m) .<= 0.2, -13., -10.)),
        #                    δ2Η_init  = ((Δz_m)->ifelse.(cumsum(Δz_m) .<= 0.2, -95., -70.))))
        # keys(Δz) == (:thickness_m, :functions)
        # # TODO

    # Solve ODE:
    simulate!(simulation);
    # plot(simulation.ODESolution)

    #u_soil1_amt_d18_d2 = reduce(hcat, [simulation.ODESolution[i].x[simulation.ODESolution.prob.p.row_idx_SWATI][1,:] for i = eachindex(simulation.ODESolution)])
    #u_soilX_amts       = reduce(hcat, [simulation.ODESolution[i].x[simulation.ODESolution.prob.p.row_idx_SWATI][:,1] for i = eachindex(simulation.ODESolution)])
    ####################

    ####################
    ## Benchmarking of input_path = "examples/isoBEAdense2010-18-reset-FALSE-input/"
    # @time simulation.ODESolution = solve_LWFB90(u0, tspan, p);
    # @time simulation.ODESolution = LWFBrook90.simulate!(simulation);
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
                    # LWFBrook90.f_LWFBrook90!(du, integrator.u, p, 0);
                    # du
                    # step!(integrator)

    # using BenchmarkTools # for benchmarking
    # simulation.ODESolution = @btime solve(ode_LWFBrook90; dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity
    # simulation.ODESolution = @btime solve(ode_LWFBrook90; saveat = tspan[1]:tspan[2], dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity
    ####################

    # ####################
    # # Postprocess simulation object (using provided functions)
    # timesteps = simulation.ODESolution_datetime
    # # Read out belowground variables
    # # get_δsoil # or alternatively: get_deltasoil
    # get_δsoil(simulation).d18O
    # get_δsoil(simulation).d18O
    # get_δsoil(simulation).d2H
    # get_δsoil(simulation; depths_to_read_out_mm = [100, 150]).d2H
    # get_δsoil(simulation; depths_to_read_out_mm = [100, 150], days_to_read_out_d = [1.0, 2.0]).d18O
    # # get_θ     # or alternatively: get_theta
    # theta_all_layers = get_θ(simulation); # returns a 2D matrix with soil layers as rows and time steps as columns (NOTE transpose to plot)
    # get_θ(simulation; depths_to_read_out_mm = [100, 150])
    # get_θ(simulation; depths_to_read_out_mm = [100, 150])
    # get_θ(simulation; depths_to_read_out_mm = [100, 150])
    # get_θ(simulation; depths_to_read_out_mm = [100, 150], days_to_read_out_d = [0.0, 1.0, 2.0])
    # get_θ(simulation; depths_to_read_out_mm = [100, 150], days_to_read_out_d = 1:1.0:100)
    # # get_ψ()     # or alternatively: get_psi
    # # get_W()
    # # get_SWATI()
    # # get_K()

    # plot(timesteps, theta_all_layers') # Note transpose to plot, so that each series to plot is a column

    # # Read out results for aboveground/scalar variables
    # timesteps_to_read_out = [0.0, 1.0, 2.0]
    # dates_to_read_out = RelativeDaysFloat2DateTime.(timesteps_to_read_out, simulation.continuous_SPAC.reference_date)
    # # LWFBrook90.DateTime2RelativeDaysFloat.( dates_to_read_out, simulation.continuous_SPAC.reference_date)
    # # LWFBrook90.DateTime2RelativeDaysFloat.( DateTime("2020-10-10", "yyyy-mm-dd"), simulation.continuous_SPAC.reference_date)
    # # LWFBrook90.DateTime2RelativeDaysFloat.( DateTime("2020-10-10_10:12", "yyyy-mm-dd_HH:MM"), simulation.continuous_SPAC.reference_date)
    # get_aboveground(simulation)
    # get_aboveground(simulation; days_to_read_out_d = timesteps_to_read_out)
    # dat_aboveground = get_aboveground(simulation; days_to_read_out_d = timesteps_to_read_out)
    # plot(timesteps_to_read_out, Matrix(dat_aboveground), labels = reshape(names(dat_aboveground),1,:))
    # plot(dates_to_read_out,     Matrix(dat_aboveground), labels = reshape(names(dat_aboveground),1,:))

    # dat_delta_aboveground = get_δ(simulation) # or alternatively: get_delta
    # plot(simulation.ODESolution_datetime, Matrix(dat_delta_aboveground[:,1:7]),
    #     labels = reshape(names(dat_delta_aboveground[:,1:7]),1,:))
    # plot(simulation.ODESolution_datetime, Matrix(dat_delta_aboveground[:,8:14]),
    #     labels = reshape(names(dat_delta_aboveground[:,8:14]),1,:))

    # # Postprocess simulation object (with your own functions)
    # # # Simulation object contains
    # # #   Type `simulation.` and wait for the autocomplete!
    # # simulation.soil_discretization
    # # simulation.continuous_SPAC.canopy_evolution
    # # simulation.continuous_SPAC.soil_horizons
    # # simulation.continuous_SPAC.params

    # # # simulation.ODESolution is a ODESolution object: documentation how to access under:
    # # # https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/
    # # # and more generally:
    # # # https://docs.sciml.ai/Overview/stable/

    # # simulation.ODESolution.t
    # # simulation.ODESolution.u
    # # simulation.ODESolution.retcode
    # # simulation.ODESolution.destats
    # ####################


    # ####################
    # # Note, that it is possible to use R Code from within Julia, e.g ggplot:
    # # https://stackoverflow.com/a/70073193/3915004
    # ####################

    ####################
    # Plotting simulation object (using provided functions)
    out_dir = "gitignored/out/$(today())"
    mkpath(out_dir)
    fname = joinpath(
        out_dir,
        input_prefix*"_NLAYER+" * string(simulation.ODESolution.prob.p.p_soil.NLAYER)*
        # "-git+"*chomp(Base.read(`git rev-parse --short HEAD`, String))*
        # ifelse(length(Base.read(`git status --porcelain`, String))==0, "+gitclean","+gitdirty")*
        "-iso+"*string(simulation.ODEProblem.p.simulate_isotopes))

    # 0a) Amounts
    pl1 = plotamounts(simulation, :above_and_belowground, :showRWUcentroid)
    # plotamounts(simulation, :θ) # errors as not yet implemented
    # plotamounts(simulation, :ψ) # errors as not yet implemented
    # plotamounts(simulation, :aboveground) # errors as not yet implemented
    # plotamounts(simulation, :mm) # errors as not yet implemented
    pl1
    savefig(pl1, fname*"_plotRecipe_AMT.png")

    # 0b) Isotopes
    pl2 = plotisotopes(simulation, :d18O, (d18O = :auto, d2H = :auto), :showRWUcentroid)
    pl2
    # plotisotopes(simulation, :d2H)
    # plotisotopes(simulation, :d18O)
    # plotisotopes(simulation;
    #              xlim = (DateTime("2010-01-01"), DateTime("2010-03-31")))
    # plotisotopes(simulation, :d2H, (d18O = :auto, d2H = :auto), :showRWUcentroid;
    #              xlim = (DateTime("2010-01-01"), DateTime("2010-08-30")))

    savefig(pl2, fname*"_plotRecipe_ISO.png")

    # 0c) Forcing and states: internal check:
    pl3 = plotforcingandstates(simulation)
    savefig(pl3, fname*"_plotRecipe_CHECK.png")
    # ####################

    ####################
    ## Plotting
    p = simulation.ODESolution.prob.p;
    θ_limits = (minimum(simulation.ODESolution.prob.p.p_soil.p_θr), maximum(simulation.ODESolution.prob.p.p_soil.p_THSAT));
    ψ_pF_limits = (-2, 7);

    # plot RWU time series
    # idx_u_scalar_isotopes_d18O = simulation.ODESolution.prob.p[1][4][8     ]
    #     idx_u_vector_isotopes_d18O = simulation.ODESolution.prob.p[1][4][9     ]
    #     idx_u_scalar_isotopes_d2H  = simulation.ODESolution.prob.p[1][4][10     ]
    #     idx_u_vector_isotopes_d2H  = simulation.ODESolution.prob.p[1][4][11     ]
    if simulate_isotopes
        # all states:                        [simulation.ODESolution[t] for t = eachindex(simulation.ODESolution)]
        # all states as matrix: reduce(hcat, [simulation.ODESolution[t] for t = eachindex(simulation.ODESolution)])
        # reduce(hcat, [simulation.ODESolution[t].x[1][:,1] for t = eachindex(simulation.ODESolution)])
        # reduce(hcat, [simulation.ODESolution[t].x[2][:,1] for t = eachindex(simulation.ODESolution)])
        # row_RWU_d18O = reduce(hcat, [simulation.ODESolution[t].x[simulation.ODESolution.prob.p.row_idx_scalars.RWU][:,simulation.ODESolution.prob.p.col_idx_d18O] for t = eachindex(simulation.ODESolution)])
        # row_XYL_d18O = reduce(hcat, [simulation.ODESolution[t].x[simulation.ODESolution.prob.p.row_idx_scalars.XYLEM][:,simulation.ODESolution.prob.p.col_idx_d18O] for t = eachindex(simulation.ODESolution)])
        row_RWU_d18O  = reduce(hcat, [simulation.ODESolution[t_idx].RWU.d18O   for t_idx = eachindex(simulation.ODESolution)])
        row_XYL_d18O  = reduce(hcat, [simulation.ODESolution[t_idx].XYLEM.d18O for t_idx = eachindex(simulation.ODESolution)])
        plot([transpose(row_RWU_d18O) transpose(row_XYL_d18O)], labels = ["δ_RWU" "δ_XylemV"])
    end

    if input_prefix == "BEA2016-reset-FALSE" || input_prefix == "isoBEA2016-reset-FALSE" || input_prefix == "isoBEAdense2016-reset-FALSE"
        ref_aboveground =
            DataFrame(File(
            # "test/test-assets-external/BEA-2016/BEA-IntegrationTests-LWFBrook90/output_LWFBrook90R/BEA2016-reset-FALSE_NLAYER14_LWFBrook90R-0.4.5daily_output.csv"))[:,
            "out/BEA2016-reset-FALSE_NLAYER14_LWFBrook90R-0.4.5daily_output.csv"))[:,
                [:yr, :mo, :da, :doy, :intr, :ints, :snow, :gwat]]
        pl_ab_3 = plot(simulation.ODESolution; vars = [(2 -1)*3+1, (3 -1)*3+1, (4 -1)*3+1], #vars = [2, 3, 4],
            label=["INTS (mm)" "INTR (mm)" "SNOW (mm)"])
        plot!(pl_ab_3,
                [ref_aboveground.intr,
                ref_aboveground.ints,
                ref_aboveground.snow], label = "LWFBrook90R", line = :dash, color = :black)
        savefig(fname*"_plot-INTS_INTR_SNOW.png")
        pl_ab_4 = plot(simulation.ODESolution; vars = [(2 -1)*3+1, (3 -1)*3+1],
            label=["INTS (mm)" "INTR (mm)"])
        plot!(pl_ab_4,
                [ref_aboveground.intr,
                ref_aboveground.ints], label = "LWFBrook90R", line = :dash, color = :black)
        savefig(fname*"_plot-INTS_INTR.png")
    end


    if (true)
        simulate_isotopes = true
        # theme(:default) # theme(:sand)
        PREC_color = :black
        depth_to_read_out_mm = [150, 500, 800, 1500]
        if simulate_isotopes
            df_δsoil      = get_soil_([:SWATI, :W, :ψ, :θ, :K, :δ18O, :δ2H], example_result; days_to_read_out_d = t_out, depths_to_read_out_mm = depth_to_read_out_mm)
            δ_resultsSoil = (d18O = Matrix(permutedims(select(df_δsoil, r"δ18O_"))), d2H = Matrix(permutedims(select(df_δsoil, r"δ2H_"))))
            δ_results = get_δ(simulation)
        end
        θ_limits = (minimum(simulation.ODESolution.prob.p.p_soil.p_θr), maximum(simulation.ODESolution.prob.p.p_soil.p_THSAT))
        ψ_pF_limits = (-2, 7)

        df_θψ = get_soil_([:θ, :ψ], simulation, depths_to_read_out_mm = depth_to_read_out_mm)

        pl_θ = plot(simulation.ODESolution_datetime,
            Matrix(select(df_θψ, r"θ_")),
            labels = permutedims(names(select(df_θψ, r"θ_"))),
            xlabel = "Date",
            ylabel = "θ\n[-]", ylims = θ_limits,
            legend = :outerright)
        pl_ψ = plot(simulation.ODESolution_datetime,
            # -Matrix(select(df_θψ, r"ψ_")) .+ 1, yaxis = :log, yflip = true,
            # Matrix(select(df_θψ, r"ψ_")),  ylabel = "ψ\n[kPa]",
            log10.(-10 .* Matrix(select(df_θψ, r"ψ_"))),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
            ylims = ψ_pF_limits,
            labels = permutedims(names(select(df_θψ, r"ψ_"))),
            xlabel = "Date",
            legend = :outerright);
        # if simulate_isotopes
            pl_δ18O = plot(simulation.ODESolution_datetime,
                Matrix(select(df_δsoil, r"δ18O_")),
                labels = permutedims(names(select(df_δsoil, r"δ18O_"))),
                xlabel = "Date",
                ylabel = "δ¹⁸O soil\n[‰]",
                legend = :outerright);
            pl_δ2H = plot(simulation.ODESolution_datetime,
                Matrix(select(df_δsoil, r"δ2H_")),
                labels = permutedims(names(select(df_δsoil, r"δ2H_"))),
                xlabel = "Date",
                ylabel = "δ²H soil\n[‰]",
                legend = :outerright);
            # add precipitation to soil δ
            plot!(pl_δ2H,
                simulation.ODESolution_datetime,δ_results[:,:PREC_d2H],
                labels = "PREC", color = PREC_color, linestyle = :dot);
            plot!(pl_δ18O,
                simulation.ODESolution_datetime,δ_results[:,:PREC_d18O],
                labels = "PREC", color = PREC_color, linestyle = :dot);
                # # add raw data values to check interpolation
                # data_days_to_plot = (input_meteoiso.days .> tspan[1] .&& input_meteoiso.days .< tspan[2])
                # scatter!(LWFBrook90.RelativeDaysFloat2DateTime.(input_meteoiso.days[data_days_to_plot], input_meteoveg_reference_date),
                #         input_meteoiso.delta18O_permil[data_days_to_plot])

        # else
        #     pl_δ18O = plot();
        #     pl_δ2H = plot();
        # end
        pl_PREC = plot(
            simulation.ODESolution_datetime,
            simulation.ODESolution.prob.p.p_PREC.(simulation.ODESolution.t),
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

    # BALERD_SWAT  = [simulation.ODESolution[idx].accum.BALERD_SWAT for idx in eachindex(simulation.ODESolution)]
    # BALERD_total = [simulation.ODESolution[idx].accum.BALERD_total for idx in eachindex(simulation.ODESolution)]

    # plot(simulation.ODESolution_datetime,
    #             [BALERD_SWAT BALERD_total],
    #             legend = :outerright, labels = ["BALERD_SWAT" "BALERD_total"],
    #             ylabel = "Water balance error [mm]")
    # savefig(fname*"_plot-water-balance-error.png")

    if (true)
        # Depth-vs-time heatmap plots of θ, ψ, RWU
        # For all heatmaps
        t_ref = simulation.parametrizedSPAC.reference_date
        x = simulation.ODESolution_datetime
        y = cumsum(simulation.ODESolution.prob.p.p_soil.p_THICK) - simulation.ODESolution.prob.p.p_soil.p_THICK/2
        optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
        y_soil_ticks = optim_ticks(0., round(maximum(cumsum(simulation.ODESolution.prob.p.p_soil.p_THICK))))[1]
        y_ticks= y_soil_ticks
        y_labels=round.(y_soil_ticks; digits=0)
        # a) for θ and ψ:
        (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
                        LWFBrook90.get_auxiliary_variables(simulation.ODESolution)
        t_to_plot = simulation.ODESolution.t
        pl_θ = heatmap(x, y,
                    u_aux_θ, colorbar_title = "θ [-]", clims = θ_limits, c=cgrad(:inferno, rev=true),
                    yflip = true,
                    xlabel = "",#"Date",
                    ylabel = "Depth [mm]");
        pl_ψ = heatmap(x, y,
                    # u_aux_PSIM', colorbar_title = "ψ_m [kPa]", clims = ψ_pF_limits, c=cgrad(:inferno, rev=false),
                    log10.(-10 .* u_aux_PSIM),  colorbar_title = "pF = log₁₀(-ψ hPa)", clims = ψ_pF_limits, c=cgrad(:inferno, rev=false),
                    yflip = true,
                    xlabel = "",#"Date",
                    ylabel = "Depth [mm]");
        # -LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm) .+ 1, yaxis = :log, yflip = true,
            # LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm),  ylabel = "ψ\n[kPa]",



        # plot(pl_θ,pl_ψ,layout = (2,1))
        # b) for RWU
        # b1) RWU heatmap (depths vs time)
        rows_RWU_mmDay = reduce(hcat, [simulation.ODESolution(t).TRANI.mmday   for t in t_to_plot])
        # plot(reduce(hcat, [simulation.ODESolution(t).RWU.d18O for t in t_to_plot])')
        # plot(reduce(hcat, [simulation.ODESolution(t).RWU.d2H for t in t_to_plot])')
        scalar_RWU_mmDay = reduce(hcat, [simulation.ODESolution(t).RWU.mmday   for t in t_to_plot])#simulation.ODESolution[simulation.ODESolution.prob.p.row_idx_RWU, 1, :]
        rows_RWU = rows_RWU_mmDay ./ simulation.ODESolution.prob.p.p_soil.p_THICK
        pl_RWU = heatmap(x,y,rows_RWU; #yticks = (y_ticks, y_labels),
                colorbar_title = "RWU [mm water/day per mm soil depth]",
                c = :diverging_bwr_20_95_c54_n256, clim = maximum(abs.(rows_RWU)) .* (-1, 1),
                yflip = true,
                ylabel = "Depth [mm]",
                size=(1400,800), dpi = 300, leftmargin = 15mm);
        # b2) distribution of RWU and roots
        plot_avgRWU = plot(sum(rows_RWU,dims=2), y,
            yflip = true, ylabel = "Depth [mm]", #xlabel = "Total RWU over\nsimulation period [mm]")
            #xlabel = "RWU (mm)"
            legend = :bottomright
        )
        t_toPlot = range(extrema(simulation.ODESolution.prob.tspan)..., length=5)
        root_data_start_end = [simulation.ODESolution.prob.p.p_fT_RELDEN.(t, 1:simulation.ODESolution.prob.p.p_soil.NLAYER) for t in t_toPlot]
        pl_roots = plot(root_data_start_end, y,
            labels = (format.(permutedims(RelativeDaysFloat2DateTime.(t_toPlot, t_ref)), "yyyy-mm")),
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
        # ode_LWFBrook90_2, unstable_check_function_2 = define_LWFB90_ODE(simulation.ODESolution.u[130], (130,140), p);
        # @time simulation.ODESolution_2 = solve(ode_LWFBrook90_2, Tsit5(); progress = true,
        #     unstable_check = unstable_check_function_2, # = (dt,u,p,t) -> false, #any(isnan,u),
        #     days_to_read_out_d = 130:140,
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
        #                 LWFBrook90.get_auxiliary_variables(simulation.ODESolution_2)
        # heatmap(RelativeDaysFloat2DateTime.(simulation.ODESolution_2.t, t_ref),
        #         cumsum(simulation.ODESolution_2.prob.p.p_soil.p_THICK) - simulation.ODESolution_2.prob.p.p_soil.p_THICK/2,
        #         u_aux_θ, yflip = true)
    end

    # # # 1) very basic
    # # using Plots # Install plot package at first use with `]` and then `add Plots`
    # # # plot(p.p_δ18O_PREC(1:300))
    # # # plot(p.p_δ2H_PREC(1:300))

    # # # # # Plot 1
    # # plot(simulation.ODESolution; vars = [1, 2, 3, 4, 5, 6],
    # #      label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

    # # idx_u_vector_amounts       = simulation.ODESolution.prob.p.row_idx_SWATI
    # # idx_u_scalar_isotopes_d18O = simulation.ODESolution.prob.p.row_idx_accum
    # # idx_u_vector_isotopes_d18O = simulation.ODESolution.prob.p[1][4][6  ]
    # # idx_u_scalar_isotopes_d2H  = simulation.ODESolution.prob.p.names_accum
    # # idx_u_vector_isotopes_d2H  = simulation.ODESolution.prob.p[1][4][8  ]
    # # idx_u_vector_accumulators  = simulation.ODESolution.prob.p[1][4][9  ]


    # # a = plot(simulation.ODESolution; vars = idx_u_scalar_isotopes_d18O,
    # #     ylabel = "δ18O [‰]",
    # #     label=["GWAT (‰)" "INTS (‰)" "INTR (‰)" "SNOW (‰)"])
    # # plot!(simulation.ODESolution.t, simulation.ODESolution.prob.p.p_δ18O_PREC.(simulation.ODESolution.t), label = "PREC (‰)")
    # # # plot!(ylims = (-50,50))

    # # plot(plot(simulation.ODESolution; vars = idx_u_scalar_isotopes_d18O[2], ylabel = "δ18O [‰]", label=["INTS (‰)"]),
    # #     plot(simulation.ODESolution; vars = [2],label=["INTS (mm)"]),
    # #     layout=(2,1))
    # # # plot!(simulation.ODESolution.t, simulation.ODESolution.prob.p.p_δ18O_PREC.(simulation.ODESolution.t), label = "PREC (‰)")
    # # plot(plot(simulation.ODESolution; vars = idx_u_scalar_isotopes_d18O[3], ylabel = "δ18O [‰]", label=["INTR (‰)"]),
    # #     plot(simulation.ODESolution; vars = [3],label=["INTR (mm)"]),
    # #     layout=(2,1))
    # # # plot!(simulation.ODESolution.t, simulation.ODESolution.prob.p.p_δ18O_PREC.(simulation.ODESolution.t), label = "PREC (‰)")
    # # plot(plot(simulation.ODESolution; vars = idx_u_scalar_isotopes_d18O[4], ylabel = "δ18O [‰]", label=["SNOW (‰)"]),
    # #     plot(simulation.ODESolution; vars = [4],label=["SNOW (mm)"]),#, ylims = (0,0.01)),
    # #     layout=(2,1))
    # # # plot!(simulation.ODESolution.t, simulation.ODESolution.prob.p.p_δ18O_PREC.(simulation.ODESolution.t), label = "PREC (‰)")
    # # plot(plot(simulation.ODESolution; vars = idx_u_scalar_isotopes_d18O[1], ylabel = "δ18O [‰]", label=["GWAT (‰)"]),
    # #     plot(simulation.ODESolution; vars = [1],label=["GWAT (mm)"]),
    # #     layout=(2,1))
    # # # plot!(simulation.ODESolution.t, simulation.ODESolution.prob.p.p_δ18O_PREC.(simulation.ODESolution.t), label = "PREC (‰)")


    # # # simulation.ODESolution[26,1,1:10]
    # # plot(simulation.ODESolution; vars = idx_u_scalar_isotopes_d2H,
    # #     ylabel = "δ2H [‰]",
    # #     label=["GWAT (‰)" "INTS (‰)" "INTR (‰)" "SNOW (‰)"])
    # # # plot!(simulation.ODESolution.t, simulation.ODESolution.prob.p.p_δ2H_PREC.(simulation.ODESolution.t), label = "PREC (‰)")
    # # scatter([
    # #         (simulation.ODESolution[13 + 1,:], simulation.ODESolution[24 + 1,:]),
    # #         (simulation.ODESolution[13 + 2,:], simulation.ODESolution[24 + 2,:]),
    # #         (simulation.ODESolution[13 + 3,:], simulation.ODESolution[24 + 3,:]),
    # #         (simulation.ODESolution[13 + 4,:], simulation.ODESolution[24 + 4,:])
    # #     ];
    # #     xlabel = "δ18O [‰]", ylabel = "δ2H [‰]",
    # #     label = ["GWAT" "INTS" "INTR" "SNOW"])
    # # # plot!(ylims = (-100,-60), xlims = (-15, -8))
    # # # plot!(simulation.ODESolution.prob.p.p_δ18O_PREC.(simulation.ODESolution.t), simulation.ODESolution.prob.p.p_δ2H_PREC.(simulation.ODESolution.t), label = "PREC")

    # # # # Plot 2
    # # # # http://docs.juliaplots.org/latest/generated/gr/#gr-ref43
    # # x = LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, input_meteoveg_reference_date)
    # # y = cumsum(simulation.ODESolution.prob.p.p_soil.p_THICK)
    # # n = simulation.ODESolution.prob.p.p_soil.NLAYER
    # # y_centers = [ 0; cumsum(simulation.ODESolution.prob.p.p_soil.p_THICK)[1:(n-1)] ] +
    # #     simulation.ODESolution.prob.p.p_soil.p_THICK / 2

    # # z = simulation.ODESolution[7 .+ (0:simulation.ODESolution.prob.p.p_soil.NLAYER-1),
    # #                     1,
    # #                     :]./simulation.ODESolution.prob.p.p_soil.p_THICK;
    # # z2 = simulation.ODESolution[idx_u_vector_isotopes_d18O,1,:];
    # # z3 = simulation.ODESolution[idx_u_vector_isotopes_d2H,1,:];

    # # heatmap(x, y_centers, z, yflip = true,
    # #         xlabel = "Date",
    # #         ylabel = "Depth [mm]",
    # #         colorbar_title = "θ [-]")
    # # hline!([0; cumsum(simulation.ODESolution.prob.p.p_soil.p_THICK)], yflip = true, xticks = false,
    # #     color = :black, linestyle = :dot
    # #     #title = "N_layer = "*string(simulation.ODESolution.prob.[1][1].NLAYER)
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
    # # # ###################

    # # plot(x,
    # #      z[LWFBrook90.find_soilDiscr_indices(simulation; depth_to_read_out_mm = depth_to_read_out_mm), :]',
    # #      labels = string.(permutedims(depth_to_read_out_mm)) .* "mm",
    # #      xlabel = "Date",
    # #      ylabel = "θ [-]",
    # #      legend = :bottomright)
    # # savefig(input_prefix*"_θ-feinUndSteinErde_depths_NLAYER"*string(simulation.ODESolution.prob.p.p_soil.NLAYER)*".png")





    # ###################################
    # ## More complex plots (comparing with LWFBrook90R results if available):
    # using Statistics
    # using Plots, Dates, StatsPlots, DataFramesMeta
    # using Gadfly, Cairo #,Fontconfig
    # using CSV

    # input_meteoveg = @chain DataFrame(File(path_meteoveg;
    #     skipto=3, delim=',', ignorerepeated=false,
    #     # Be strict about loading NA's -> error if NA present
    #     types=parsing_types, missingstring = nothing, strict=true))
    #     transform(:dates = DateTime.(:dates))
    #     end

    if (input_prefix == "DAV2010-2021" || input_prefix == "DAV_LW1_def");
        # File("../../../LWF-Brook90.jl-calibration/Calibration-Isotope-Data/2-combined-deposition/outputs/DAV2010-2021_meteoiso.csv")
        # File("../../../LWF-Brook90.jl-calibration/Input-Meteo-Data/DAV2010-2021/DAV2010-2021_meteoiso.csv")
        dat_δ_soilSol   = DataFrame(File("../../../LWF-Brook90.jl-calibration/Calibration-Isotope-Data/3-LWF-soil-solution/outputs/DAV-2022-07_soilsolutioniso.csv";header=1,skipto=3))
        dat_δ_soilbulk  = DataFrame(File("../../../LWF-Brook90.jl-calibration/Calibration-Isotope-Data/4-LWF-xylem-bulksoil/outputs/DAV-2022-07_bulksoiliso.csv";header=1,skipto=3))
        dat_δ_xylem     = DataFrame(File("../../../LWF-Brook90.jl-calibration/Calibration-Isotope-Data/4-LWF-xylem-bulksoil/outputs/DAV-2022-07_xylemiso.csv";header=1,skipto=3))
        dat_ψ_tensiometer= DataFrame(File("../../../LWF-Brook90.jl-calibration/Calibration-Soil-Data/1-LWF-tensiometer/DAV-2022-07_tensiometer.csv";header=1,skipto=3))
        dat_θ_EC5Waldner = DataFrame(File("../../../LWF-Brook90.jl-calibration/Calibration-Soil-Data/2-LWF-EC5-Waldner/outputs/DAV-2022-07_ec5-Waldner.csv";header=1,skipto=3))
        dat_ψ_tensiomark = DataFrame(File("../../../LWF-Brook90.jl-calibration/Calibration-Soil-Data/3-LWF-Tensiomark-EC5-Meusburger/outputs/DAV-2022-07_tensiomark-Meusburger.csv";header=1,skipto=3))
        dat_ψ_tensiomarkAVG= DataFrame(File("../../../LWF-Brook90.jl-calibration/Calibration-Soil-Data/3-LWF-Tensiomark-EC5-Meusburger/outputs/DAV-2022-07_tensiomarkAVG-Meusburger.csv";header=1,skipto=3))

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

        depth_to_read_out_mm = [150, 500, 800, 1500]

        @chain dat_ψ_tensiometer[1:100,:]   scatter(convert.(DateTime, _.dates), _.psim_median_kPa, group = _.depth_cm .* 10)
        pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, simulation.parametrizedSPAC.reference_date),
            # -LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm) .+ 1, yaxis = :log, yflip = true,
            LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm),  ylabel = "ψ\n[kPa]",
                         # TODO: replace get_θ(...) by get_soil_(:ψ, ...)
            # log10.(-10 .* LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm)),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
            labels = string.(permutedims(depth_to_read_out_mm)) .* "mm",
            xlabel = "Date",
            legend = :outerright)

        @chain dat_θ_EC5Waldner[1:900,:]    scatter(convert.(DateTime, _.dates), _.theta_m3m3,      group = _.depth_cm)
        df_θψ = get_soil_([:θ, :ψ], simulation, depths_to_read_out_mm = depth_to_read_out_mm)
        pl_θ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, simulation.parametrizedSPAC.reference_date),
            Matrix(select(df_θψ, r"θ_")),
            labels = permutedims(names(select(df_θψ, r"θ_"))),
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
        pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, simulation.parametrizedSPAC.reference_date),
            # -LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm) .+ 1, yaxis = :log, yflip = true,
            LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm)',  ylabel = "ψ\n[kPa]",
             # TODO: replace get_θ(...) by get_soil_(:ψ, ...)
            # log10.(-10 .* LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm)),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
            labels = reshape(string.(permutedims(depth_to_read_out_mm)) .* "mm",1,:),
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
        pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, simulation.parametrizedSPAC.reference_date),
            log10.(-10 .* LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm))',  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
            # TODO: replace get_θ(...) by get_soil_(:ψ, ...)
            labels = reshape(string.(permutedims(depth_to_read_out_mm)) .* "mm",1,:),
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

        pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, simulation.parametrizedSPAC.reference_date),
            # -LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm) .+ 1, yaxis = :log, yflip = true,
            LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm)',  ylabel = "ψ\n[kPa]",
             # TODO: replace get_θ(...) by get_soil_(:ψ, ...)
            # log10.(-10 .* LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm)),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
            labels = reshape(string.(permutedims(depth_to_read_out_mm)) .* "mm",1,:),
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
        pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, simulation.parametrizedSPAC.reference_date),
            log10.(-10 .* LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm))',  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
            # TODO: replace get_θ(...) by get_soil_(:ψ, ...)
            labels = reshape(string.(permutedims(depth_to_read_out_mm)) .* "mm",1,:),
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
        pl_ψ = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, simulation.parametrizedSPAC.reference_date),
            # log10.(-10 .* LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm)),  yflip = true, ylabel = "pF = \nlog₁₀(-ψ hPa)", #ylabel = "pF = \nlog₁₀(-ψₕₚₐ)",
            # TODO: replace get_θ(...) by get_soil_(:ψ, ...)
            LWFBrook90.get_ψ(simulation; depths_to_read_out_mm = depth_to_read_out_mm)',  ylabel = "ψ [kPa]",
            labels = reshape(string.(permutedims(depth_to_read_out_mm)) .* "mm",1,:),
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
        df_δsoil      = get_soil_([:δ18O, :δ2H], example_result; depths_to_read_out_mm = depth_to_read_out_mm  .+ 1)
        pl_δsoil = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, simulation.parametrizedSPAC.reference_date),
            Matrix(select(df_δsoil, r"δ18O_")),
            ylabel = "δ18O [permil]",
            labels = permutedims(names(select(df_δsoil, r"δ18O_"))),
            xlabel = "Date", legend = :outerright, color_palette = col_palette,
            linewidth = 3, alpha = 0.5)
        # add precipitation
        pl_δsoil = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, simulation.parametrizedSPAC.reference_date),
                        simulation.ODESolution.prob.p.p_δ18O_PREC.(simulation.ODESolution.t),
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
            LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, simulation.parametrizedSPAC.reference_date),
            permutedims([LWFBrook90.get_δ(simulation.ODESolution).RWU.d18O;
                         LWFBrook90.get_δ(simulation.ODESolution).XYL.d18O]),
            ylabel = "δ18O [permil]",
            labels = ["RWU" "Xylem"],
            xlabel = "Date", legend = :outerright,
            linewidth = 3, alpha = 0.5)
        # add precipitation
        pl_δsoil = plot!(LWFBrook90.RelativeDaysFloat2DateTime.(simulation.ODESolution.t, simulation.parametrizedSPAC.reference_date),
                        simulation.ODESolution.prob.p.p_δ18O_PREC.(simulation.ODESolution.t),
                        color = "grey", alpha = 0.5, labels = "Precipiation")

        # ylims!(-15,0)
        xlims!(extrema(convert.(DateTime, dat_δ_xylem.dates)) .+ (-Day(100), Day(50)))
        plot!(size=(800,500), dpi = 200, leftmargin = 5mm)
        savefig(fname*"_vs-δ-xylem.png")
    end
end



# input_prefix = "isoBEA2016-reset-FALSE"
# input_path = "examples/isoBEA2016-reset-FALSE-input/"
# input_prefix = "isoBEA2010-18-reset-FALSE";
# input_path = "examples/isoBEA2010-18-reset-FALSE-input/";
input_prefix = "isoBEAdense2010-18-reset-FALSE";
input_path = "examples/isoBEAdense2010-18-reset-FALSE-input/";
run_main_with_isotopes(;input_prefix = input_prefix, input_path = input_path)
input_prefix = "BEA2010-2021";
input_path = "../../../LWF-Brook90.jl-calibration/Input-Meteo-Data/BEA2010-2021/";
run_main_with_isotopes(;input_prefix = input_prefix, input_path = input_path)
# input_prefix = "LAU2010-2021";
# input_path = "../../../LWF-Brook90.jl-calibration/Input-Meteo-Data/LAU2010-2021/";
# input_prefix = "LAE2010-2021";
# input_path = "../../../LWF-Brook90.jl-calibration/Input-Meteo-Data/LAE2010-2021/";
# input_prefix = "DAV2010-2021";
# input_path = "../../../LWF-Brook90.jl-calibration/Input-Meteo-Data/DAV2010-2021/";
# input_prefix = "WaldLab";
# input_path = "../../../LWF-Brook90.jl-calibration/Input-Meteo-Data/WaldLab/";
input_prefix = "DAV_LW1_def";
input_path = "../../../LWF-Brook90.jl-calibration/Input-Meteo-Data/DAV_LW1_def/";
run_main_with_isotopes(;input_prefix = input_prefix, input_path = input_path)
input_prefix = "SCH-2023-03";
input_path = "../../../LWF-Brook90.jl-calibration/pgm-calibrate-HYPERION/input/";
run_main_with_isotopes(;input_prefix = input_prefix, input_path = input_path)

input_prefix = "infiltrationSaturationINFEXP0";
input_path = "examples/infiltrationSaturationINFEXP0/";
run_main_with_isotopes(;input_prefix = input_prefix, input_path = input_path)
input_prefix = "infiltrationSaturationINFEXP1";
input_path = "examples/infiltrationSaturationINFEXP1/";
run_main_with_isotopes(;input_prefix = input_prefix, input_path = input_path)






####### Test case for cumulative daily water balance error (Ireson 2023)
using LWFBrook90
input_path = "examples/DAV2020-full/"; input_prefix = "DAV2020-full";
model = loadSPAC(input_path, input_prefix; simulate_isotopes = true);
simulation          = setup(model)
simulate!(simulation, saveat = range(simulation.ODEProblem.tspan..., step = 30), save_everystep = false)
simulate!(simulation)
using Plots, Measures; gr();
pl1 = plotamounts(simulation, :above_and_belowground, :showRWUcentroid)
pl1

t_ref = simulation.ODESolution.prob.p.REFERENCE_DATE
# t_plot = simulation.ODESolution.t
t_plot = range(extrema(simulation.ODESolution.t)..., step=1.0) # we want daily values exactly!
x_plot = RelativeDaysFloat2DateTime.(t_plot, t_ref);

plot(pl1[4])
plot!(x_plot, [simulation.ODESolution(t).accum.StorageSWAT  for t in t_plot])
plot!(x_plot, [simulation.ODESolution(t).accum.StorageWATER for t in t_plot])
# plot!(x_plot, [simulation.ODESolution(t).accum.ε_prev_StorageSWAT  for t in t_plot])
# plot!(x_plot, [simulation.ODESolution(t).accum.ε_prev_StorageWATER for t in t_plot])

[simulation.ODESolution(t).accum.BALERD_SWAT for t in t_plot]

plot(pl1[5])

plot!(x_plot, [simulation.ODESolution(t).accum.BALERD_SWAT         for t in t_plot])
plot!(x_plot, cumsum([simulation.ODESolution(t).accum.BALERD_SWAT  for t in t_plot]))
plot!(x_plot, [simulation.ODESolution(t).accum.BALERD_total        for t in t_plot])
plot!(x_plot, cumsum([simulation.ODESolution(t).accum.BALERD_total for t in t_plot]))
