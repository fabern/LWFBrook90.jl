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
using Test: @testset, @test, @test_throws, @test_broken, @test_skip, @test_logs
using Random
using Printf
using Logging
using Dates: today

# Note the git hash-string and repository status (can be optionally used in filenames for plots etc.)
git_status_string = "$(today())-git+"*
            "__"*chomp(Base.read(`git branch --show-current`, String))*
            "__"     *chomp(Base.read(`git rev-parse --short HEAD`, String))*
            ifelse(length(Base.read(`git status --porcelain`, String))==0, "+gitclean","+gitdirty")*
            "__"

# A macro for timing that also prints out the git commit hash:
macro githash_time(variable)
    quote
        #model = replace(chomp(Base.read(`sysctl hw.model`, String)), "hw.model: " => "")
        model = "amberMBP"
        hash = chomp(Base.read(`git rev-parse --short HEAD`, String))
        print(model*"-git-"*hash*":")
        @time $(esc(variable))
    end
end
# Alternatively: think about measuring performance along the lines discussed in
# https://discourse.julialang.org/t/benchmarking-tests-to-ensure-prs-dont-introduce-regressions/8630/6
# or then https://github.com/maxbennedich/julia-regression-analysis or ...)

# A flag that determines if tests are run on a CI system
is_a_CI_system = issubset(["GITHUB_ACTION"], collect(keys(ENV))) # checks if ENV["GITHUB_ACTION"] exists
@show is_a_CI_system


Random.seed!(1234)

# cd("test")

# if !is_a_CI_system; include("00-plot-if-not-CI-system.jl"); end

        @show pwd()
        parametrizedSPAC = loadSPAC("examples/DAV2020-full/", "DAV2020-full"; simulate_isotopes = true);
        Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 5)]; # grid spacing (heterogenous), meter (N=20)
        parametrizedSPAC = loadSPAC("examples/DAV2020-full/", "DAV2020-full";
            simulate_isotopes = true,
            Δz_thickness_m = Δz_m,
            root_distribution = (beta = 0.77, z_rootMax_m = -0.5),
            IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -10.11111, delta2H_init_permil = -91.1111),
        );

        soil_output_depths_m = zeros(Float64, 0)
        parametrizedSPAC.tspan
        DS = LWFBrook90.setup(parametrizedSPAC; requested_tspan = (0., 10.));

        # Test the f and cb() functions (LWFBrook90.f_LWFBrook90!, )
        t0 = 0.0
        u0 = DS.ODEProblem.u0
        p = DS.ODEProblem.p
        du = copy(u0)

        LWFBrook90.f_LWFBrook90!(du, u0, p, 1.0)
        integrator = init(DS.ODEProblem, Tsit5();)
        # @time LWFBrook90.LWFBrook90R_updateAmounts_INTS_INTR_SNOW_CC_SNOWLQ!(integrator)
        @time LWFBrook90.LWFBrook90R_updateIsotopes_INTS_INTR_SNOW!(integrator)
        # # LWFBrook90.LWFBrook90R_updateIsotopes_GWAT_SWAT!(u0, t0, integrator)
        @btime LWFBrook90.LWFBrook90R_updateIsotopes_GWAT_SWAT_AdvecDiff!(u0, t0, integrator)
        # # LWFBrook90.LWFBrook90R_check_balance_errors!(integrator)
        # # # @run simulate!(DS) # with ArrayParition this simulating 10 days takes about 0.373707 seconds (1.45 M allocations: 105.396 MiB, 14.36% gc time)
        # DS.ODEProblem.tspan
        # simulate!(DS); # with ArrayParition this simulating 10 days takes about 0.373707 seconds (1.45 M allocations: 105.396 MiB, 14.36% gc time) Time steps for solving: 580 (580 accepted out of 612 total)
        # # simulate!(DS); # with ArrayParition this simulating 10 days takes about 0.373707 seconds (1.45 M allocations: 105.396 MiB, 14.36% gc time) Time steps for solving: 580 (580 accepted out of 612 total)
        # DS = LWFBrook90.setup(parametrizedSPAC::SPAC; Δz = Δz);
        # DS.ODEProblem.tspan
        # # simulate!(DS); # with ArrayParition this simulating all the days takes about 167.115807 seconds (762.30 M allocations: 53.365 GiB, 6.18% gc time, 12.43% compilation time) Time steps for solving: 295894 (295894 accepted out of 314371 total)





    # TRY OUT WITH CALIBRATION RUNS:
        # julia --project=. --threads=auto -- script7-calibration.jl "LAU" "with-isotopes" 10

        site = "LAU"
        simulate_isotopes = true
        N_forwardRuns = 10

        parametrizedSPAC = loadSPAC("examples/DAV2020-full/", "DAV2020-full"; simulate_isotopes = true);


input_path            = "../../../LWF-Brook90.jl-calibration/pgm-calibrate-HYPERION/input-2023-07";
calibration_data_path = "../../../LWF-Brook90.jl-calibration/pgm-calibrate-HYPERION/input-2023-07";

input_prefix = (
    site == "SCH" ? "SCH-2023-07-input" :
    site == "LAU" ? "LAU-2023-07-input" :
    site == "CEL" ? "CEL-2023-07-input" :
    # error if unsupported site:
    error("Unsupported site requested. Requested: '$site'."))
@show input_prefix
(fname_ψ_TMBoWa,
fname_θ_EC5BoWa,
fname_δ_xylem,
fname_δ_soilSol,
fname_δ_soilbulk) = (
    site == "SCH" ? ("SCH-2023-07-calibrate_BoWaTensiomarkAVG.csv",
                        "SCH-2023-07-calibrate_BoWaEC5AVG.csv",
                        "SCH-2023-07-calibrate_XylemIso.csv",
                        "SCH-2023-07-calibrate_SoilSolutionIso.csv",
                        "SCH-2023-07-calibrate_BulkSoilIso.csv") :
    site == "LAU" ? ("LAU-2023-07-calibrate_BoWaTensiomarkAVG.csv",
                        "LAU-2023-07-calibrate_BoWaEC5AVG.csv",
                        "LAU-2023-07-calibrate_XylemIso.csv",
                        "LAU-2023-07-calibrate_SoilSolutionIso.csv",
                        "LAU-2023-07-calibrate_BulkSoilIso.csv") :
    site == "CEL" ? ("CEL-2023-07-calibrate_BoWaTensiomarkAVG.csv",
                        "CEL-2023-07-calibrate_BoWaEC5AVG.csv",
                        "CEL-2023-07-calibrate_XylemIso.csv",
                        "CEL-2023-07-calibrate_SoilSolutionIso.csv",
                        "CEL-2023-07-calibrate_BulkSoilIso.csv") :
    # error if unsupported site:
    error("Requested site: $site. But no calibration data configured."))



############ Load data
# include("../../../LWF-Brook90.jl-calibration/pgm-calibrate-HYPERION/helper_load-data-script7.jl")
include("../../../../LWF-Brook90.jl-calibration/pgm-calibrate-HYPERION/helper_load-data-script7.jl")

# prepare observed values
obs_combined_df =
    load_calibration_data(
        calibration_data_path,
        fname_ψ_TMBoWa,
        fname_θ_EC5BoWa,
        fname_δ_xylem,
        fname_δ_soilSol,
        fname_δ_soilbulk)
# jldsave(joinpath(out_dir, fname*"_reproducibility-00-obs_combined_df.jld2"), true, IOStream; obs_combined_df)

            # for development use:
            obs_summary = @chain obs_combined_df begin
                groupby([:site, :variable, :depth, :species, :treeID])
                @combine begin
                    :N_obs_available = length(:obs)
                    :min_Date = minimum(:dates)
                    :max_Date = maximum(:dates)
                end
            end
            @show obs_summary;

            @orderby obs_summary :variable :depth

# Derive simulation span and when to store simulation values
minmax_dates = extrema(@subset(obs_combined_df, :variable .== "θ_m3m3" .|| :variable .== "ψ_kPa").dates)

# Simulate maximally to 2022-12-31
minmax_dates = (minmax_dates[1], min(minmax_dates[2], DateTime("2022-12-31")))
@show minmax_dates

# Derive simulation tspan
spinup_days = 365
simulation_tspan_dates = minmax_dates .- Day.([spinup_days, 0])

@show spinup_days
@show simulation_tspan_dates




####### Setup base simulation ##############################################################
# Define simulation model by reading in system definition and input data
# manually specify resolution of simulation:
Δz = (
    site == "SCH" ? #[fill(0.02, 10); fill(0.10, 16)] :
                    # [fill(0.04, 5);  fill(0.10, 16)] : # 47 (8.3) seconds for 1000 days (1862. 2860.) with isotopes, 9.2 (XX) s without (in brackets v0.9.6)
                    [fill(0.05, 4);  fill(0.10, 14)] : # 18 (4.9) seconds for 1000 days (1862. 2860.) with isotopes, 4.1 (2.5) s without (in brackets v0.9.6)
    site == "LAU" ? [fill(0.05, 4);  fill(0.10, 14)] :
    site == "CEL" ? [fill(0.05, 4);  fill(0.10, 10)] :
    # error if unsupported site:
    error("Unspecified Δz for requested site. Requested: '$site'."))

# generate base model (by specifiying options here, less input files are needed):
loadSPAC_args = (input_path = input_path, input_prefix = input_prefix,
        simulate_isotopes = simulate_isotopes,
        Δz_thickness_m    = Δz,
        root_distribution = (beta = 0.98, z_rootMax_m=-0.9),
        IC_soil           = (PSIM_init_kPa = -6.0,
                            delta18O_init_permil = -9.0,
                            delta2H_init_permil = -11.0),
        canopy_evolution  = (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel    = 100,
                             LAI_rel = (DOY_Bstart = 110,    Bduration  = 20,
                                        DOY_Cstart = 270,    Cduration  = 60,
                                        LAI_perc_BtoC = 100, LAI_perc_CtoB = 0)),
        storm_durations_h = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
        IC_scalar         = (amount = (u_GWAT_init_mm = 0.0,         u_INTS_init_mm = 0.0,
                                       u_INTR_init_mm = 0.0,         u_SNOW_init_mm = 0.0,
                                       u_CC_init_MJ_per_m2 = 0.0001, u_SNOWLQ_init_mm =  0.),
                            d18O    = (u_GWAT_init_permil = -11.111, u_INTS_init_permil = -12.222,
                                       u_INTR_init_permil = -13.333, u_SNOW_init_permil = -14.444),
                            d2H     = (u_GWAT_init_permil = -95.111, u_INTS_init_permil = -95.222,
                                       u_INTR_init_permil = -95.333, u_SNOW_init_permil = -95.444)
));
base_model = loadSPAC(
    loadSPAC_args[:input_path],
    loadSPAC_args[:input_prefix];
    loadSPAC_args[Not([:input_path, :input_prefix])]...)

LWFBrook90.DateTime2RelativeDaysFloat.(simulation_tspan_dates, base_model.reference_date)
base_model_tspan_dates = LWFBrook90.RelativeDaysFloat2DateTime.(base_model.tspan, base_model.reference_date)
simulation_tspan = LWFBrook90.DateTime2RelativeDaysFloat.(simulation_tspan_dates, base_model.reference_date)

# does_x_cover_y(base_model_tspan_dates, simulation_tspan_dates) || @warn "Available period of input data $base_model_tspan_dates, does not cover requested simulation period $simulation_tspan_dates."

# Note, transformed tspan might not start at 0, but that's ok.
@show  simulation_tspan_dates

# soil_output_depths_m = -[0.150, 0.500, 0.800, 1.500]
base_simulation = setup(base_model, requested_tspan = simulation_tspan_dates);
display(base_simulation)


# Solve ODE:
@time simulate!(base_simulation);
@time simulate!(base_simulation);

# plot(simulation.ODESolution)
plotamounts(base_simulation)
pl2 = plotisotopes(base_simulation, :d18O, (d18O = :auto, d2H = :auto), :showRWUcentroid)
pl2

u = base_simulation.ODESolution.u
Plots.plot([ut.INTS.d18O for ut in u], label = "INTS")
Plots.plot!([ut.INTR.d18O for ut in u], label = "INTR")
Plots.plot!([ut.SNOW.d18O for ut in u], label = "SNOW")

# Test the f and cb() functions (LWFBrook90.f_LWFBrook90!, )
using BenchmarkTools
t0 = 0.0
u0 = base_simulation.ODEProblem.u0
p = base_simulation.ODEProblem.p
du = copy(u0)
LWFBrook90.f_LWFBrook90!(du, u0, p, 1.0)
integrator = init(base_simulation.ODEProblem, Tsit5();)
@time LWFBrook90.LWFBrook90R_updateIsotopes_INTS_INTR_SNOW!(integrator)
@btime LWFBrook90.LWFBrook90R_updateIsotopes_GWAT_SWAT_AdvecDiff!(u0, t0, integrator)
# @enter LWFBrook90.LWFBrook90R_updateIsotopes_GWAT_SWAT_AdvecDiff!(u0, t0, integrator)

@time LWFBrook90.LWFBrook90R_updateIsotopes_INTS_INTR_SNOW!(integrator)
@time LWFBrook90.LWFBrook90R_updateIsotopes_INTS_INTR_SNOW!(integrator)


@btime LWFBrook90.ISO.compute_isotope_U_of_INTS_INTR_SNOW_and_SLFL(
                $integrator.p.p_δ2H_PREC($integrator.t), $integrator.p.p_δ18O_PREC($integrator.t), $integrator.p.p_fT_TADTM[1], $integrator.p.p_VAPPRES($integrator.t),
                # for INTS (in: SINT; out: ISVP):
                $integrator.u.INTS.mm, $integrator.p.aux_du_SINT[1], $integrator.p.aux_du_ISVP[1], $integrator.p.p_DTP, $integrator.u.INTS.d2H, $integrator.u.INTS.d18O,
                # for INTR (in: RINT; out: IRVP):
                $integrator.u.INTR.mm, $integrator.p.aux_du_RINT[1], $integrator.p.aux_du_IRVP[1], $integrator.u.INTR.d2H, $integrator.u.INTR.d18O,
                # for SNOW (in: STHR, RSNO (both δ_PREC); out: SMLT, SNVP (δ_SNOW and fractionated)):
                $integrator.u.SNOW.mm, $integrator.p.p_fu_STHR[1], $integrator.p.aux_du_RSNO[1], $integrator.p.aux_du_SMLT[1], $integrator.p.aux_du_SNVP[1], $integrator.u.SNOW.d2H, $integrator.u.SNOW.d18O,
                # to compute isotopic signature of soil infiltration: SLFL
                $integrator.p.p_fu_RNET[1])
