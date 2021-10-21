using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5
# example = LWFBrook90.run_example()

# Read in input data
# input_prefix = "isoBEA2016-reset-FALSE"
# input_path = "example/BEA2016-reset-FALSE-input/"
input_prefix = "isoBEA2010-18-reset-FALSE";
input_path = "example/isoBEA2010-18-reset-FALSE-input/";

####################
simulate_isotopes = true
(input_meteoveg,
    input_meteoiso,
    input_meteoveg_reference_date,
    input_param,
    input_storm_durations,
    input_initial_conditions,
    input_soil_horizons,
    input_soil_discretization,
    simOption_FLAG_MualVanGen) = read_inputData(input_path, input_prefix;
                                                simulate_isotopes = simulate_isotopes);
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
####################

####################
# Define parameters for differential equation
uSoil_initial, p = define_LWFB90_p(
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
    simulate_isotopes = simulate_isotopes#,
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
#TODO:simulate_isotopes, input_soil_discretization
u0 = define_LWFB90_u0(p, input_initial_conditions,
                      uSoil_initial,
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
tspan = (0., 700.)

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
@time sol_LWFBrook90 = solve(ode_LWFBrook90, Tsit5(); progress = true,
    unstable_check = unstable_check_function, # = (dt,u,p,t) -> false, #any(isnan,u),
    saveat = tspan[1]:tspan[2], dt=1e-3, adaptive = false);
    # 700 days in: 40 seconds dt=1e-3, adaptive = false (isoBea default spacing: NLAYER = 7)
    # 700 days in: 90 seconds dt=1e-3, adaptive = false (isoBea dense spacing (0.05m): NLAYER = 26)
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

using Plots, Measures
using LWFBrook90
optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
# pl1 = LWFBrook90.ISO.plotisotopes(sol_LWFBrook90);
pl2 = LWFBrook90.ISO.plotisotopes(
    sol_LWFBrook90, optim_ticks;
    layout = grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4]),
    size=(1000,1400), dpi = 300, leftmargin = 15mm);
plot!(pl2, link = :x);
savefig(pl2, input_prefix*"_plotRecipe.png")

# # 1) very basic
# using Plots # Install plot package at first use with `]` and then `add Plots`
# # plot(p[2][15](1:300))
# # plot(p[2][16](1:300))

# # # # Plot 1
# plot(sol_LWFBrook90; vars = [1, 2, 3, 4, 5, 6],
#      label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

# idx_u_vector_amounts       = sol_LWFBrook90.prob.p[1][4][4]
# idx_u_scalar_isotopes_d18O = sol_LWFBrook90.prob.p[1][4][5]
# idx_u_vector_isotopes_d18O = sol_LWFBrook90.prob.p[1][4][6]
# idx_u_scalar_isotopes_d2H  = sol_LWFBrook90.prob.p[1][4][7]
# idx_u_vector_isotopes_d2H  = sol_LWFBrook90.prob.p[1][4][8]
# idx_u_vector_accumulators  = sol_LWFBrook90.prob.p[1][4][9]


# a = plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d18O,
#     ylabel = "δ18O [mUr]",
#     label=["GWAT (mUr)" "INTS (mUr)" "INTR (mUr)" "SNOW (mUr)"])
# plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), label = "PREC (mUr)")
# # plot!(ylims = (-50,50))

# plot(plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d18O[2], ylabel = "δ18O [mUr]", label=["INTS (mUr)"]),
#     plot(sol_LWFBrook90; vars = [2],label=["INTS (mm)"]),
#     layout=(2,1))
# # plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), label = "PREC (mUr)")
# plot(plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d18O[3], ylabel = "δ18O [mUr]", label=["INTR (mUr)"]),
#     plot(sol_LWFBrook90; vars = [3],label=["INTR (mm)"]),
#     layout=(2,1))
# # plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), label = "PREC (mUr)")
# plot(plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d18O[4], ylabel = "δ18O [mUr]", label=["SNOW (mUr)"]),
#     plot(sol_LWFBrook90; vars = [4],label=["SNOW (mm)"]),#, ylims = (0,0.01)),
#     layout=(2,1))
# # plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), label = "PREC (mUr)")
# plot(plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d18O[1], ylabel = "δ18O [mUr]", label=["GWAT (mUr)"]),
#     plot(sol_LWFBrook90; vars = [1],label=["GWAT (mm)"]),
#     layout=(2,1))
# # plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), label = "PREC (mUr)")


# # sol_LWFBrook90[26,1,1:10]
# plot(sol_LWFBrook90; vars = idx_u_scalar_isotopes_d2H,
#     ylabel = "δ2H [mUr]",
#     label=["GWAT (mUr)" "INTS (mUr)" "INTR (mUr)" "SNOW (mUr)"])
# # plot!(sol_LWFBrook90.t, sol_LWFBrook90.prob.p[2][16].(sol_LWFBrook90.t), label = "PREC (mUr)")
# scatter([
#         (sol_LWFBrook90[13 + 1,:], sol_LWFBrook90[24 + 1,:]),
#         (sol_LWFBrook90[13 + 2,:], sol_LWFBrook90[24 + 2,:]),
#         (sol_LWFBrook90[13 + 3,:], sol_LWFBrook90[24 + 3,:]),
#         (sol_LWFBrook90[13 + 4,:], sol_LWFBrook90[24 + 4,:])
#     ];
#     xlabel = "δ18O [mUr]", ylabel = "δ2H [mUr]",
#     label = ["GWAT" "INTS" "INTR" "SNOW"])
# # plot!(ylims = (-100,-60), xlims = (-15, -8))
# # plot!(sol_LWFBrook90.prob.p[2][15].(sol_LWFBrook90.t), sol_LWFBrook90.prob.p[2][16].(sol_LWFBrook90.t), label = "PREC")

# # # Plot 2
# # # http://docs.juliaplots.org/latest/generated/gr/#gr-ref43
# x = LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date)
# y = cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK)
# z = sol_LWFBrook90[7 .+ (0:sol_LWFBrook90.prob.p[1][1].NLAYER-1),
#                     1,
#                     :]./sol_LWFBrook90.prob.p[1][1].p_THICK;
# z2 = sol_LWFBrook90[idx_u_vector_isotopes_d18O,1,:];
# z3 = sol_LWFBrook90[idx_u_vector_isotopes_d2H,1,:];
# heatmap(x, y, z, yflip = true,
#         xlabel = "Date",
#         ylabel = "Depth [mm]",
#         colorbar_title = "θ [-]")
# heatmap(x, y, z2, yflip = true,
#         xlabel = "Date",
#         ylabel = "Depth [mm]",
#         colorbar_title = "δ18O [mUr]")
# heatmap(x, y, z3, yflip = true,
#         xlabel = "Date",
#         ylabel = "Depth [mm]",
#         colorbar_title = "δ2H [mUr]")
# ####################
