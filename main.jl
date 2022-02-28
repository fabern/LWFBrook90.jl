using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5
# example = LWFBrook90.run_example()

# Read in input data
input_prefix = "BEA2016-reset-FALSE"
input_path = "example/" * input_prefix * "-input/"

####################
(input_meteoveg,
    input_meteoveg_reference_date,
    input_param,
    input_storm_durations,
    input_initial_conditions,
    input_soil_horizons,
    simOption_FLAG_MualVanGen) = read_inputData(input_path, input_prefix)

input_soil_discretization = discretize_soil(input_path, input_prefix)

####################

####################
# Define solver options
Reset = false                          # currently only Reset = 0 implemented
compute_intermediate_quantities = true # Flag whether ODE containes additional quantities than only states

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
    input_meteoveg_reference_date,
    input_param,
    input_storm_durations,
    input_soil_horizons,
    input_soil_discretization,
    simOption_FLAG_MualVanGen;
    Reset = Reset,
    # soil_output_depths = collect(-0.05:-0.05:-1.1),
    compute_intermediate_quantities = compute_intermediate_quantities)
####################

####################
# Define initial states of differential equation
# state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI
# Create u0 for DiffEq.jl
u0 = define_LWFB90_u0(p, input_initial_conditions,
    uSoil_initial,
    compute_intermediate_quantities)
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
tspan = (minimum(input_meteoveg[:, "days"]),
    maximum(input_meteoveg[:, "days"])) # simulate all available days
# tspan = (LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1980,1,1), reference_date),
#          LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1985,1,1), reference_date)) # simulates selected period

# Define ODE:
ode_LWFBrook90 = define_LWFB90_ODE(u0, tspan, p)
####################

####################
## Solve ODE:
sol_LWFBrook90 = solve(ode_LWFBrook90, Tsit5();
    progress = true,
    saveat = tspan[1]:tspan[2], dt = 1e-6, adaptive = true); # dt is initial dt, but adaptive
####################

####################
## Benchmarking
# @time sol_LWFBrook90 = solve(ode_LWFBrook90, progress = true;
#     saveat = tspan[1]:tspan[2], dt=1e-1, adaptive = false);
# @time sol_LWFBrook90 = solve(ode_LWFBrook90, progress = true, Euler(); # Note: Euler sometimes hangs
#     saveat = tspan[1]:tspan[2], dt=1e-1, adaptive = false);
# using BenchmarkTools # for benchmarking
# sol_LWFBrook90 = @btime solve(ode_LWFBrook90; dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity
# sol_LWFBrook90 = @btime solve(ode_LWFBrook90; saveat = tspan[1]:tspan[2], dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity
####################


####################
## Plotting
mkpath("out")
# 0) using plotting recipe defined in LWFBrook90.jl
using Plots, Measures
using LWFBrook90
optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)

pl_final = LWFBrook90.plotlwfbrook90(sol_LWFBrook90, optim_ticks)
savefig(pl_final,
    joinpath("out", input_prefix * "_plotRecipe_NLAYER" * string(sol_LWFBrook90.prob.p[1][1].NLAYER) * ".png")
    )

# 1) very basic
# using Plots # Install plot package at first use with `]` and then `add Plots`
# # Plot 1
# plot(sol_LWFBrook90; vars = [1, 2, 3, 4, 5, 6],label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])
# # Plot 2
# # http://docs.juliaplots.org/latest/generated/gr/#gr-ref43
# x = LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date)
# y = cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK)
# z = sol_LWFBrook90[7 .+ (0:sol_LWFBrook90.prob.p[1][1].NLAYER-1), 1, :]./sol_LWFBrook90.prob.p[1][1].p_THICK
# heatmap(x, y, z, yflip = true,
#         xlabel = "Date",
#         ylabel = "Depth [mm]",
#         colorbar_title = "θ [-]")
####################
