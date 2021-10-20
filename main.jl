# This is an example script. See others in folder "examples/scripts/".

using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5, init
# example = LWFBrook90.run_example()

# Read in input data
input_prefix = "BEA2016-reset-FALSE"
input_path = "examples/" * input_prefix * "-input/"

####################
(input_meteoveg,
    _input_meteoiso, # unused
    input_meteoveg_reference_date,
    input_param,
    input_storm_durations,
    input_initial_conditions,
    input_soil_horizons,
    simOption_FLAG_MualVanGen) = read_inputData(input_path, input_prefix)

# Define grid for spatial discretization as well as initial conditions and root densities
# a) either read the discretization from a file `soil_discretization.csv`
unused = discretize_soil(input_path, input_prefix)

# b) or define them manually
ψ_init = unused.uAux_PSIM_init_kPa[1]
# Δz_m = round.(-diff([unused.Upper_m[1]; unused.Lower_m]), digits=5)
# Δz_m = [fill(0.005, 4); fill(0.02, 99)]                          # grid spacing (heterogenous), meter (N=103)
# Δz_m = round.(diff(0.0:0.05:-minimum(unused.Lower_m)), digits=5) # grid spacing (heterogenous), meter

# As a test: subsequently increase resolution of the top layers.
Δz_m = [0.04, 0.04, 0.12, 0.25, 0.3, 0.35, 0.1]                            # grid spacing (heterogenous), meter (N=7)
# Δz_m = [0.04, 0.04, 0.04, 0.08, 0.25, 0.3, 0.35, 0.1]                    # grid spacing (heterogenous), meter (N=8)
# Δz_m = [0.04, 0.04, 0.04, 0.04, 0.04, 0.25, 0.3, 0.35, 0.1]              # grid spacing (heterogenous), meter (N=9)
# Δz_m = [fill(0.04, 5); fill(0.05, 5); 0.3; 0.35; 0.1]                    # grid spacing (heterogenous), meter (N=13)
# Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); 0.35; 0.1]          # grid spacing (heterogenous), meter (N=17)
# Δz_m = [fill(0.04, 5); fill(0.05, 5); fill(0.06, 5); fill(0.07, 5); 0.1] # grid spacing (heterogenous), meter (N=21)

f1 = (Δz_m) -> LWFBrook90.Rootden_beta_(0.97, Δz_m = Δz_m)  # function for root density as f(Δz)
f2 = (Δz_m) -> fill(-6.3, length(Δz_m))          # function for initial conditions as f(Δz)
input_soil_discretization = discretize_soil(;Δz_m = Δz_m, Rootden_ = f1, uAux_PSIM_init_kPa = f2)
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
ψM_initial, p = define_LWFB90_p(
    input_meteoveg,
    _input_meteoiso, # unused
    input_meteoveg_reference_date,
    input_param,
    input_storm_durations,
    input_soil_horizons,
    input_soil_discretization,
    simOption_FLAG_MualVanGen;
    Reset = Reset,
    # soil_output_depths = collect(-0.05:-0.05:-1.1),
    # soil_output_depths = [-0.1, -0.5, -1.0, -1.5, -1.9],
    compute_intermediate_quantities = compute_intermediate_quantities,
    simulate_isotopes = false)
####################

####################
# Define initial states of differential equation
# state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI
# Create u0 for DiffEq.jl
u0 = define_LWFB90_u0(p, input_initial_conditions,
    ψM_initial,
    compute_intermediate_quantities,
    simulate_isotopes = false)
####################

####################
# Define ODE problem which consists of
#   - definition of right-hand-side (RHS) function f
#   - definition of callback function cb
#   - u0:     initial condition of states
#   - tspan:  definition of simulation time span
#   - p:      parameters

# Define simulation time span:
# tspan = (0.0, 10.0) # simulate 5 days
# tspan = (0.,  100.) # simulate 100 days # NOTE: KAU bugs in "branch 005-" when at least 3*365
tspan = (minimum(input_meteoveg[:, "days"]),
         maximum(input_meteoveg[:, "days"])) # simulate all available days
# tspan = (LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1980,1,1), reference_date),
#          LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1985,1,1), reference_date)) # simulates selected period

# Define ODE:
ode_LWFBrook90, unstable_check_function = define_LWFB90_ODE(u0, tspan, p);
####################

####################
## Solve ODE:
@time sol_LWFBrook90 = solve(ode_LWFBrook90, Tsit5();
    progress = true,
    saveat = tspan[1]:tspan[2], dt = 1e-3, adaptive = true); # dt is initial dt, but adaptive
####################

####################
## Benchmarking
# @time sol_LWFBrook90 = solve(ode_LWFBrook90, progress = true;
#     saveat = tspan[1]:tspan[2], dt=1e-1, adaptive = false);
# @time sol_LWFBrook90 = solve(ode_LWFBrook90, progress = true, Euler(); # Note: Euler sometimes hangs
#     saveat = tspan[1]:tspan[2], dt=1e-1, adaptive = false);
using BenchmarkTools # for benchmarking
# sol_LWFBrook90 = @btime solve(ode_LWFBrook90; dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity
# sol_LWFBrook90 = @btime solve(ode_LWFBrook90; saveat = tspan[1]:tspan[2], dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity

## Performance optimizing
@time LWFBrook90.f_LWFBrook90R(copy(u0), u0, p, 1.0) # 0.000039 seconds (59 allocations: 5.047 KiB)
# @btime LWFBrook90.f_LWFBrook90R(copy(u0), u0, p, 1.0)
# @code_warntype LWFBrook90.f_LWFBrook90R(copy(u0), u0, p, 1.0)
@enter LWFBrook90.f_LWFBrook90R(copy(u0), u0, p, 1.0)
# integrator = init(ode_LWFBrook90, Tsit5();
#     progress = true,
#     saveat = tspan[1]:tspan[2], dt = 1e-3, adaptive = true)
# @time LWFBrook90.LWFBrook90R_update_INTS_INTR_SNOW_CC_SNOWLQ!(integrator) # 0.000048 seconds (119 allocations: 9.391 KiB)
# @enter LWFBrook90.LWFBrook90R_update_INTS_INTR_SNOW_CC_SNOWLQ!(integrator)
    # @time LWFBrook90.KPT.derive_auxiliary_SOILVAR(u_SWATI, p_soil) # 0.000010 seconds (10 allocations: 1.406 KiB)
# @time LWFBrook90.KPT.derive_auxiliary_SOILVAR(integrator.u[integrator.p[1][4][4]], integrator.p[1][1]) # 0.000013 seconds (11 allocations: 1.547 KiB)
####################


####################
## Plotting and exporting to CSV
mkpath("out")

ref_date = sol_LWFBrook90.prob.p[2][15]
time_to_plot = RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, ref_date) # define simulation times as dates

using Plots # Install plot package at first use with `]` and then `add Plots`

##### Plot 0 using plotting recipe defined in LWFBrook90.jl
using Plots, Measures
optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)

pl_final = LWFBrook90.plotlwfbrook90(sol_LWFBrook90, optim_ticks)
savefig(pl_final, joinpath("out", input_prefix *
    "_plotRecipe_NLAYER" * string(sol_LWFBrook90.prob.p[1][1].NLAYER) * ".png"))

##### Plot 1 very basic
# Plot a)
plot(sol_LWFBrook90; vars = [1, 2, 3, 4, 5, 6],label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

##### Plot 2
# http://docs.juliaplots.org/latest/generated/gr/#gr-ref43 # Heatmap with DateTime axis
(u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
    LWFBrook90.get_auxiliary_variables(sol_LWFBrook90)

n = sol_LWFBrook90.prob.p[1][1].NLAYER
# time_to_plot = #(see above)
y_cell_centers = [0; cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK)[1:(n-1)]] + sol_LWFBrook90.prob.p[1][1].p_THICK / 2
z = sol_LWFBrook90[7 .+ (0:sol_LWFBrook90.prob.p[1][1].NLAYER-1), 1, :]./sol_LWFBrook90.prob.p[1][1].p_THICK

heatmap(time_to_plot, y_cell_centers, u_aux_θ', yflip = true,
        xlabel = "Date",
        ylabel = "Depth [mm]",
        colorbar_title = "θ [-]")
heatmap(time_to_plot, y_cell_centers, u_aux_PSIM', yflip = true,
        xlabel = "Date",
        ylabel = "Depth [mm]",
        colorbar_title = "ψ [kPa]")

##### Plot 3
# aux_indices = sol_LWFBrook90.prob.p[1][4][5]
# aux_names = sol_LWFBrook90.prob.p[1][4][7]
# plot(sol_LWFBrook90; vars = aux_indices, label = aux_names)
# plot(sol_LWFBrook90; vars = aux_indices[17:25], label = aux_names[:, 17:25])

##### Plot 4: illustration how to extract certain depths and do manual calculations
depth_to_read_out_mm = [10 150 500 1000 1150]
pl_θ_fineSoil = plot(time_to_plot,
    LWFBrook90.get_θ(depth_to_read_out_mm, sol_LWFBrook90),
    labels = string.(depth_to_read_out_mm) .* "mm",
    title = "Plot using get_θ()",
    xlabel = "Date",
    ylabel = "θ [-]",
    legend = :bottomright)

# Some manual computations to include the gravel fraction
layer_indices = LWFBrook90.find_indices(depth_to_read_out_mm, sol_LWFBrook90)

# Variant 1: compute θ and correct with gravel fraction
θ_layers = u_aux_θ[:,layer_indices]
gravel_fraction_layers = sol_LWFBrook90.prob.p[1][1].p_STONEF[layer_indices]
θtotal_layers = θ_layers .* transpose(1 .- gravel_fraction_layers) # θtotal = θ*(1-gravelfrc) + 0.0*gravelfrc
pl_θ_totalVolume1 = plot(time_to_plot, θtotal_layers,
    labels = string.(depth_to_read_out_mm) .* "mm",
    title = "Plot using get_θ() and divide by p_STONEF",
    xlabel = "Date",
    ylabel = "θ total soil+rock volume [-]",
    legend = :bottomright)

# Variant 1: compute SWATI and correct with thickness of cell/layer
# u_SWATI # alternatively u_SWATI = transpose(sol_LWFBrook90[7 .+ (0:n-1), 1, :])
θtotal = u_SWATI ./ transpose(sol_LWFBrook90.prob.p[1][1].p_THICK)
θtotal_layers2 = θtotal[:, layer_indices]
pl_θ_totalVolume2 = plot(time_to_plot, θtotal_layers2,
    labels = string.(depth_to_read_out_mm) .* "mm",
    title = "Plot using u_SWATI and divide by p_THICK",
    xlabel = "Date",
    ylabel = "θ total soil+rock volume [-]",
    legend = :bottomright)

plot(pl_θ_fineSoil, pl_θ_totalVolume1, pl_θ_totalVolume2;
    layout=(3,1), size=(600,1200), left_margin = 10mm) # left_margin needed as soon as we use size
savefig(joinpath("out",input_prefix * "_different_θ_depths_NLAYER" * string(sol_LWFBrook90.prob.p[1][1].NLAYER) * ".png"))

##### Illustration how to export certain depths into a *.csv
using CSV, DataFrames
df_out = DataFrame(
    LWFBrook90.get_θ(depth_to_read_out_mm, sol_LWFBrook90),
    "θ_" .* string.(depth_to_read_out_mm[:]) .* "mm")
df_out.Date = time_to_plot
CSV.write(
    joinpath("out",input_prefix * "_θ_depths_NLAYER" * string(sol_LWFBrook90.prob.p[1][1].NLAYER) * ".csv"),
    df_out)
####################
