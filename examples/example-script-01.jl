# # Example Script 01
# This example was generated DATEOFTODAY

# Example data from Beatenberg is located in subfolder `examples/`. WSL is acknowledged for providing the input data (see section [Acknowledgments](@ref)).

# To run your own simulations: follow the step-by-step instructions below.

# ## Setup and run simulation

## How to work with a Julia script in a new environment
## create a new folder: "pgm"
## a) start `julia` -> `]` -> activate . -> add LWFBrook90
## b) in VSCode -> click at bottom and select "pgm" as Julia env
## c) run the code
## alternatively right-click on the folder in VSCode
##   -> Change to this Directory and Activate this environment

## Load packages:
using LWFBrook90

## `using LWFBrook90` brings the following exported functions into scope:
## - `loadSPAC()`                   # to define and run a simulation
## - `setup()`                      # to define and run a simulation
## - `simulate!()`                  # to define and run a simulation
## - `remakeSPAC()`                 # to define and run a simulation
## - `SPAC()`                       # to define and run a simulation
## - `DiscretizedSPAC()`            # to define and run a simulation
## - `get_states()`                 # to postprocess a simulation
## - `get_fluxes()`                 # to postprocess a simulation
## - `get_forcing()`                # to postprocess a simulation
## - `get_soil_()`                  # to postprocess a simulation
## - `get_water_partitioning()`     # to postprocess a simulation
## - `RelativeDaysFloat2DateTime()` # a helper function for time
## - `prepare_for_LWFBrook90R()`    # a helper function

# Documentation can be accessed by typing: `?loadSPAC`, `?setup`, `?simulate!`, `?plotisotopes`, etc.

# Define simulation model by reading in system definition and input data from input files.
# When printed, the generated SPAC model gives a summary.

## Read in input data
input_path = "../../../examples/DAV2020-full/"; input_prefix = "DAV2020-full"; # make sure to get the path correct relative to `pwd()``
## input_path = "examples/DAV2020-full/"; input_prefix = "DAV2020-full";
## input_path = "examples/isoBEAdense2010-18-reset-FALSE-input/"; input_prefix = "isoBEAdense2010-18-reset-FALSE";
model = loadSPAC(input_path, input_prefix; simulate_isotopes = false);
model = loadSPAC(input_path, input_prefix; simulate_isotopes = true);
model

# If wanted, arguments can be supplied to modify the loaded SPAC (e.g. for parameter estimation), hence requiring less input CSVs.
# Warnings are output (see below) to signal to the user that existing input CSVs might be ignored.
## Read in providing arguments, requiring less CSVs
model_modified =
    loadSPAC(input_path, input_prefix;
        simulate_isotopes = true,
        Δz_thickness_m    = [fill(0.04, 5); # grid spacing (heterogenous), meter (N=21)
                             fill(0.05, 5); # write Δ in VSCode by typing \Delta and hit shift
                             fill(0.06, 5);
                             fill(0.07, 5);
                             fill(0.10, 1)],
        root_distribution = (beta = 0.98, z_rootMax_m=-0.5),
        IC_soil           = (PSIM_init_kPa = -7.0,
                            delta18O_init_permil = -9.0,
                            delta2H_init_permil = -11.0),
        canopy_evolution  = (DENSEF_rel = 100,
                             HEIGHT_rel = 100,
                             SAI_rel    = 100,
                             LAI_rel = (DOY_Bstart = 120,
                                     Bduration  = 21,
                                     DOY_Cstart = 270,
                                     Cduration  = 60,
                                     LAI_perc_BtoC = 100,
                                     LAI_perc_CtoB = 55)),
        storm_durations_h = [5.44, 5.44, 5.44, 5.44, 5.44, 5.44,
                             5.44, 5.44, 5.44, 5.44, 5.44, 5.44],
        IC_scalar         = (amount = (u_GWAT_init_mm = 1.,
                                       u_INTS_init_mm = 13.7,
                                       u_INTR_init_mm = 0.,
                                       u_SNOW_init_mm = 22.222,
                                       u_CC_init_MJ_per_m2 = 0.101010,
                                       u_SNOWLQ_init_mm =  0.),
                            d18O    = (u_GWAT_init_permil = -11.111,
                                       u_INTS_init_permil = -12.222,
                                       u_INTR_init_permil = -13.333,
                                       u_SNOW_init_permil = -14.444),
                            d2H     = (u_GWAT_init_permil = -95.111,
                                       u_INTS_init_permil = -95.222,
                                       u_INTR_init_permil = -95.333,
                                       u_SNOW_init_permil = -95.444)));

# Documentation can be accessed by typing: `?loadSPAC`, `?setup`, `?simulate!`, `?plotisotopes`, etc.

# The objects `model` and `model_modified` are instances of a `SPAC()` (soil-plant-atmosphere continuum)
# and contain fully specified simulations. To run them they need to be transformed into a
# system of ODES (`setup()`) and solved (`simulate()`):

## Setup and run simulation
simulation     = setup(model)
simulation_mod = setup(model_modified)
## Inputs can be checked with:
## `plot_inputs(simulation)`. Note: NOT YET IMPLEMENTED

simulate!(simulation_mod)
simulate!(simulation_mod) # Run it a second time to showcase shorter runtime
## The computed solution is stored within the object.
## When using VSCode a progress bar is shown in the status bar.

# Notice that the very first time a model is simulated, the runtime is longer because it includes compilation.
# This compilation is not needed for any subsequent run.

# ## Postprocessing results
# These simulations can then be post-processed with predefined functions as shown below.
# Alternatively the following variables contain the simuation results for user-defined post-processing:
simulation_mod.ODESolution;
simulation_mod.ODESolution_datetime;
typeof(simulation_mod.ODESolution);
propertynames(simulation_mod.ODESolution)
## HINT: To get these names, type `simulation_mod.` and wait for the autocomplete!

# `simulation_mod.ODESolution` is a ODESolution object: documentation how to access under:
# https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/
# and more generally:
# https://docs.sciml.ai/Overview/stable/
# simulation_mod.ODESolution.t
# simulation_mod.ODESolution.u
# simulation_mod.ODESolution.retcode
# simulation_mod.ODESolution.destats

# Note that it is possible to use R Code from within Julia, e.g ggplot: https://stackoverflow.com/a/70073193/3915004


# ## Postprocessing: Extract values
# Below illustrative code how to extract states and fluxes:

## Get simulation output (states, fluxes, forcing) as DataFrames
# The three functions: `get_states()`, `get_fluxes()`, and `get_forcing()` return `DataFrames``
# with daily values of all the state variables, fluxes, and the input forcing.
get_states(simulation_mod)
get_fluxes(simulation_mod)
get_forcing(simulation_mod)

get_states(simulation_mod; days_to_read_out_d = 0:30:360) # specific days relative to simulation.parametrizedSPAC.reference_date can be specified


# Below an illustrative code snippet  how to export certain depths into a *.csv:
using CSV, DataFrames

## Single soil variables: How to get θ, or ψ, or δ18Os?
get_soil_(:θ, simulation_mod;
    depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)
depth_to_read_out_mm = [10, 150, 500, 1000, 1150]
get_soil_(:θ, simulation_mod;
    depths_to_read_out_mm = depth_to_read_out_mm, days_to_read_out_d = nothing)
get_soil_(:θ, simulation_mod;
    depths_to_read_out_mm = depth_to_read_out_mm)
get_soil_(:ψ, simulation_mod;
    depths_to_read_out_mm = depth_to_read_out_mm)
get_soil_(:δ18O, simulation_mod;
    depths_to_read_out_mm = depth_to_read_out_mm)

## How to export θ as CSV?
## Only every day, provide days_to_read_out
days_to_read_out = range(simulation_mod.ODESolution.prob.tspan...)
dates_to_read_out = LWFBrook90.RelativeDaysFloat2DateTime.(
    days_to_read_out, simulation_mod.parametrizedSPAC.reference_date)
df_out_daily = get_soil_(:θ, simulation_mod;
    depths_to_read_out_mm = depth_to_read_out_mm, days_to_read_out_d = days_to_read_out)

insertcols!(df_out_daily, 1, :dates => dates_to_read_out);
show(df_out_daily)
##plot(df_out_daily[:,:dates], Matrix(df_out_daily[:,Not([:dates, :time])]))
##CSV.write(
##    joinpath(out_dir, fname * "_θ_depths_daily.csv"),
##    df_out_daily)


# ## Postprocessing: Plotting (using provided functions)
# Below an example script using the provided plot recipes that plot a) amounts, b) isotopes, or c) forcing and states as an additional internal check:

using Plots, Measures; gr();
pl1 = plotamounts(simulation_mod, :above_and_belowground, :showRWUcentroid)
pl1
#-
pl2 = plotisotopes(simulation_mod, :d18O, (d18O = :auto, d2H = :auto), :showRWUcentroid)
pl2
#-
pl3 = plotforcingandstates(simulation_mod)
pl3

## Use of additional arguments to `plotisotopes` is illustrated here.
## Note that typing `?plotisotopes` provides the documentation to the function.
## plotisotopes(simulation_mod, :d2H)
## plotisotopes(simulation_mod, :d18O)
## plotisotopes(simulation_mod;
##              xlim = (DateTime("2010-01-01"), DateTime("2010-03-31")))
## plotisotopes(simulation_mod, :d2H, (d18O = :auto, d2H = :auto), :showRWUcentroid;
##              xlim = (DateTime("2010-01-01"), DateTime("2010-08-30")))

# Saving plots can be done with `savefig`

## using Dates: today;
## out_dir = joinpath("gitignored","out", string(today()))
## mkpath(out_dir)

## fname = input_prefix
## savefig(pl1, joinpath(out_dir, fname*"_plotRecipe_AMT.png"))
## savefig(pl2, joinpath(out_dir, fname*"_plotRecipe_ISO.png"))
## savefig(pl3, joinpath(out_dir, fname*"_plotRecipe_CHECK.png"))


# ## Postprocessing: Plotting (using your own functions)
# Below an example script using the manually written code to plot the simulation_mod.#

# sim_states = get_states(simulation) # this gives an error since it has not been simulated!
sim_states = get_states(simulation_mod)
sim_fluxes = get_fluxes(simulation_mod)

names(sim_states) # show column names
names(sim_fluxes) # show column names

# plot some of the scalar states
states_to_plot = ["INTS_mm", "INTR_mm", "SNOW_mm", "GWAT_mm", "SWAT_mm"]
sim_states_to_plot = sim_states[:, states_to_plot]
plot(sim_states[:,"dates"], Matrix(sim_states_to_plot),
    labels=permutedims(names(sim_states_to_plot)), legend=:topleft, ylabel = "-")


# plot the fates of water partitioning
fluxes_to_plot = ["flow", "seep", "evap"]
sim_fluxes_to_plot = sim_fluxes[:, fluxes_to_plot]
plot(sim_fluxes[:,"dates"], Matrix(sim_fluxes_to_plot),
    labels=permutedims(names(sim_fluxes_to_plot)), legend=:topleft, ylabel = "mm/day")

# plot internal fluxes
fluxes_to_plot = [
    "cum_d_prec", "cum_d_rnet", "cum_d_smlt",
    "cum_d_tran",
    "srfl", "slfl", "byfl", "dsfl", "gwfl", "vrfln"]
sim_fluxes_to_plot = sim_fluxes[:, fluxes_to_plot]
plot(sim_fluxes[:,"dates"], Matrix(sim_fluxes_to_plot),
    labels=permutedims(names(sim_fluxes_to_plot)), legend=:topleft, ylabel = "mm/day")

# plot isotope signatures of fluxes
fluxes_to_plot = ["RWU_d18O", "RWU_d2H", "PREC_d18O", "PREC_d2H"]
sim_fluxes_to_plot = sim_fluxes[:, fluxes_to_plot]
plot(sim_fluxes[:,"dates"], Matrix(sim_fluxes_to_plot),
    labels=permutedims(names(sim_fluxes_to_plot)), legend=:topleft, ylabel = "mm/day")


# Belowground quantities (θ,ψ,δ of soil water)

# if you want some specific depths: use get_soil_(): e.g
# depth_to_read_out_mm = [150, 500, 800, 1500]
# df_δsoil = get_soil_([:δ18O, :δ2H], simulation_mod; depths_to_read_out_mm = depth_to_read_out_mm)
# df_θψ = get_soil_([:θ, :ψ],
#     simulation_mod, depths_to_read_out_mm = depth_to_read_out_mm)

# if we use depths that are already contained in sim_states we can use those:
@show propertynames(sim_states[:,r"(d18O)|(d2H)"]);
df_δsoil = sim_states[:, ["d18O_permil_160mm", "d18O_permil_510mm",
                          "d18O_permil_820mm", "d18O_permil_1200mm",
                          "d2H_permil_160mm", "d2H_permil_510mm",
                          "d2H_permil_820mm", "d2H_permil_1200mm"]]
df_θψ    = sim_states[:, ["θ_m3m3_160mm", "θ_m3m3_510mm",
                          "θ_m3m3_820mm", "θ_m3m3_1200mm",
                          "ψ_kPa_160mm", "ψ_kPa_510mm",
                          "ψ_kPa_820mm", "ψ_kPa_1200mm"]]
pl_θ = plot(sim_states.dates,
    Matrix(select(df_θψ, r"θ_")),
    labels = permutedims(names(select(df_θψ, r"θ_"))),
    xlabel = "Date",
    ylabel = "θ\n[-]",
    legend = :outerright)
pl_ψ = plot(sim_states.dates,
    ## Matrix(select(df_θψ, r"ψ_")),
    log10.(-Matrix(select(df_θψ, r"ψ_"))), yflip = true,
    labels = permutedims(names(select(df_θψ, r"ψ_"))),
    xlabel = "Date",
    ylabel = "log10(ψ\n[kPa])",
    legend = :outerright);
pl_δ18O = plot(sim_states.dates,
    Matrix(select(df_δsoil, r"d18O_")),
    labels = permutedims(names(select(df_δsoil, r"d18O_"))),
    xlabel = "Date",
    ylabel = "δ¹⁸O soil\n[‰]",
    legend = :outerright);
pl_δ2H = plot(sim_states.dates,
    Matrix(select(df_δsoil, r"d2H_")),
    labels = permutedims(names(select(df_δsoil, r"d2H_"))),
    xlabel = "Date",
    ylabel = "δ²H soil\n[‰]",
    legend = :outerright);
## add precipitation to soil δ
PREC_color = :black
plot!(pl_δ18O, sim_fluxes.dates, sim_fluxes.PREC_d18O,
    labels = "PREC", color = PREC_color, linestyle = :dot);
plot!(pl_δ2H, sim_fluxes.dates, sim_fluxes.PREC_d2H,
    labels = "PREC", color = PREC_color, linestyle = :dot);


pl_PREC = plot(
    sim_fluxes.dates,
    sim_fluxes.cum_d_prec,
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
    size=(600,500), dpi = 300, margin = 5mm)
# ####################