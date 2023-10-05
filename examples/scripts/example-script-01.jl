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

#jl ## `using LWFBrook90` brings the following exported functions into scope:
#jl ## - `SPAC()`
#jl ## - `DiscretizedSPAC()`
#jl ## - `discretize()`
#jl ## - `simulate!()`
#jl ## - `get_soil_()`
#jl ## - `get_aboveground()`
#jl ## - `get_δ()`,

# Define simulation model by reading in system definition and input data from input files

## Read in input data
println(pwd())
@info pwd() # to debug Documenter.jl #hide
@warn pwd() # to debug Documenter.jl #hide
println(@__DIR__)
@info @__DIR__ # to debug Documenter.jl #hide
@warn @__DIR__ # to debug Documenter.jl #hide
println(readdir(pwd()))
@info readdir(pwd()) # to debug Documenter.jl #hide
@warn readdir(pwd()) # to debug Documenter.jl #hide
println(readdir("../"))
@info readdir("../") # to debug Documenter.jl #hide
@warn readdir("../") # to debug Documenter.jl #hide

# #src input_path = "examplesisoBEAdense2010-18-reset-FALSE-input/"; input_prefix = "isoBEAdense2010-18-reset-FALSE";
# input_path = "examples/DAV2020-full/"; input_prefix = "DAV2020-full";
# ## input_path = "examples/isoBEAdense2010-18-reset-FALSE-input/"; input_prefix = "isoBEAdense2010-18-reset-FALSE";
# model = loadSPAC(input_path, input_prefix; simulate_isotopes = false);
# model = loadSPAC(input_path, input_prefix; simulate_isotopes = true);
# model

# # If wanted, arguments can be supplied to modify the loaded SPAC (e.g. for parameter estimation), hence requiring less input CSVs.
# ## Read in providing arguments, requiring less CSVs
# model_modified =
#     loadSPAC(input_path, input_prefix;
#         simulate_isotopes = true,
#         Δz_thickness_m    = [fill(0.04, 5); # grid spacing (heterogenous), meter (N=21)
#                              fill(0.05, 5); # write Δ in VSCode by typing \Delta and hit shift
#                              fill(0.06, 5);
#                              fill(0.07, 5);
#                              fill(0.10, 1)],
#         root_distribution = (beta = 0.98, z_rootMax_m=-0.5),
#         IC_soil           = (PSIM_init_kPa = -7.0,
#                             delta18O_init_permil = -9.0,
#                             delta2H_init_permil = -11.0),
#         canopy_evolution  = (DENSEF_rel = 100,
#                              HEIGHT_rel = 100,
#                              SAI_rel    = 100,
#                              LAI_rel = (DOY_Bstart = 120,
#                                      Bduration  = 21,
#                                      DOY_Cstart = 270,
#                                      Cduration  = 60,
#                                      LAI_perc_BtoC = 100,
#                                      LAI_perc_CtoB = 55)),
#         storm_durations_h = [5.44, 5.44, 5.44, 5.44, 5.44, 5.44,
#                              5.44, 5.44, 5.44, 5.44, 5.44, 5.44],
#         IC_scalar         = (amount = (u_GWAT_init_mm = 1.,
#                                        u_INTS_init_mm = 13.7,
#                                        u_INTR_init_mm = 0.,
#                                        u_SNOW_init_mm = 22.222,
#                                        u_CC_init_MJ_per_m2 = 0.101010,
#                                        u_SNOWLQ_init_mm =  0.),
#                             d18O    = (u_GWAT_init_permil = -11.111,
#                                        u_INTS_init_permil = -12.222,
#                                        u_INTR_init_permil = -13.333,
#                                        u_SNOW_init_permil = -14.444),
#                             d2H     = (u_GWAT_init_permil = -95.111,
#                                        u_INTS_init_permil = -95.222,
#                                        u_INTR_init_permil = -95.333,
#                                        u_SNOW_init_permil = -95.444)));

# # Documentation can be accessed by typing: `?loadSPAC`, `?setup`, `?simulate!`, `?plotisotopes`, etc.

# # The objects `model` and `model_modified` are instances of a `SPAC()` (soil-plant-atmosphere continuum)
# # and contain fully specified simulations. To run them they need to be transformed into a
# # system of ODES (`setup()`) and solved (`simulate()`):

# ## Setup and run simulation
# simulation          = setup(model)
# simulation_modified = setup(model_modified)
# ## Inputs can be checked with:
# ## `plot_inputs(simulation)`. Note: NOT YET IMPLEMENTED

# simulate!(simulation)
# simulate!(simulation_modified)
# ## The computed solution is stored within the object.
# ## When using VSCode a progress bar is shown in the status bar.

# # ## Postprocessing results
# # These simulations can then be post-processed with predefined functions.
# # Alternatively the followig variables contain the simuation results for user-defined post-processing:
# simulation.ODESolution;
# simulation.ODESolution_datetime;
# typeof(simulation.ODESolution);
# ## HINT: To get these names, type `simulation.` and wait for the autocomplete!

# # `simulation.ODESolution` is a ODESolution object: documentation how to access under:
# # `https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/`
# # and more generally:
# # `https://docs.sciml.ai/Overview/stable/`
# # simulation.ODESolution.t
# # simulation.ODESolution.u
# # simulation.ODESolution.retcode
# # simulation.ODESolution.destats

# # Note that it is possible to use R Code from within Julia, e.g ggplot:
# # `https://stackoverflow.com/a/70073193/3915004`

# # ### Plotting (using provided functions)
# # Below an example script using the provided plot recipes that plot a) amounts, b) isotopes, or c) forcing and states as an additional internal check:

# using Plots, Measures; gr();

# pl1 = plotamounts(simulation_modified, :above_and_belowground, :showRWUcentroid)
# pl1

# pl2 = plotisotopes(simulation_modified, :d18O, (d18O = :auto, d2H = :auto), :showRWUcentroid)
# pl2

# pl3 = plotforcingandstates(simulation_modified)
# pl3

# ## Use of additional arguments to `plotisotopes` is illustrated here.
# ## Note that typing `?plotisotopes` provides the documentation to the function.
# ## plotisotopes(simulation_modified, :d2H)
# ## plotisotopes(simulation_modified, :d18O)
# ## plotisotopes(simulation_modified;
# ##              xlim = (DateTime("2010-01-01"), DateTime("2010-03-31")))
# ## plotisotopes(simulation_modified, :d2H, (d18O = :auto, d2H = :auto), :showRWUcentroid;
# ##              xlim = (DateTime("2010-01-01"), DateTime("2010-08-30")))

# # Saving plots can be done with `savefig`

# using Dates: today;
# out_dir = joinpath("gitignored","out", string(today()))
# mkpath(out_dir)

# fname = input_prefix
# savefig(pl1, joinpath(out_dir, fname*"_plotRecipe_AMT.png"))
# savefig(pl2, joinpath(out_dir, fname*"_plotRecipe_ISO.png"))
# savefig(pl3, joinpath(out_dir, fname*"_plotRecipe_CHECK.png"))


# # ### Extract values
# # Below an illustrative code snippet  how to export certain depths into a *.csv:

# using CSV, DataFrames

# ## How to get θ?
# get_soil_(:θ, simulation_modified;
#     depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)
# depth_to_read_out_mm = [10, 150, 500, 1000, 1150]
# get_soil_(:θ, simulation_modified;
#     depths_to_read_out_mm = depth_to_read_out_mm, days_to_read_out_d = nothing)
# get_soil_(:θ, simulation_modified;
#     depths_to_read_out_mm = depth_to_read_out_mm)

# ## How to export θ as CSV?
# ## Only every day, provide days_to_read_out
# days_to_read_out = range(simulation_modified.ODESolution.prob.tspan...)
# dates_to_read_out = LWFBrook90.RelativeDaysFloat2DateTime.(
#     days_to_read_out, simulation_modified.parametrizedSPAC.reference_date)
# df_out_daily = get_soil_(:θ, simulation_modified;
#     depths_to_read_out_mm = depth_to_read_out_mm, days_to_read_out_d = days_to_read_out)

# insertcols!(df_out_daily, 1, :dates => dates_to_read_out)

# plot(df_out_daily[:,:dates], Matrix(df_out_daily[:,Not([:dates, :time])]))
# CSV.write(
#     joinpath(out_dir, fname * "_θ_depths_daily.csv"),
#     df_out_daily)
# ####################