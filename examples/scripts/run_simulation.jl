# Run this on the cmd line like:
#  julia --project=. run_simulation.jl [] []
#  julia --project=. run_simulation.jl "../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl"
# "Hammel_loam-NLayer-27-RESET=FALSE"
#
# Note: I am aware that this isn't efficient as LWFBrook90 gets recompiled each time. See:
#       https://stackoverflow.com/questions/50608970/if-a-julia-script-is-run-from-the-command-line-does-it-need-to-be-re-compiled-e
#       See: specifically https://stackoverflow.com/a/42040763 from
#       https://stackoverflow.com/questions/42030150/how-to-exectute-julia-code-in-command-line
#       See:
#       https://discourse.julialang.org/t/how-to-pass-command-line-arguments-to-include-julia-file/4731

using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5, init
using Dates: now
using CSV: read, File
using DataFrames: DataFrame, rename, unstack, Not# ,select

# Parse command line arguments and call run_simulation() (if called from command line)
sim_result = LWFBrook90.run_simulation(ARGS)
plot_and_save_results(sim_result...)

# Run this script interactively many simulations with run_simulation()
# The ones with Reset=FALSE
# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_sand-NLayer-27-RESET=FALSE"])
# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_sand-NLayer-103-RESET=FALSE"])
# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_sand-NLayer-400-RESET=FALSE"])

# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_loam-NLayer-27-RESET=FALSE"])
# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_loam-NLayer-103-RESET=FALSE"])
# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_loam-NLayer-400-RESET=FALSE"])

# The ones with Reset=TRUE (even if this is not implemented in Julia, but just for the comparison)
# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_sand-NLayer-27-RESET=TRUE"])
# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_sand-NLayer-103-RESET=TRUE"])
# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_sand-NLayer-400-RESET=TRUE"])

# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_loam-NLayer-27-RESET=TRUE"])
# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_loam-NLayer-103-RESET=TRUE"])
# run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_loam-NLayer-400-RESET=TRUE"])

# sol1 = run_simulation(["examples/isoBEAdense2010-18-reset-FALSE-input/" "isoBEAdense2010-18-reset-FALSE" "true"])
# # sol1[1] # solution
# # sol1[2] # input_prefix
# # sol1[3] # input_path
# using Plots, Measures
# optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
# # pl1 = LWFBrook90.ISO.plotisotopes(sol_LWFBrook90);
# pl2 = LWFBrook90.ISO.plotisotopes(
#     sol1[1], optim_ticks;
#     layout = grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4]),
#     size=(1000,1400), dpi = 300, leftmargin = 15mm);
# plot!(pl2, link = :x)

# Test the plotting and saving
# sol1 = run_simulation(["../Unit_Tests/Hammel-IntegrationTests-LWFBrook90/input_LWFBrook90.jl/" "Hammel_sand-NLayer-27-RESET=FALSE"])
# plot_and_save_results(sol1...)
