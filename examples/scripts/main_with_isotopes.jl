# # Example isotope simulation
# This script illustrates a script that runs a isotope simulation.

# ## Loading of the package LWFBrook90 and other dependencies:

using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5, init

# ## Define input data and reading it in:
input_path = "examples/isoBEAdense2010-18-reset-FALSE-input/"; input_prefix = "isoBEAdense2010-18-reset-FALSE"
model = loadSPAC(input_path, input_prefix; simulate_isotopes = true);

# ## Define grid for spatial discretization as well as initial conditions and root densities
root_distribution = (beta = 0.97, )
Δz_m = [0.04, 0.04, 0.12, 0.25, 0.3, 0.35, 0.1]; # grid spacing (heterogenous), meter (N=7)

# f1 = (Δz_m) -> LWFBrook90.Rootden_beta_(root_distribution.beta, Δz_m = Δz_m)
# f2 = (Δz_m) -> fill(-6.3, length(Δz_m))
# input_soil_discretization = discretize_soil(;
#     Δz_m = Δz_m,
#     Rootden_ = f1,
#     uAux_PSIM_init_kPa = f2,
#     u_delta18O_init_permil = ifelse.(cumsum(Δz_m) .<= 0.2, -13., -10.),
#     u_delta2H_init_permil  = ifelse.(cumsum(Δz_m) .<= 0.2, -95., -70.));

# If wanted, arguments can be supplied to modify the loaded SPAC (e.g. for parameter estimation):
model_modified =
    loadSPAC(input_path, input_prefix;
        simulate_isotopes = true,
        Δz_thickness_m    = Δz_m,
        root_distribution = root_distribution,
        IC_soil           = (PSIM_init_kPa = -7.0, delta18O_init_permil = -9.0, delta2H_init_permil = -11.0));

simulation_modified = setup(model_modified)
simulate!(simulation_modified)

# ## Plotting
using DataFrames
using CSV: File
using Plots, Measures
## fname = joinpath(
##     "out",
##     input_prefix*"_NLAYER+" * string(simulation_modified.ODESolution.prob.p[1][1].NLAYER)*
##     "-git+"*chomp(Base.read(`git rev-parse --short HEAD`, String)))*
##     ifelse(length(read(`git status --porcelain`, String))==0, "+gitclean","+gitdirty")*
##     "-iso+"*string(simulate_isotopes)

# Resulting plot of some aboveground quantities:
# ref_aboveground =
#     DataFrame(File(
#     "../../test/test-assets-external/BEA-2016/BEA-IntegrationTests-LWFBrook90/output_LWFBrook90R/BEA2016-reset-FALSE_NLAYER14_LWFBrook90R-0.4.5daily_output.csv"))[:,
#     # "../../gitignored/out/BEA2016-reset-FALSE_NLAYER14_LWFBrook90R-0.4.5daily_output.csv"))[:,
#         [:yr, :mo, :da, :doy, :intr, :ints, :snow, :gwat]]
pl_ab_3 = plot(simulation_modified.ODESolution; vars = [2, 3, 4],
    label=["INTS (mm)" "INTR (mm)" "SNOW (mm)"])
plot!(pl_ab_3,
        [ref_aboveground.intr,
        ref_aboveground.ints,
        ref_aboveground.snow], label = "LWFBrook90R", line = :dash, color = :black)
## savefig(fname*"_plot-INTS_INTR_SNOW.png")

# Resulting plot of other aboveground quantities:
pl_ab_4 = plot(simulation_modified.ODESolution; vars = [2, 3],
    label=["INTS (mm)" "INTR (mm)"])
plot!(pl_ab_4,
        [ref_aboveground.intr,
        ref_aboveground.ints], label = "LWFBrook90R", line = :dash, color = :black)
## savefig(fname*"_plot-INTS_INTR.png")

pl2 = plotisotopes(simulation_modified, :d18O, (d18O = :auto, d2H = :auto), :showRWUcentroid)
pl2 = plotisotopes(simulation_modified, :d18O_and_d2H, (d18O = :auto, d2H = :auto), :showRWUcentroid)


# Belowground quantities (θ,ψ,δ of soil water)
PREC_color = :black
depth_to_read_out_mm = [150 500 800 1500]
if simulate_isotopes
    δ_resultsSoil = get_δsoil(simulation_modified, depths_to_read_out_mm = depth_to_read_out_mm)
    δ_results = get_δ(simulation_modified)
end

pl_θ = plot(simulation_modified.ODESolution_datetime,
    # TODO: replace get_θ(...) by get_soil_(:θ, ...)
    get_θ(simulation_modified, depths_to_read_out_mm = depth_to_read_out_mm)',  # TODO
    labels = string.(depth_to_read_out_mm) .* "mm",
    xlabel = "Date",
    ylabel = "θ\n[-]",
    legend = :outerright);
pl_ψ = plot(simulation_modified.ODESolution_datetime,
    ## -LWFBrook90.get_ψ(depth_to_read_out_mm, simulation_modified.ODESolution) .+ 1, yaxis = :log, yflip = true,
    get_ψ(simulation_modified, depths_to_read_out_mm = depth_to_read_out_mm)',
    # TODO: replace get_θ(...) by get_soil_(:ψ, ...)
    labels = string.(depth_to_read_out_mm) .* "mm",
    xlabel = "Date",
    ylabel = "ψ\n[kPa]",
    legend = :outerright);

if simulate_isotopes
    pl_δ18O = plot(simulation_modified.ODESolution_datetime,
        δ_resultsSoil.d18O',
        labels = string.(depth_to_read_out_mm) .* "mm",
        xlabel = "Date",
        ylabel = "δ¹⁸O soil\n[‰]",
        legend = :outerright);
    pl_δ2H = plot(simulation_modified.ODESolution_datetime,
        δ_resultsSoil.d2H',
        labels = string.(depth_to_read_out_mm) .* "mm",
        xlabel = "Date",
        ylabel = "δ²H soil\n[‰]",
        legend = :outerright);
    ## add precipitation to soil δ
    plot!(pl_δ2H,
        simulation_modified.ODESolution_datetime,
        δ_results.PREC_d2H, labels = "PREC", color = PREC_color, linestyle = :dot);
    plot!(pl_δ18O,
        simulation_modified.ODESolution_datetime,
        δ_results.PREC_d18O, labels = "PREC", color = PREC_color, linestyle = :dot);
else
    pl_δ18O = plot();
    pl_δ2H = plot();
end

pl_PREC = plot(
    simulation_modified.ODESolution_datetime,
    simulation_modified.ODESolution.prob.p.p_PREC.(simulation_modified.ODESolution.t),
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