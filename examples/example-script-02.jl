# # Example Script 02
# This example was generated DATEOFTODAY

# ## Soil evaporation depth and isotope fractionation
# Compare a 10 cm evaporation source zone without fractionation, the same zone with
# Craig–Gordon fractionation, and a 30 cm zone with fractionation. Comparing cases 1
# and 2 isolates fractionation; comparing cases 2 and 3 isolates the source depth.
# See [Soil evaporation depth and isotope fractionation](@ref soil-evaporation-settings) for the settings and
# assumptions. `SLVPDEPTH_m` defaults to `0.0` (the uppermost input layer), and soil
# evaporation fractionation defaults to `false`.

# ## Load the common inputs
# Use the bundled Davos data and a common grid with 3 cm and 7 cm surface layers,
# followed by 10 cm layers. This resolves both the 10 cm and 30 cm source zones.
# The input path is independent of the working directory.

using LWFBrook90
using Plots

input_path = normpath(joinpath(dirname(pathof(LWFBrook90)), "..", "examples", "DAV2020-full"))
model = loadSPAC(input_path, "DAV2020-full";
    simulate_isotopes = true,
    simulate_evaporation_fractionation = false,
    Δz_thickness_m = [0.03, 0.07, fill(0.10, 10)...],
    root_distribution = (beta = 0.97, z_rootMax_m = -0.5),
    IC_soil = (PSIM_init_kPa = -6.3,
        delta18O_init_permil = -13.0, delta2H_init_permil = -95.0))

# ## Run three cases
# Use a short summer period, in days relative to the model's reference date.
# Each case starts from the same prescribed initial state at day 150; this is a
# sensitivity exercise rather than a simulation spun up from the start of the year.
# Rainfall, transpiration, and soil-water transport remain active in all cases.
# `remakeSPAC` returns a simulation ready to run, and `(SLVPDEPTH_m = depth,)`
# supplies the parameter in meters. Fractionation is a solver option and requires
# isotope simulation to be enabled.

simulation_period = (150.0, 180.0)
saved_days = range(simulation_period...; step = 1.0)
time_edges_days = [first(saved_days) - 0.5; saved_days .+ 0.5]
cases = (
    (label = "10 cm, no fractionation", depth_m = 0.10, fractionates = false),
    (label = "10 cm, fractionation", depth_m = 0.10, fractionates = true),
    (label = "30 cm, fractionation", depth_m = 0.30, fractionates = true),
)
simulations = map(cases) do case
    simulation = remakeSPAC(model;
        params = (SLVPDEPTH_m = case.depth_m,),
        solver_options = (simulate_evaporation_fractionation = case.fractionates,),
        requested_tspan = simulation_period)
    simulate!(simulation; progress = false, save_everystep = false, saveat = saved_days)
    simulation
end;

# ## Extract soil isotope signatures
# `get_soil_` returns a DataFrame with isotope signatures in per mil (‰).
# Its requested depths are in millimeters, so 15, 65, 150, and 250 mm sample the
# centers of the four layers in the top 30 cm of this grid.

soil_isotopes = get_soil_([:d18O, :d2H], simulations[2];
    depths_to_read_out_mm = [15, 65, 150, 250])
first(soil_isotopes, 5)

# ## Compare the soil profiles through time
# Plot the top 40 cm to show the source zones and the soil immediately below them.
# Rows show δ¹⁸O and δ²H; columns show the three cases. All columns in each row
# share the same color limits. Explicit depth edges preserve the unequal layer
# thicknesses, and the dashed line marks the evaporation source depth.

isotope_heatmaps = Plots.Plot[]
for (isotope, isotope_label) in ((:d18O, "δ¹⁸O (‰)"), (:d2H, "δ²H (‰)"))
    profiles = map(simulations) do simulation
        [getproperty(u.SWATI, isotope)[layer]
         for layer in eachindex(simulation.ODEProblem.p.p_soil.p_THICK),
             u in simulation.ODESolution.u]
    end
    shared_clims = extrema(Iterators.flatten(profiles))
    for (case, simulation, profile) in zip(cases, simulations, profiles)
        depth_edges_cm = [0.0; cumsum(simulation.ODEProblem.p.p_soil.p_THICK) ./ 10]
        panel = heatmap(time_edges_days, depth_edges_cm, profile;
            title = case.label, titlefontsize = 10,
            xlabel = "Time (days from reference date)", ylabel = "Soil depth (cm)",
            yflip = true, ylims = (0.0, 40.0), yticks = [0, 3, 10, 20, 30, 40],
            xlims = simulation_period, clims = shared_clims,
            color = :viridis, colorbar_title = isotope_label)
        hline!(panel, [100 * case.depth_m]; color = :white, linestyle = :dash, label = false)
        push!(isotope_heatmaps, panel)
    end
end
comparison_plot = plot(isotope_heatmaps...; layout = (2, 3), size = (1200, 700), link = :both)
display(comparison_plot) #src

# During drying, the shallow fractionating case can enrich more rapidly near the
# surface, while the deeper case distributes enrichment over the top 30 cm.
# Here rainfall and mixing may interrupt that pattern. Compare the fractionating
# 10 cm case with its non-fractionating control before attributing a change to
# evaporation. Atmospheric vapor isotopes are inferred from equilibrium with
# precipitation; they are not independent measured forcing in this example.

## To save the comparison, uncomment this line:
## savefig(comparison_plot, "soil-evaporation-isotopes.png")
