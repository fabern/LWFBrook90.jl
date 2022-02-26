# TODO: maybe move this into a module? on soil discretization?

using DataFrames: DataFrame, rename# ,select
using DataFramesMeta#: @linq, transform, DataFramesMeta

"""
    read_inputData(folder::String, prefix::String)

Load different input files for LWFBrook90:
- meteoveg.csv
- param.csv
- meteo_storm_durations.csv
- initial_conditions.csv
- soil_horizons.csv
- soil_discretization.csv

These files were created with an R script `generate_LWFBrook90Julia_Input.R` that
takes the same arguements as the R funciton `LWFBrook90R::run_LWFB90()` and generates
the corresponding input files for LWFBrook90.jl.
"""
function discretize_soil(folder::String, prefix::String; suffix::String = "")
    # Load pre-defined discretization
    path_soil_discretization = joinpath(folder, prefix * "_soil_discretization" * suffix * ".csv")
    input_soil_discretization = LWFBrook90.read_path_soil_discretization(path_soil_discretization)

    return input_soil_discretization
end



"""
    Rootden_beta_(
    β;
    h_m,
    total_effective_rooting_depth_m = maximum(cumsum(h_m)),
    z_Upper_m = 0)

Define the relative root density in each discretized soil layer based on the beta model from
(Gale and Grigal, 1987).

This function returns the instantaneous root fraction (dY/dd) (units of -/cm). Unless cropped
    by the total effective rooting depth the area under the curve sums up to 1.0 (when plotted vs cm).
    The method reduces the amount of roots in a discretization layer in case the effective maximal rooting
    depth comes to lie within that layer - it does not need to modify the discretization in that case.
"""
function Rootden_beta_(
    β;
    h_m,
    total_effective_rooting_depth_m = maximum(cumsum(h_m)),
    z_Upper_m = 0) # # Upper interface of topmost discretization cell
    # e.g. h_m = fill(0.10, 10)
    #      total_effective_rooting_depth_m = 0.22
    #      β = 0.97


    # Notation:
    # z: absolute vertical position (positive upward, i.e. soil values are negative)
    # d: depth below beginning of mineral soil
    @assert z_Upper_m <= 0 "Discretization of humus layer (>0m) is untested."
    @assert minimum(h_m) > 0 "Discretization cell height h_m must be positive"
    z_interfaces_m = z_Upper_m .+ [0; cumsum(-h_m)]

    # (Gale and Grigal, 1987)
    Y_fct(β, d_cm) = 1 .- β .^ d_cm
    dYdd_fct(β, d_cm) = -β .^ d_cm .* log(β) # unused

    # # (Gale and Grigal, 1987) figures 1 and 2
    # d_midpoints_m = -z_Upper_m + h_m[1]/2 .+
    #     cumsum(h_m) .- h_m[1] # d is soil depth d in centimeters
    # # TODO: by using d_midpoint of we make a small error.
    # #       Better would be to use integral between cell boundaries.
    # d_cm = d_midpoints_m*100
    # pl_1 = plot(
    #     Y_fct(0.92, d_cm), d_cm,
    #     xflip = true,
    #     yflip = true, ylabel = "Depth_cm", xlabel = "Cumulative Root Fraction (Y)");
    # plot!(Y_fct(0.95, d_cm), d_cm);
    # plot!(Y_fct(0.97, d_cm), d_cm);
    # pl_2 = plot(
    #     dYdd_fct(0.92, d_cm), d_cm,
    #     yflip = true, ylabel = "Depth_cm", xlabel = "Instantaneous Root Fraction (dY/dd)");
    # plot!(dYdd_fct(0.95, d_cm), d_cm);
    # plot!(dYdd_fct(0.97, d_cm), d_cm);
    # plot(pl_1, pl_2; layout=(1,2))

    # 0) Preparatory steps
    # Depth of cell interfaces at which we evaluate the cumulative root fraction
    d_interfaces_cm = -z_interfaces_m * 100
    # # Crop cumulative density at total effective rooting depth
    d_root_depth_cm = total_effective_rooting_depth_m * 100

    # 1) Compute cumulative density inside of a cell (i.e. increase in cumulative denisity from upper to lower interface)
    rootden_increase_per_cell = diff(Y_fct(β, d_interfaces_cm))

    # 2) Correct for total effective rooting depth
    #   a) total rooting depth is below lowest discretization cell                    (nothing needs changing)
    #   b) total rooting depth coincides with lower boundary of a discretization cell (remove roots further below)
    #   c) total rooting depth is inside of a discretization cell                     (remove roots below and correct cell)
    d_upper_interfaces_cm = d_interfaces_cm[Not(end),]
    d_lower_interfaces_cm = d_interfaces_cm[Not(1),]

    # remove roots from cells that are completely below rooting depth
    rootden_increase_per_cell[d_upper_interfaces_cm.>=d_root_depth_cm] .= 0
    rootden_increase_per_cell

    # correct root density in cell where total rooting depth lies
    idx = findfirst(d_root_depth_cm .<= d_lower_interfaces_cm)
    if (!isnothing(idx))
        # First computation above was:
        # Y_fct(β, d_lower_interfaces_cm[idx]) - Y_fct(β, d_upper_interfaces_cm[idx])
        # Corrected computation is:
        rootden_increase_per_cell[idx] =
            Y_fct(β, d_root_depth_cm) - Y_fct(β, d_upper_interfaces_cm[idx])
        # (note if d_lower_interfaces_cm[idx] == d_root_depth_cm, this keeps the value before)
    end

    # 3) Normalize the increase in root density over the cell by the cell height
    #    note that this uses h_m also for the one cell which contains the total rooting depth
    rootden_increase_per_cmDepth = rootden_increase_per_cell ./ (h_m .* 100)
end

# Rootden_beta_(0.97, h_m = fill(0.10, 10))
# Rootden_beta_(0.97, h_m = fill(0.10, 10), total_effective_rooting_depth_m = 0.2)

# Plot Rootden_beta_
# h_m = fill(0.10, 10) # (m)
# h_m = diff(collect(0:0.01:1)) # (m)
# total_effective_rooting_depth_m = 0.225
# plot(Rootden_beta_(0.97, h_m = h_m),
#     cumsum(h_m * 100),
#     ylabel = "Depth (cm)", xlabel = "Instantaneous root fraction (dY/dd) (Area of 1.0, unless cropped)",
#     yflip = true, seriestype = :path, legend = :bottomright,
#     label = "No max Root")
# # plot!(Rootden_beta_(0.97, h_m = h_m, total_effective_rooting_depth_m = total_effective_rooting_depth_m),
# #     cumsum(h_m * 100), seriestype = [:scatter], label = "maxRoot: $(total_effective_rooting_depth_m)m")
# # plot!(Rootden_beta_(0.97, h_m = h_m, total_effective_rooting_depth_m = 1.2),
# #     cumsum(h_m * 100), seriestype = [:scatter], label = "maxRoot: $(1.2)m")
# # plot!(Rootden_beta_(0.97, h_m = h_m, total_effective_rooting_depth_m = 0.5),
# #     cumsum(h_m * 100), seriestype = [:scatter], label = "maxRoot: $(0.5)m",)
# plot!(Rootden_beta_(0.97, h_m = h_m, total_effective_rooting_depth_m = total_effective_rooting_depth_m),
#     cumsum(h_m * 100), seriestype = [:path], label = "maxRoot: $(total_effective_rooting_depth_m)m")
# plot!(Rootden_beta_(0.97, h_m = h_m, total_effective_rooting_depth_m = 0.9),
#     cumsum(h_m * 100), seriestype = [:path], label = "maxRoot: $(0.9)m")
# plot!(Rootden_beta_(0.97, h_m = h_m, total_effective_rooting_depth_m = 0.5),
#     cumsum(h_m * 100), seriestype = [:path], label = "maxRoot: $(0.5)m",)