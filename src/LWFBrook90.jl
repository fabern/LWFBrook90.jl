module LWFBrook90

using OrdinaryDiffEq  # instead of loading the full DifferentialEquations
using DiffEqCallbacks # instead of loading the full DifferentialEquations
using RecipesBase

export read_inputData, discretize_soil
export define_LWFB90_p, define_LWFB90_u0, define_LWFB90_ODE
export KPT_SOILPAR_Mvg1d, KPT_SOILPAR_Ch1d
export RelativeDaysFloat2DateTime, plot_LWFBrook90

# on modules: https://discourse.julialang.org/t/large-programs-structuring-modules-include-such-that-to-increase-performance-and-readability/29102/5
include("module_CONSTANTS.jl");  # to bring into scope: using .CONSTANTS
include("module_KPT.jl");        using .KPT # using to bring exports into scope
include("module_WAT.jl");        using .WAT # using to bring exports into scope
include("module_SUN.jl");        # to bring into scope: using .SUN
include("module_PET.jl");        # to bring into scope: using .PET
include("module_SNO.jl");        # to bring into scope: using .SNO
include("module_EVP.jl");        # to bring into scope: using .SNO

include("func_read_inputData.jl")
include("func_discretize_soil_domain.jl")
include("func_DiffEq_definition_u0.jl")
include("func_DiffEq_definition_p.jl")
include("func_DiffEq_definition_cb.jl")
include("func_DiffEq_definition_f.jl")
include("func_DiffEq_definition_ode.jl")
include("func_MSB_functions.jl")

include("../example/BEA2016-reset-FALSE-input/func_run_example.jl") # defines RelativeDaysFloat2DateTime

function find_indices(depths_to_read_out_mm, solution)
    # depths and lower_boundaries must all be positive numbers
    @assert all(depths_to_read_out_mm .> 0)

    lower_boundaries = cumsum(solution.prob.p[1][1].p_THICK)

    # Only read out values that are within the simulation domain
    depths_to_read_out_mm = depths_to_read_out_mm[depths_to_read_out_mm .<= maximum(lower_boundaries)]

    idx_to_read_out = []
    for curr_depth_mm in depths_to_read_out_mm
        # idx_to_read_out = findfirst(curr_depth_mm .<= y)
        append!(idx_to_read_out, findfirst(curr_depth_mm .<= lower_boundaries))
    end
    return idx_to_read_out
end

function get_auxiliary_variables(solution)
    p_soil = solution.prob.p[1][1]
    NLAYER = p_soil.NLAYER
    idx_u_vector_amounts = solution.prob.p[1][4][4]
    u_SWATI = [solution.u[i][idx_u_vector_amounts] for i = 1:length(solution.u)]
    # (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
    #         LWFBrook90.KPT.derive_auxiliary_SOILVAR.(u_SWATI, Ref(p_soil)) # Ref fixes scalar argument for broadcasting "."
    u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK =
        (fill(NaN, (size(u_SWATI, 1), size(u_SWATI[1], 1))) for i in 1:5)
    for t in 1:length(u_SWATI)
        (u_aux_WETNES[t, :], u_aux_PSIM[t, :], u_aux_PSITI[t, :], u_aux_θ[t, :], p_fu_KK[t, :]) =
            LWFBrook90.KPT.derive_auxiliary_SOILVAR(u_SWATI[t], p_soil)
    end
    # returns arrays of dimenstion (t,z) where t is number of timesteps and z number of computational layers
    return (transpose(hcat(u_SWATI...)), u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK)
end

function get_θ(depths_to_read_out_mm, solution)
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        get_auxiliary_variables(solution)

    idx = find_indices(depths_to_read_out_mm, solution)

    return u_aux_θ[:, idx]
end

############################################################################################
############################################################################################
############################################################################################

@userplot PlotLWFBrook90

@recipe function f(h::PlotLWFBrook90)
    # This plot recipe generates a function LWFBrook90.plotlwfbrook90()
    # i.e. the function name is the type of the argument in lower case.

    # 0) parse input arguments
    # if length(h.args) != 1 # || !(typeof(h.args[1]) <: SciMLBase.ODESolution)
    #     error("Plot should be given two arguments of type ODESolution.  Got: $(typeof(h.args))")
    # end
    if length(h.args) == 1 # || !(typeof(h.args[1]) <: SciMLBase.ODESolution)
        # all okay, use default tick_function
        sol = h.args[1]
        tick_function = (x1, x2) -> [range(x1, x2; length = 4)][1]
    elseif length(h.args) == 2
        # all okay, tick_function was provided
        # check that second argument is indeed a function
        if isempty(methods(h.args[2]))
            error("Plot' second argument should be a function for the ticks.")
        end
        sol = h.args[1]
        tick_function = h.args[2]
    else
        error("Plot should be given two arguments of type ODESolution.  Got: $(typeof(h.args))")
    end

    # 1) data to plot

    # 1a) extract data from solution object `sol`
    # idx_u_vector_amounts       = sol.prob.p[1][4][4]
    # idx_u_vector_accumulators  = sol.prob.p[1][4][5]

    t_ref = sol.prob.p[2][15]
    x = RelativeDaysFloat2DateTime.(sol.t, t_ref)
    y = cumsum(sol.prob.p[1][1].p_THICK) - sol.prob.p[1][1].p_THICK / 2
    n = sol.prob.p[1][1].NLAYER
    y_centers = [0; cumsum(sol.prob.p[1][1].p_THICK)[1:(n-1)]] +
                sol.prob.p[1][1].p_THICK / 2
    # y_extended = [y; (maximum(y) .+ [10, 200])]
    y_soil_ticks = tick_function(0.0, round(maximum(cumsum(sol.prob.p[1][1].p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
    # y_labels   = [round.(y_soil_ticks; digits=0)]

    # row_PREC_amt = sol.prob.p[2][8].(sol.t)
    rows_SWAT_amt = sol[7 .+ (0:sol.prob.p[1][1].NLAYER-1),
        1,
        :] ./ sol.prob.p[1][1].p_THICK

    # row_NaN = fill(NaN, 1, length(x))

    # ###########################################################
    # # THIS PLOTTING RECIPE WILL REPRODUCE THE FOLLOWING PLOTS
    # # Plot 1
    # pl1 = plot(x,
    #     [sol_LWFBrook90[1, 1, :],
    #         sol_LWFBrook90[2, 1, :],
    #         sol_LWFBrook90[3, 1, :],
    #         sol_LWFBrook90[4, 1, :],
    #         sol_LWFBrook90[5, 1, :],
    #         sol_LWFBrook90[6, 1, :]];
    #     label = ["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])
    # # alternative: # pl1 = plot(x, sol_LWFBrook90[1, 1, :], label = "GWAT (mm)")
    # # alternative: # plot!(pl1, x, sol_LWFBrook90[2, 1, :], label = "INTS (mm)")
    # # alternative: # plot!(pl1, x, sol_LWFBrook90[3, 1, :], label = "INTR (mm)")
    # # alternative: # plot!(pl1, x, sol_LWFBrook90[4, 1, :], label = "SNOW (mm)")
    # # alternative: # plot!(pl1, x, sol_LWFBrook90[5, 1, :], label = "CC (MJ/m2)")
    # # alternative: # plot!(pl1, x, sol_LWFBrook90[6, 1, :], label = "SNOWLQ (mm)")
    # # alternative: # f_t2Dates = (t,x) -> (RelativeDaysFloat2DateTime(t, input_meteoveg_reference_date), x)
    # # alternative: # plot( f_t2Dates.(sol_LWFBrook90.t, sol_LWFBrook90[1,1,:]), label = "GWAT (mm)")
    # # alternative: # plot!(f_t2Dates.(sol_LWFBrook90.t, sol_LWFBrook90[2,1,:]), label = "INTS (mm)")
    # # alternative: # plot!(f_t2Dates.(sol_LWFBrook90.t, sol_LWFBrook90[3,1,:]), label = "INTR (mm)")
    # # alternative: # plot!(f_t2Dates.(sol_LWFBrook90.t, sol_LWFBrook90[4,1,:]), label = "SNOW (mm)")
    # # alternative: # plot!(f_t2Dates.(sol_LWFBrook90.t, sol_LWFBrook90[5,1,:]), label = "CC (MJ/m2)")
    # # alternative: # plot!(f_t2Dates.(sol_LWFBrook90.t, sol_LWFBrook90[6,1,:]), label = "SNOWLQ (mm)")
        # # Using dates (but not interpolated)
        # plot(example["solutionDates"],
        #     example["solution"][[1, 2, 3, 4, 5, 6], :]',
        #     label = ["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

        # # Using simple plot recipe that interpolates, but without dates
        # plot(example["solution"];
        #     vars = [1, 2, 3, 4, 5, 6],
        #     label = ["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

    # # Plot 2
    # # http://docs.juliaplots.org/latest/generated/gr/#gr-ref43
    # y = cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK)
    # y_extended = [y; (maximum(y) .+ [10, 200])]
    # y_soil_ticks = tick_function(0.0, round(maximum(cumsum(sol_LWFBrook90.prob.p[1][1].p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
    # # y_labels   = [round.(y_soil_ticks; digits=0)]
    # z = sol_LWFBrook90[7 .+ (0:sol_LWFBrook90.prob.p[1][1].NLAYER-1), 1, :] ./ sol_LWFBrook90.prob.p[1][1].p_THICK

    # pl2 = heatmap(x, y, z,
    #     yflip = true, yticks = y_soil_ticks,#(y_soil_ticks, y_labels),
    #     xlabel = "Date",
    #     ylabel = "Depth [mm]",
    #     colorbar_title = "θ [-]")
    # pl_final = plot(pl1, pl2; layout = (2, 1), link = :x)#@layout [a{0.75h}; b{0.8w} c])
    # ###########################################################

    # 2) set up the common arguments for all plots below
    # NOTE: --> sets attributes only when they don't already exist
    # NOTE: :=  sets attributes even when they already exist
    # link := :x #TODO(bernhard): link doesn't seem to work...
    layout --> (2, 1)
    # using layout because @layout is unsupported: https://github.com/JuliaPlots/RecipesBase.jl/issues/15
    # TODO(bernhard): find an easy way to do: l = @layout [grid(2, 1, heights=[0.2, 0.8]) a{0.055w}]
    # Possibly it needs to be provided when calling `plot_LWFBrook90_isotopes(... ; layout = ...)`
    # Then make sure to number the subplots correctly and include the colorbars from 3c)


    # 3) generate plots
    # NOTE: --> sets attributes only when they don't already exist
    # NOTE: :=  sets attributes even when they already exist

    # 3a) Scalar variables
    @series begin # pl1
        title := "Scalar state variables"
        seriestype := :path
        # # linecolor := :match
        label := ["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"]
        # fill_z := reshape(row_PREC_d18O, 1, :)
        # clims := clims_d18O
        # colorbar_title := "δ18O [mUr]"
        # colorbar := true_to_check_colorbar # TODO: define this once for all plots (except force it for colorbar plot representing common legend)
        # yguide := "PREC [mm]"
        # legend := false
        subplot := 1

        # and other arguments:
        x, transpose(sol[1:6, 1, :]) #transpose(sol[1:6, 1, :])
    end
    # # 3b) Heatmap (containing SWATI)
    @series begin # pl2
        title := "Soil (distributed state)"
        seriestype := :heatmap
        yflip := true
        yticks := y_soil_ticks #(y_ticks, y_labels)
        colorbar := true #true_to_check_colorbar
        yguide := "Depth [mm]"
        colorbar_title := "θ [-]"
        # clims := clims_d18O
        subplot := 2

        # and other arguments:
        # x, y_extended, rows_SWAT_amt
        # x, y, rows_SWAT_amt
        x, y_centers, rows_SWAT_amt
    end

    # TODO: edges of cells in heatmap are not entirely correct. Find a way to override heatmap()
    #       where we provide cell edges (n+1) instead of cell centers (n)
    # TODO: e.g. plots_heatmap_edges: @recipe function f(::Type{Val{:plots_heatmap_edges}}, xe, ye, z)
    # TODO: e.g. plots_heatmap_edges:     m, n = size(z.surf)
    # TODO: e.g. plots_heatmap_edges:     x_pts, y_pts = fill(NaN, 6 * m * n), fill(NaN, 6 * m * n)
    # TODO: e.g. plots_heatmap_edges:     fz = zeros(m * n)
    # TODO: e.g. plots_heatmap_edges:     for i in 1:m # y
    # TODO: e.g. plots_heatmap_edges:         for j in 1:n # x
    # TODO: e.g. plots_heatmap_edges:             k = (j - 1) * m + i
    # TODO: e.g. plots_heatmap_edges:             inds = (6 * (k - 1) + 1):(6 * k - 1)
    # TODO: e.g. plots_heatmap_edges:             x_pts[inds] .= [xe[j], xe[j + 1], xe[j + 1], xe[j], xe[j]]
    # TODO: e.g. plots_heatmap_edges:             y_pts[inds] .= [ye[i], ye[i], ye[i + 1], ye[i + 1], ye[i]]
    # TODO: e.g. plots_heatmap_edges:             fz[k] = z.surf[i, j]
    # TODO: e.g. plots_heatmap_edges:         end
    # TODO: e.g. plots_heatmap_edges:     end
    # TODO: e.g. plots_heatmap_edges:     ensure_gradient!(plotattributes, :fillcolor, :fillalpha)
    # TODO: e.g. plots_heatmap_edges:     fill_z := fz
    # TODO: e.g. plots_heatmap_edges:     line_z := fz
    # TODO: e.g. plots_heatmap_edges:     x := x_pts
    # TODO: e.g. plots_heatmap_edges:     y := y_pts
    # TODO: e.g. plots_heatmap_edges:     z := nothing
    # TODO: e.g. plots_heatmap_edges:     seriestype := :shape
    # TODO: e.g. plots_heatmap_edges:     label := ""
    # TODO: e.g. plots_heatmap_edges:     widen --> false
    # TODO: e.g. plots_heatmap_edges:     ()
    # TODO: e.g. plots_heatmap_edges: end
    # TODO: e.g. plots_heatmap_edges: @deps plots_heatmap_edges shape
    # TODO: e.g. plots_heatmap_edges: @shorthands plots_heatmap_edges
    # TODO: e.g. plots_heatmap_edges:
    # TODO: e.g. plots_heatmap_edges: Plots.heatmap(x[1:100], y_centers, z[:,1:100])
    # TODO: e.g. plots_heatmap_edges: Plots.heatmap(x[1:100], y_centers, z[:,1:100])
    # TODO: e.g. plots_heatmap_edges: plot(t = :heatmap, x[1:50], y_centers, z[:,1:50]) # works
    # TODO: e.g. plots_heatmap_edges: plot(t = :plots_heatmap, x[1:50], y_centers, z[:,1:50]) # doesn't work
    # TODO: e.g. plots_heatmap_edges: plot(t = :plots_heatmap_edges, x[1:50], y_centers, z[:,1:50]) # doesn't work either


    # 3c) Colorbar
    # ....
    # TODO: this code needs to reproduce below steps:
    # pl_colorbar_δ18O = plot([0,0], [0,1], zcolor=[0,1], t=:scatter, xlims=(1,1.1), # set xlims that no marker is shown
    #                         clims=clims_d18O, colorbar_title="δ180 [mUr]",
    #                         grid=false, showaxis=false, ticks=false, label=false);
    # pl_colorbar_δ2H  = plot([0,0], [0,1], zcolor=[0,1], t=:scatter, xlims=(1,1.1), # set xlims that no marker is shown
    #                         clims=clims_d18O, colorbar_title="δ2H [mUr]",
    #                         grid=false, showaxis=false, ticks=false, label=false);

    # l = @layout [grid(2, 1, heights=[0.2, 0.8]) a{0.055w}]
    # pl_final_δ18O = plot(ts_PREC_δ18O, pl_δ18O, pl_colorbar_δ18O,
    #     clims = clims_d18O,
    #     layout = l, link = :x);
    # l = @layout [grid(2, 1, heights=[0.2, 0.8]) a{0.055w}]
    # pl_final_δ2H = plot(ts_PREC_δ2H,  pl_δ2H, pl_colorbar_δ2H,
    #     clims = clims_d2H,
    #     layout = l, link = :x);
end

end # module
