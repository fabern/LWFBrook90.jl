function get_RWU_centroid(rows_RWU_mmDay, y_center)
    # old solution, that gives bad values when some RWU is negative:
    RWU_percent = rows_RWU_mmDay ./ sum(rows_RWU_mmDay; dims = 1)
    RWUcentroidLabel = "mean RWU depth"
    if (any(RWU_percent .< 0))
        @warn "Some root water outfluxes detected. Centroid of RWU is  based only on uptakes."
        rows_RWU_mmDay_onlyUptake = ifelse.(rows_RWU_mmDay.>0,rows_RWU_mmDay, 0)
        RWU_percent_onlyUptake = rows_RWU_mmDay_onlyUptake ./ sum(rows_RWU_mmDay_onlyUptake; dims = 1)
        RWU_percent = RWU_percent_onlyUptake
        # RWU_percent = min.(0,rows_RWU_mmDay) ./ sum(min.(0,rows_RWU_mmDay); dims = 1)

        RWUcentroidLabel = "mean RWU depth\n(based on uptake only)"
    end

    return (row_RWU_centroid_mm = sum(RWU_percent .* y_center; dims=1),
            RWUcentroidLabel)
end

"""
    plotamounts(simulation::DiscretizedSPAC)
    plotamounts(simulation::DiscretizedSPAC, compartments::Symbol)
    plotamounts(simulation::DiscretizedSPAC, compartments::Symbol, RWUcentroid::Symbol)

Plots the amount results of a SPAC Simulation. By default both above and belowground.
The user can override this with the second argument isotope as one of `:aboveground`, `:belowground`, or `:above_and_belowground`.
RWUcentroid can have values of either `:dontShowRWUcentroid` or `:showRWUcentroid`.
"""
@userplot PlotAmounts
@recipe function f(plam::PlotAmounts)
    # 0) parse input arguments

    if length(plam.args) == 3
        simulation = plam.args[1]
        compartments = plam.args[2]
        RWUcentroid = plam.args[3]
    elseif length(plam.args) == 2
        simulation = plam.args[1]
        compartments = plam.args[2]
        RWUcentroid = :dontShowRWUcentroid
    elseif length(plam.args) == 1
        simulation = plam.args[1]
        compartments = :above_and_belowground
        RWUcentroid= :dontShowRWUcentroid
    else
        error("plotamounts requires an unnamed first argument of type DiscretizedSPAC, and optional unnamed second/third arguments (:aboveground, :belowground, or :above_and_belowground) and (:dontShowRWUcentroid, :showRWUcentroid). Other arguments to plot() should be separated by `;`.")
    end
    if (compartments != :above_and_belowground)
        error("TODO: currently selecting a subset of compartments are not supported by plotamounts. Please only do: `plotamounts(simulation)` or `plotamounts(simulation, :above_and_belowground, :showRWUcentroid)`")
    end
    if !(RWUcentroid == :dontShowRWUcentroid || RWUcentroid == :showRWUcentroid)
        error("Third unnamed argument to plotamounts should be one of (:dontShowRWUcentroid, :showRWUcentroid). Got: $(RWUcentroid)")
    end
    if !(compartments == :aboveground || compartments == :belowground || compartments == :above_and_belowground)
        error("Second unnamed argument to plotamounts should be one of (:above_and_belowground, :aboveground, or :belowground). Got: $( compartments)")
    end
    if !(simulation isa DiscretizedSPAC)
        error("First unnamed argument to plotamounts should be of type DiscretizedSPAC. Got: $(typeof(simulation))")
    end
    if isnothing(simulation.ODESolution)
        error("plotamounts requires a solved system. Please `simulate!()` the DiscretizedSPAC first.")
    end

    # 1) prepare data to plot
    solu = simulation.ODESolution

    # 1a) extract data from solution object `solu`
        # # Two ways to extract data from soil object: using `[]` or `()`
        # u_SWATI = reduce(hcat, [solu[t_idx].SWATI.mm  for t_idx = eachindex(solu)])
        # u_SWATI = reduce(hcat, [solu(t_days).SWATI.mm for t_days = days_to_read_out_d])

        # days_to_read_out_d decides which points to use:
        # days_to_read_out_d = nothing # read out all simulation steps
        days_to_read_out_d = :daily  # read out only daily values
        if isnothing(days_to_read_out_d)
            days_to_read_out_d = solu.t # warning this can
            error("Too many output times requested. This easily freezes the program...")
        elseif days_to_read_out_d == :daily
            days_to_read_out_d = unique(round.(solu.t))
        end
        t_ref = solu.prob.p.REFERENCE_DATE
        x = RelativeDaysFloat2DateTime.(days_to_read_out_d, t_ref);
        y_center = cumsum(solu.prob.p.p_soil.p_THICK) - solu.prob.p.p_soil.p_THICK/2

    # Some hardcoded options:
    xlimits = RelativeDaysFloat2DateTime.(solu.prob.tspan, t_ref)
    # color_scheme = :default # https://docs.juliaplots.org/latest/generated/colorschemes
    # color_scheme = :heat
    color_scheme = :blues

    true_to_check_colorbar = true; # set this flag to false for final plot, true for debugging.
    tick_function = (x1, x2) -> PlotUtils.optimize_ticks(x1, x2; k_min = 4)
    y_soil_ticks = tick_function(0., round(maximum(cumsum(solu.prob.p.p_soil.p_THICK))))[1]
    y_ticks    = [-500;       -300;       -200;       -100;          y_soil_ticks;             (maximum(cumsum(solu.prob.p.p_soil.p_THICK)) .+ [    100;      250;     400])]
    y_labels   = ["PREC";   "INTS";     "INTR";     "SNOW";     round.(y_soil_ticks; digits=0);                                                "GWAT";    "RWU";     "XYLEM"]


    # Results
    rows_SWAT_amt  = reduce(hcat, [solu(t).SWATI.mm    for t in days_to_read_out_d]) ./ solu.prob.p.p_soil.p_THICK
    rows_RWU_mmDay = reduce(hcat, [solu(t).TRANI.mmday for t in days_to_read_out_d])
    row_NaN       = fill(NaN, 1,length(x))
    row_PREC_amt = reshape(solu.prob.p.p_PREC.(days_to_read_out_d), 1, :)
    col_INTS_amt = [solu(t).INTS.mm for t in days_to_read_out_d]
    col_INTR_amt = [solu(t).INTR.mm   for t in days_to_read_out_d]
    col_SNOW_amt = [solu(t).SNOW.mm   for t in days_to_read_out_d]
    col_GWAT_amt = [solu(t).GWAT.mm   for t in days_to_read_out_d]
    col_RWU_amt  = [solu(t).RWU.mmday for t in days_to_read_out_d]
    col_XYL_amt  = [solu(t).XYLEM.mm  for t in days_to_read_out_d]

    if (RWUcentroid == :showRWUcentroid)
        row_RWU_centroid_mm, RWUcentroidLabel = get_RWU_centroid(rows_RWU_mmDay, y_center)
    end

    # reduce how deep to plot soil:
    soil_discr_to_plot = simulation.parametrizedSPAC.soil_discretization.df#[simulation.soil_discretization.Lower_m .>= -0.2, :]

    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        get_auxiliary_variables(simulation, days_to_read_out_d = days_to_read_out_d);

    # compute total amount of soil water: a) in whole domain, b) per physical soil layer
    # a)
    col_sumSWATI_amt = sum(u_SWATI; dims = [1])[:]

    # b)
    # sum up SWATI to get total SWAT for each of the physical soil horizons
    # terminology: soil horizon refers to a physical horizon, discretization refers to the computational layers/cells
    # soil_horizons = unique(soil_discr_to_plot, :HorizonNr)
    soil_horizons = unique(soil_discr_to_plot, :HorizonNr)[:,[:HorizonNr, :Horizon_Upper_m, :Horizon_Lower_m]]
    cols_sumSWATIperLayer_amt =
        [sum(u_SWATI[soil_discr_to_plot.HorizonNr .== horizon.HorizonNr, :]; dims = [1])'
            for horizon in eachrow(soil_horizons)]
    horizon_labels =
        ["Horizon-$(horizon.HorizonNr): ($(round(horizon.Horizon_Upper_m;digits=2)) to $(round(horizon.Horizon_Lower_m;digits=2)) m)"
        for horizon in eachrow(soil_horizons)]

    # 2) set up the common arguments for all plots below
    # if (compartments == :aboveground || compartments == :belowground)
    #     layout --> RecipesBase.grid(2, 1, heights=[0.2 ,0.8])           #layout --> (2,1)
    #     size --> (1000,700)
    #     idx_d18O_PREC = 1
    #     idx_d18O_SWAT = 2
    #     idx_d2H_PREC  = 1
    #     idx_d2H_SWAT  = 2
    # elseif (compartments == :above_and_belowground)
        # lay = RecipesBase.@layout [RecipesBase.grid(11, 1)]
        lay = RecipesBase.@layout([ °{0.80w} _ ; #1 have no colorbar
                                    °{0.80w} _ ; #2 have no colorbar
                                    °{0.80w} _ ; #3 have no colorbar
                                    ° ; #4
                                    ° ; #5
                                    ° ; #6
                                    ° ;]) #7
        layout --> lay

        size --> (1000,2100)
    #     idx_d18O_PREC = 1
    #     idx_d18O_SWAT = 2
    #     idx_d2H_PREC  = 3
    #     idx_d2H_SWAT  = 4
    # else
    #     # nothing
    # end
    dpi --> 300
    xlim --> xlimits
    leftmargin --> 15mm
    rightmargin --> 15mm

    # # 3) generate plots
    # # NOTE: --> sets attributes only when they don't already exist
    # # NOTE: :=  sets attributes even when they already exist

    # 3a) Precipitation
    @series begin
        # title := "Precipitation"
        seriestype := :bar
        # linecolor := :match
        # line_z := reshape(row_PREC_d18O,1,:)
        # fill_z := reshape(row_PREC_d18O,1,:)
        # clims := clims_d18O
        # colorbar_title := "δ18O [‰]"
        # colorbar := true_to_check_colorbar # TODO: define this once for all plots (except force it for colorbar plot representing common legend)
        yguide := "PREC [mm]"
        legend := false
        subplot := 1
        # and other arguments:
        reshape(x,1,:), reshape(row_PREC_amt,1,:)
    end
    if (compartments == :aboveground || compartments == :above_and_belowground)
        # plot(x, [col_INTS_amt col_INTR_amt col_SNOW_amt col_GWAT_amt col_XYL_amt],
        #      labels = ["INTS" "INTR" "SNOW" "GWAT" "XYL"], ylab = "Amount [mm]")
        @series begin
            title := "Aboveground"
            labels := ["INTS" "INTR" "SNOW" "GWAT" "XYL"]
            ylab := "Amount [mm]"
            seriestype := :line
            # legend := false
            subplot := 2
            bg_legend --> colorant"rgba(100%,100%,100%,0.8)"; legend := :topright
            # and other arguments:
            x, [col_INTS_amt col_INTR_amt col_SNOW_amt col_GWAT_amt col_XYL_amt]
        end
    end
    if (compartments == :belowground || compartments == :above_and_belowground)
        @series begin
            title := "Belowground"
            labels := ["total SWAT" reduce(hcat, horizon_labels)]
            ylab := "Amount [mm]"
            seriestype := :line
            subplot := 3
            bg_legend --> colorant"rgba(100%,100%,100%,0.8)"; legend := :topright
            # and other arguments:
            x, hcat(col_sumSWATI_amt[:], reduce(hcat, cols_sumSWATIperLayer_amt))
        end

        # rows_SWAT_amt0 = u_aux_θ
        rows_SWAT_amt1 = u_SWATI ./solu.prob.p.p_soil.p_THICK  # mm per mm of soil thickness
        rows_SWAT_amt2 = u_SWATI ./ solu.prob.p.p_soil.p_THICK ./ (1 .- solu.prob.p.p_soil.p_STONEF)
        # mm per mm of fine soil thickness (assuming gravel fraction contains no water)
        rows_ψₘpF  = log10.(u_aux_PSIM  .* -10) #(from kPa to log10(hPa))
        rows_ψₜₒₜpF = log10.(u_aux_PSITI .* -10) #(from kPa to log10(hPa))
        # rows_ψₘpF= log10.(u_aux_PSIM  .* -10 + 0.1) #(from kPa to log10(hPa) and small shift to start at pF = -1 instead of pF = -∞)

        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar; # clims := clims_d2H
        #     yguide := "Depth [mm]"; colorbar_title := "ψₘ [kPa]"
        #     c := cgrad(color_scheme,  rev = false)
        #     subplot := 4
        #     # and other arguments:
        #     x, y_center, u_aux_PSIM;
        # end
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar; # clims := clims_d2H
        #     yguide := "Depth [mm]"; colorbar_title := "ψₜ [kPa]"
        #     c := cgrad(color_scheme,  rev = true)
        #     subplot := 5
        #     # and other arguments:
        #     x, y_center, u_aux_PSITI;
        # end
        @series begin
            seriestype := :heatmap
            yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
            colorbar := true_to_check_colorbar; # clims := clims_d2H
            yguide := "Depth [mm]"; colorbar_title := "pF = log₁₀(-ψₘ hPa)"
            c := cgrad(color_scheme,  rev = true)
            subplot := 4
            # and other arguments:
            x, y_center, rows_ψₘpF
        end
        if (RWUcentroid == :showRWUcentroid)
            @series begin
                #plot!(x, row_RWU_centroid_mm', yflip=true, color=:white, label = "")
                color := :white
                label := RWUcentroidLabel
                #bg_legend --> colorant"rgba(100%,100%,100%,0.8)";legend := :bottomright
                bg_legend --> colorant"rgba(100%,100%,100%,0.0)"; legend := :bottomright; fg_legend --> :transparent; legendfontcolor := :white
                yflip := true; yticks := y_soil_ticks
                yguide := "Depth [mm]"; colorbar_title := "pF = log₁₀(-ψₘ hPa)"
                subplot := 4
                x, row_RWU_centroid_mm'
            end
        end
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar; # clims := clims_d2H
        #     yguide := "Depth [mm]"; colorbar_title := "pF = log₁₀(-ψₜₒₜ hPa)" #colorbar_title := "pF = \nlog₁₀(-ψ hPa)"
        #     c := cgrad(color_scheme,  rev = true)
        #     subplot := 5
        #     # and other arguments:
        #     x, y_center, rows_ψₜₒₜpF
        # end
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar; # clims := clims_d2H
        #     yguide := "Depth [mm]"; colorbar_title := "SWATI [mm]"
        #     c := cgrad(color_scheme,  rev = false)
        #     subplot := 5
        #     # and other arguments:
        #     x, y_center, u_SWATI;                 # deactivated u_SWATI as it is resolution dependent!
        # end
        @series begin
            seriestype := :heatmap
            yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
            colorbar := true_to_check_colorbar; # clims := clims_d2H
            yguide := "Depth [mm]"; colorbar_title := "Wetness [-]"
            c := cgrad(color_scheme,  rev = false)
            subplot := 5
            # and other arguments:
            x, y_center, u_aux_WETNES;
        end
        @series begin
            seriestype := :heatmap
            yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
            colorbar := true_to_check_colorbar; # clims := clims_d2H
            yguide := "Depth [mm]"; colorbar_title := "θ [m3/m3]\n(of fine soil volume)" #colorbar_title := "θ [m3/m3]\n(fine soil)" # "θ [m3/m3]"
            c := cgrad(color_scheme,  rev = false)
            subplot := 6
            # and other arguments:
            x, y_center, u_aux_θ;
        end
        @series begin
            seriestype := :heatmap
            yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
            colorbar := true_to_check_colorbar; # clims := clims_d2H
            yguide := "Depth [mm]"; colorbar_title := "θ [m3/m3]\n(of total volume, incl stonef)" #colorbar_title := "θ [-]\n(total, incl stonef)"
            c := cgrad(color_scheme,  rev = false)
            subplot := 7
            # and other arguments:
            x, y_center, rows_SWAT_amt1;
        end
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar
        #     yguide := "Depth [mm]"; colorbar_title := "θ [-] (fine soil 2)"#colorbar_title := "θ [-]\n(fine soil 2)"
        #     c := cgrad(color_scheme,  rev = false)
        #     subplot := 8
        #     # and other arguments:
        #     x, y_center, rows_SWAT_amt2           # deactivated: as it was same as `x, y_center, u_aux_θ;`
        # end
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar; # clims := clims_d2H
        #     yguide := "Depth [mm]"; colorbar_title := "K [mm/day]"
        #     c := cgrad(color_scheme,  rev = false)
        #     subplot := 8
        #     # and other arguments:
        #     x, y_center, p_fu_KK;
        #     # x, y_center, log10.(p_fu_KK);
        # end

        # TODO: edges of cells in heatmap are not entirely correct. Find a way to override heatmap()
        #       where we provide cell edges (n+1) instead of cell centers (n)
        #       e.g. plots_heatmap_edges: @recipe function f(::Type{Val{:plots_heatmap_edges}}, xe, ye, z)
        #       e.g. plots_heatmap_edges:     m, n = size(z.surf)
        #       e.g. plots_heatmap_edges:     x_pts, y_pts = fill(NaN, 6 * m * n), fill(NaN, 6 * m * n)
        #       e.g. plots_heatmap_edges:     fz = zeros(m * n)
        #       e.g. plots_heatmap_edges:     for i in 1:m # y
        #       e.g. plots_heatmap_edges:         for j in 1:n # x
        #       e.g. plots_heatmap_edges:             k = (j - 1) * m + i
        #       e.g. plots_heatmap_edges:             inds = (6 * (k - 1) + 1):(6 * k - 1)
        #       e.g. plots_heatmap_edges:             x_pts[inds] .= [xe[j], xe[j + 1], xe[j + 1], xe[j], xe[j]]
        #       e.g. plots_heatmap_edges:             y_pts[inds] .= [ye[i], ye[i], ye[i + 1], ye[i + 1], ye[i]]
        #       e.g. plots_heatmap_edges:             fz[k] = z.surf[i, j]
        #       e.g. plots_heatmap_edges:         end
        #       e.g. plots_heatmap_edges:     end
        #       e.g. plots_heatmap_edges:     ensure_gradient!(plotattributes, :fillcolor, :fillalpha)
        #       e.g. plots_heatmap_edges:     fill_z := fz
        #       e.g. plots_heatmap_edges:     line_z := fz
        #       e.g. plots_heatmap_edges:     x := x_pts
        #       e.g. plots_heatmap_edges:     y := y_pts
        #       e.g. plots_heatmap_edges:     z := nothing
        #       e.g. plots_heatmap_edges:     seriestype := :shape
        #       e.g. plots_heatmap_edges:     label := ""
        #       e.g. plots_heatmap_edges:     widen --> false
        #       e.g. plots_heatmap_edges:     ()
        #       e.g. plots_heatmap_edges: end
        #       e.g. plots_heatmap_edges: @deps plots_heatmap_edges shape
        #       e.g. plots_heatmap_edges: @shorthands plots_heatmap_edges
        #       e.g. plots_heatmap_edges:
        #       e.g. plots_heatmap_edges: Plots.heatmap(x[1:100], y_center, z[:,1:100])
        #       e.g. plots_heatmap_edges: Plots.heatmap(x[1:100], y_center, z[:,1:100])
        #       e.g. plots_heatmap_edges: plot(t = :heatmap, x[1:50], y_center, z[:,1:50]) # works
        #       e.g. plots_heatmap_edges: plot(t = :plots_heatmap, x[1:50], y_center, z[:,1:50]) # doesn't work
        #       e.g. plots_heatmap_edges: plot(t = :plots_heatmap_edges, x[1:50], y_center, z[:,1:50]) # doesn't work either

    end
end

"""
    plotisotopes(simulation::DiscretizedSPAC)
    plotisotopes(simulation::DiscretizedSPAC, isotope::Symbol)
    plotisotopes(simulation::DiscretizedSPAC, isotope::Symbol, RWUcentroid::Symbol)

Plots the isotope results of a SPAC Simulation. By default both δ18O and δ2H.
The user can override this with the second argument isotope as one of `:d18O`, `:d2H`, or `:d18O_and_d2H`.
RWUcentroid can have values of either `:dontShowRWUcentroid` or `:showRWUcentroid`.
"""
@userplot PlotIsotopes
@recipe function f(pliso::PlotIsotopes)
    # 0) parse input arguments
    if length(pliso.args) == 3
        simulation = pliso.args[1]
        isotope    = pliso.args[2]
        RWUcentroid= pliso.args[3]
    elseif length(pliso.args) == 2
        simulation = pliso.args[1]
        isotope    = pliso.args[2]
        RWUcentroid= :dontShowRWUcentroid
    elseif length(pliso.args) == 1
        simulation = pliso.args[1]
        isotope    = :d18O_and_d2H
        RWUcentroid= :dontShowRWUcentroid
    else
        error("plotisotopes requires an unnamed first argument of type DiscretizedSPAC, and optional unnamed second/third arguments (:d18O_and_d2H, :d18O, or :d2H) and (:dontShowRWUcentroid, :showRWUcentroid). Other arguments to plot() should be separated by `;`.")
    end
    if !(RWUcentroid == :dontShowRWUcentroid || RWUcentroid == :showRWUcentroid)
        error("Third unnamed argument to plotisotopes should be one of (:dontShowRWUcentroid, :showRWUcentroid). Got: $(RWUcentroid)")
    end
    if !(isotope == :d18O || isotope == :d2H || isotope == :d18O_and_d2H)
        error("Second unnamed argument to plotisotopes should be one of (:d18O_and_d2H, :d18O, or :d2H). Got: $(isotope)")
    end
    if !(simulation isa DiscretizedSPAC)
        error("First unnamed argument to plotisotopes should be of type DiscretizedSPAC. Got: $(typeof(simulation))")
    end
    if isnothing(simulation.ODESolution)
        error("plotisotopes requires a solved system. Please `simulate!()` the DiscretizedSPAC first.")
    end

    # 1) prepare data to plot
    solu = simulation.ODESolution

    # 1a) extract data from solution object `solu`
        # # Two ways to extract data from soil object: using `[]` or `()`
        # u_SWATI = reduce(hcat, [solu[t_idx].SWATI.mm  for t_idx = eachindex(solu)])
        # u_SWATI = reduce(hcat, [solu(t_days).SWATI.mm for t_days = days_to_read_out_d])

        # days_to_read_out_d decides which points to use:
        # days_to_read_out_d = nothing # read out all simulation steps
        days_to_read_out_d = :daily  # read out only daily values
        if isnothing(days_to_read_out_d)
            days_to_read_out_d = solu.t # warning this can
            error("Too many output times requested. This easily freezes the program...")
        elseif days_to_read_out_d == :daily
            days_to_read_out_d = unique(round.(solu.t))
        end
        t_ref = solu.prob.p.REFERENCE_DATE
        x = RelativeDaysFloat2DateTime.(days_to_read_out_d, t_ref);
        y_center = cumsum(solu.prob.p.p_soil.p_THICK) - solu.prob.p.p_soil.p_THICK/2

    simulate_isotopes = solu.prob.p.simulate_isotopes
    @assert simulate_isotopes "Provided DiscretizedSPAC() did not simulate isotopes"

    # Some hardcoded options:
    xlimits = RelativeDaysFloat2DateTime.(solu.prob.tspan, t_ref)

    clims_d18O = (-16, -6)
    clims_d2H  = (-125, -40)
    true_to_check_colorbar = false; # set this flag to false for final plot, true for debugging.
    tick_function = (x1, x2) -> PlotUtils.optimize_ticks(x1, x2; k_min = 4)
    # color_scheme = :default # https://docs.juliaplots.org/latest/generated/colorschemes
    # color_scheme = :blues
    color_scheme = :heat

    # Results
    row_PREC_amt = solu.prob.p.p_PREC.(days_to_read_out_d)
    # rows_SWAT_amt = solu[solu.prob.p.row_idx_SWATI, 1, :]./solu.prob.p.p_soil.p_THICK;
    rows_SWAT_amt  = reduce(hcat, [solu(t).SWATI.mm   for t in days_to_read_out_d]) ./ solu.prob.p.p_soil.p_THICK
    rows_SWAT_d18O = reduce(hcat, [solu(t).SWATI.d18O for t in days_to_read_out_d])
    rows_SWAT_d2H  = reduce(hcat, [solu(t).SWATI.d2H  for t in days_to_read_out_d])
    row_NaN       = fill(NaN, 1,length(x))
    row_PREC_d18O = reshape(solu.prob.p.p_δ18O_PREC.(days_to_read_out_d), 1, :)
    row_INTS_d18O = reduce(hcat, [solu(t).INTS.d18O for t in days_to_read_out_d])
    row_INTR_d18O = reduce(hcat, [solu(t).INTR.d18O for t in days_to_read_out_d])
    row_SNOW_d18O = reduce(hcat, [solu(t).SNOW.d18O for t in days_to_read_out_d])
    row_GWAT_d18O = reduce(hcat, [solu(t).GWAT.d18O for t in days_to_read_out_d])
    row_RWU_d18O  = reduce(hcat, [solu(t).RWU.d18O for t in days_to_read_out_d])
    row_XYL_d18O  = reduce(hcat, [solu(t).XYLEM.d18O for t in days_to_read_out_d])
    row_PREC_d2H  = reshape(solu.prob.p.p_δ2H_PREC.(days_to_read_out_d), 1, :)
    row_INTS_d2H  = reduce(hcat, [solu(t).INTS.d2H for t in days_to_read_out_d])
    row_INTR_d2H  = reduce(hcat, [solu(t).INTR.d2H for t in days_to_read_out_d])
    row_SNOW_d2H  = reduce(hcat, [solu(t).SNOW.d2H for t in days_to_read_out_d])
    row_GWAT_d2H  = reduce(hcat, [solu(t).GWAT.d2H for t in days_to_read_out_d])
    row_RWU_d2H   = reduce(hcat, [solu(t).RWU.d2H for t in days_to_read_out_d])
    row_XYL_d2H   = reduce(hcat, [solu(t).XYLEM.d2H for t in days_to_read_out_d])

    rows_RWU_mmDay  = reduce(hcat, [solu(t).TRANI.mmday   for t in days_to_read_out_d])
    if (RWUcentroid == :showRWUcentroid)
        row_RWU_centroid_mm, RWUcentroidLabel = get_RWU_centroid(rows_RWU_mmDay, y_center)
    end

    # # 1b) define some plot arguments based on the extracted data
    # # color scheme:
    # all_d18O_values = [rows_SWAT_d18O; row_PREC_d18O; row_INTS_d18O; row_INTR_d18O; row_SNOW_d18O; row_GWAT_d18O; row_RWU_d18O][:]
    # all_d2H_values  = [rows_SWAT_d2H ; row_PREC_d2H ; row_INTS_d2H ; row_INTR_d2H ; row_SNOW_d2H ; row_GWAT_d2H ; row_RWU_d2H ][:]
    # if clims_d18O == :auto
    #     clims_d18O = extrema(filter(!isnan, all_d18O_values))
    # else
    #     clims_d18O = (-16, -6)
    # end
    # if clims_d2H == :auto
    #     clims_d2H = clims_d2H  = extrema(filter(!isnan, all_d2H_values))
    # else
    #     clims_d2H  = (-125, -40)
    # end
    # true_to_check_colorbar = true; # set this flag to false for final plot, true for debugging.

    # 2) set up the common arguments for all plots below
    if (isotope == :d18O || isotope == :d2H)
        lay = RecipesBase.@layout([ °{0.91w} _ ;   #1 have no colorbar
                                    °          ;]) #2
        layout --> lay
        size --> (1000,700)
        idx_d18O_PREC = 1
        idx_d18O_SWAT = 2
        idx_d2H_PREC  = 1
        idx_d2H_SWAT  = 2
    elseif (isotope == :d18O_and_d2H)
        lay = RecipesBase.@layout([ °{0.86w} _ ;   #1 have no colorbar
                                    °          ;   #2
                                    °{0.86w} _ ;   #3 have no colorbar
                                    °          ;]) #4
        layout --> lay
        size --> (1000,1400)
        idx_d18O_PREC = 1
        idx_d18O_SWAT = 2
        idx_d2H_PREC  = 3
        idx_d2H_SWAT  = 4
    else
        # nothing
    end
    dpi --> 300
    xlim --> xlimits
    leftmargin --> 15mm

    # # 3) generate plots
    # # NOTE: --> sets attributes only when they don't already exist
    # # NOTE: :=  sets attributes even when they already exist

    # 3a) Precipitation
    # We reproduce the following plot:
    # Workaround for barplot # https://github.com/JuliaPlots/Plots.jl/issues/3880
    # ts_PREC_δ18O = plot(reshape(x,1,:), reshape(row_PREC_amt,1,:), t = :bar, line_z = reshape(row_PREC_d18O,1,:), fill_z = reshape(row_PREC_d18O,1,:),
    #                 clims=clims_d18O, colorbar = true_to_check_colorbar, legend = false,
    #                 ylabel = "PREC [mm]"); # TODO(bernhard): make this work with a barplot
    # ts_PREC_δ2H = ...
    if (isotope == :d18O || isotope == :d18O_and_d2H)
        @series begin # ts_PREC_δ18O
            title --> "δ18O"
            seriestype := :bar
            # linecolor := :match
            line_z := reshape(row_PREC_d18O,1,:)
            fill_z := reshape(row_PREC_d18O,1,:)
            clims := clims_d18O
            colorbar_title := "δ18O [‰]"
            yguide := "PREC [mm]"
            colorbar := true_to_check_colorbar #https://stackoverflow.com/a/59257011
            legend := false;
            subplot := idx_d18O_PREC

            # and other arguments:
            reshape(x,1,:), reshape(row_PREC_amt,1,:)
        end
    end
    if (isotope == :d2H || isotope == :d18O_and_d2H)
        @series begin # ts_PREC_δ2H
            title --> "δ2H"
            seriestype := :bar
            # linecolor := :match
            line_z := reshape(row_PREC_d2H,1,:)
            fill_z := reshape(row_PREC_d2H,1,:)
            clims := clims_d2H
            colorbar_title := "δ2H [‰]"
            yguide := "PREC [mm]"
            colorbar := true_to_check_colorbar #https://stackoverflow.com/a/59257011
            legend := false;
            subplot := idx_d2H_PREC

            # and other arguments:
            reshape(x,1,:), reshape(row_PREC_amt,1,:)
        end
    end

    # 3b) Heatmap (containing SWATI and other compartments)
    # y_labels   = ["INTS"; ""; "INTR"; ""; "SNOW"; ""; round.(y_center); "";             "GWAT"]
    # y_soil_ticks = optimize_ticks(extrema(y_center)...; k_min = 4)[1]
    # y_soil_ticks = optimize_ticks(0., round(maximum(cumsum(solu.prob.p.p_soil.p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
    y_extended = [-500; -350; -300; -250; -200; -150; -100; -50;         y_center;             (maximum(cumsum(solu.prob.p.p_soil.p_THICK)) .+ [50; 100; 150; 250; 300;400])]
    y_soil_ticks = tick_function(0., round(maximum(cumsum(solu.prob.p.p_soil.p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
    y_ticks    = [-500;       -300;       -200;       -100;          y_soil_ticks;             (maximum(cumsum(solu.prob.p.p_soil.p_THICK)) .+ [    100;      250;     400])]
    y_labels   = ["PREC";   "INTS";     "INTR";     "SNOW";     round.(y_soil_ticks; digits=0);                                                "GWAT";    "RWU";     "XYLEM"]
    z2_extended = [row_PREC_d18O; row_NaN; row_INTS_d18O; row_NaN; row_INTR_d18O; row_NaN; row_SNOW_d18O; row_NaN; rows_SWAT_d18O; row_NaN; row_GWAT_d18O; row_NaN; row_RWU_d18O; row_NaN; row_XYL_d18O]
    z3_extended = [row_PREC_d2H;  row_NaN; row_INTS_d2H;  row_NaN; row_INTR_d2H;  row_NaN; row_SNOW_d2H;  row_NaN; rows_SWAT_d2H;  row_NaN; row_GWAT_d2H;  row_NaN; row_RWU_d2H;  row_NaN; row_XYL_d2H ]

    if (isotope == :d18O || isotope == :d18O_and_d2H)
        @series begin # pl_δ18O
            seriestype := :heatmap
            yflip := true
            yticks := (y_ticks, y_labels)
            yguide := "Depth [mm]"
            colorbar_title := "δ18O [‰]"
            clims := clims_d18O
            colorbar := true # colorbar := true_to_check_colorbar
            subplot := idx_d18O_SWAT

            # and other arguments:
            x, y_extended, z2_extended;
        end
        if (RWUcentroid == :showRWUcentroid)
            @series begin
                color := :white
                label := "mean RWU depth"
                bg_legend --> colorant"rgba(100%,100%,100%,0.0)"; legend := :right;       fg_legend --> :transparent; legendfontcolor := :white
                # bg_legend --> colorant"rgba(100%,100%,100%,0.0)"; legend := :bottomright; fg_legend --> :transparent; legendfontcolor := :white
                yflip := true; yticks := (y_ticks, y_labels)
                yguide := "Depth [mm]"; colorbar_title := "δ18O [‰]"
                subplot := idx_d18O_SWAT
                x, row_RWU_centroid_mm'
            end
        end
    end

    if (isotope == :d2H || isotope == :d18O_and_d2H)
        @series begin # pl_δ2H
            seriestype := :heatmap
            yflip := true
            yticks := (y_ticks, y_labels)
            yguide := "Depth [mm]"
            colorbar_title := "δ2H [‰]"
            clims := clims_d2H
            colorbar := true # colorbar := true_to_check_colorbar
            subplot := idx_d2H_SWAT

            # and other arguments:
            x, y_extended, z3_extended;
        end
        if (RWUcentroid == :showRWUcentroid)
            @series begin
                color := :white
                label := RWUcentroidLabel
                bg_legend --> colorant"rgba(100%,100%,100%,0.0)";legend := :right; fg_legend --> :transparent; legendfontcolor := :white
                #bg_legend --> colorant"rgba(100%,100%,100%,0.0)"; legend := :right; fg_legend --> :transparent; legendfontcolor := :white
                yflip := true; yticks := (y_ticks, y_labels)
                yguide := "Depth [mm]"; colorbar_title := "δ2H [‰]"
                subplot := idx_d2H_SWAT
                x, row_RWU_centroid_mm'
            end
        end
    end

    # 3c) Colorbar
    # ....
    # TODO: this code needs to reproduce below steps:
    # pl_colorbar_δ18O = plot([0,0], [0,1], zcolor=[0,1], t=:scatter, xlims=(1,1.1), # set xlims that no marker is shown
    #                         clims=clims_d18O, colorbar_title="δ180 [‰]",
    #                         grid=false, showaxis=false, ticks=false, label=false);
    # pl_colorbar_δ2H  = plot([0,0], [0,1], zcolor=[0,1], t=:scatter, xlims=(1,1.1), # set xlims that no marker is shown
    #                         clims=clims_d18O, colorbar_title="δ2H [‰]",
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



"""
    plotforcingandstates(simulation::DiscretizedSPAC)

Plots the forcing, states and major fluxes as results of a SPAC Simulation.
"""
@userplot PlotForcingAndStates
@recipe function f(plam::PlotForcingAndStates)
    # 0) parse input arguments
    if length(plam.args) == 1
        simulation = plam.args[1]
    else
        error("plotforcingandstates requires an unnamed first argument of type DiscretizedSPAC. Other arguments to plot() should be separated by `;`.")
    end
    if !(simulation isa DiscretizedSPAC)
        error("First unnamed argument to plotamounts should be of type DiscretizedSPAC. Got: $(typeof(simulation))")
    end
    if isnothing(simulation.ODESolution)
        error("plotamounts requires a solved system. Please `simulate!()` the DiscretizedSPAC first.")
    end

    # 1) prepare data to plot
    solu = simulation.ODESolution;

    # 1a) extract data from solution object `solu`
    # 1a) i) forcing (wind, vappres, globrad, tmax, tmin, prec, lai, )
    x1  = simulation.ODESolution_datetime;
    y11 = simulation.ODESolution.prob.p.p_WIND.(simulation.ODESolution.t);      lbl11 = "p_WIND [m/s]";
    y12 = simulation.ODESolution.prob.p.p_VAPPRES.(simulation.ODESolution.t);   lbl12 = "p_VAPPRES [kPa]";
    y13 = simulation.ODESolution.prob.p.p_GLOBRAD.(simulation.ODESolution.t);   lbl13 = "p_GLOBRAD [MJ/day/m2]"
    y14 = hcat(simulation.ODESolution.prob.p.p_TMIN.(simulation.ODESolution.t),
                    simulation.ODESolution.prob.p.p_TMAX.(simulation.ODESolution.t)); lbl14 = ["p_TMIN [°C]" "p_TMAX [°C]"]
    y15 = simulation.ODESolution.prob.p.p_PREC.(simulation.ODESolution.t);      lbl15 = "p_PREC [mm]"
    # plot_forcing = plot(layout = (:,1),
    #     plot(x1, y11; labels = lbl11),
    #     plot(x1, y12; labels = lbl12),
    #     plot(x1, y13; labels = lbl13),
    #     plot(x1, y14; labels = lbl14),
    #     plot(x1, y15; labels = lbl15))

    x2 = simulation.ODESolution_datetime;
    y21 = hcat(simulation.ODESolution.prob.p.p_DENSEF.(simulation.ODESolution.t),
                        simulation.ODESolution.prob.p.p_SAI.(simulation.ODESolution.t),
                        simulation.ODESolution.prob.p.p_LAI.(simulation.ODESolution.t)); lbl21 = ["p_DENSEF [-]" "p_SAI [-]" "p_LAI [-]"];
    y22 = simulation.ODESolution.prob.p.p_HEIGHT.(simulation.ODESolution.t); lbl22 = "p_HEIGHT [-]";
    # plot_vegetation = plot(layout = (2,1), title = "Vegetation",
    #     plot(x2, y21; labels = lbl21),
    #     plot(x2, y22;labels = lbl22))
    #     # plot!(twinx(),100*rand(10))
    #     # plot!(twinx(), simulation.ODESolution_datetime, simulation.ODESolution.prob.p.p_HEIGHT.(simulation.ODESolution.t))#;labels = "p_HEIGHT [-]")

    # 1a) ii) states (scalar and vector)
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) = LWFBrook90.get_auxiliary_variables(simulation.ODESolution)

    x3 = simulation.ODESolution_datetime;
    y31 = hcat(reduce(hcat, [simulation.ODESolution(t).INTR.mm   for t in simulation.ODESolution.t])',
                reduce(hcat, [simulation.ODESolution(t).INTS.mm   for t in simulation.ODESolution.t])',
                reduce(hcat, [simulation.ODESolution(t).SNOW.mm   for t in simulation.ODESolution.t])',
                reduce(hcat, [simulation.ODESolution(t).SNOWLQ.mm   for t in simulation.ODESolution.t])',
                reduce(hcat, [simulation.ODESolution(t).XYLEM.mm   for t in simulation.ODESolution.t])');
    lbl31 = ["INTR [mm]" "INTS [mm]" "SNOW [mm]" "SNOWLQ [mm]" "XYLEM [mm]"];
    y32 = hcat(sum(u_SWATI; dims=1)',
                reduce(hcat, [simulation.ODESolution(t).GWAT.mm   for t in simulation.ODESolution.t])');
    lbl32 = ["Total Soil Water [mm]" "GWAT [mm]"];
    # pl_states1 = plot(x3,
    #         # yaxis=:log,  1 .+
    #         y31; title = "Scalar states", labels = lbl31)
    # pl_states2 = plot(x3, y32 = ; title = "Soil/ground water", labels = lbl32)
    # plot_states = plot(pl_states1, pl_states2, layout = (2,1))

    # 1a) iii) fluxes
    x4 = simulation.ODESolution_datetime;
    y41 = hcat(  reduce(hcat, [simulation.ODESolution(t).RWU.mmday   for t in simulation.ODESolution.t])',
                                        sum(reduce(hcat, [simulation.ODESolution(t).TRANI.mmday   for t in simulation.ODESolution.t]); dims=1)')
    lbl41 = ["RWU (mm/day)" "sum(TRANI)"]
    # plot_fluxes = plot(x4, y41; labels = lbl41)

    # plot(plot_forcing, plot_vegetation, plot_states, plot_fluxes,
    #     layout = grid(4, 1, heights=[0.5 ,0.2, 0.2, 0.1]),
    #     size=(800,2000), dpi = 300, leftmargin = 15mm)

    # 2) set up the common arguments for all plots below
    # lay = RecipesBase.grid(4, 1, heights=[0.5 ,0.2, 0.2, 0.1])
    lay = RecipesBase.grid(10, 1)
    layout --> lay
    size --> (800,2000)
    dpi --> 300
    # xlim --> xlimits
    leftmargin --> 15mm

    # # 3) generate plots
    # # NOTE: --> sets attributes only when they don't already exist
    # # NOTE: :=  sets attributes even when they already exist
    @series begin
        yguide := ""
        # seriestype := :bar
        # legend := false
        subplot := 1
        label := lbl11
        # and other arguments:
        x1, y11
    end
    @series begin
        yguide := ""
        # seriestype := :bar
        # legend := false
        subplot := 2
        label := lbl12
        # and other arguments:
        x1, y12
    end
    @series begin
        yguide := ""
        # seriestype := :bar
        # legend := false
        subplot := 3
        label := lbl13
        # and other arguments:
        x1, y13
    end
    @series begin
        yguide := ""
        # seriestype := :bar
        # legend := false
        subplot := 4
        label := lbl14
        # and other arguments:
        x1, y14
    end
    @series begin
        yguide := ""
        # seriestype := :bar
        # legend := false
        subplot := 5
        label := lbl15
        # and other arguments:
        x1, y15
    end
    @series begin
        yguide := ""
        # seriestype := :bar
        # legend := false
        subplot := 6
        label := lbl21
        # and other arguments:
        x2, y21
    end
    @series begin
        yguide := ""
        # seriestype := :bar
        # legend := false
        subplot := 7
        label := lbl22
        # and other arguments:
        x2, y22
    end
    @series begin
        yguide := ""
        # seriestype := :bar
        # legend := false
        subplot := 8
        label := lbl31
        # and other arguments:
        x3, y31
    end
    @series begin
        yguide := ""
        # seriestype := :bar
        # legend := false
        subplot := 9
        label := lbl32
        # and other arguments:
        x3, y32
    end
    @series begin
        yguide := ""
        # seriestype := :bar
        # legend := false
        subplot := 10
        label := lbl41
        # and other arguments:
        x4, y41
    end
end



##########################
# Functions to get values linked to soil domain:

function find_soilDiscr_indices(simulation::DiscretizedSPAC, depths_to_read_out_mm)
    solution = simulation.ODESolution
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"

    find_soilDiscr_indices(depths_to_read_out_mm, solution)
end
# TODO(bernhard): get rid of all uses of find_soilDiscr_indices(depths_to_read_out_mm, solution::ODESolution)
function find_soilDiscr_indices(depths_to_read_out_mm, solution::ODESolution)
    # depths and lower_boundaries must all be positive numbers
    @assert all(depths_to_read_out_mm .> 0)

    lower_boundaries = cumsum(solution.prob.p.p_soil.p_THICK)

    # Only read out values that are within the simulation domain
    depths_to_read_out_mm = depths_to_read_out_mm[depths_to_read_out_mm .<= maximum(lower_boundaries)]

    idx_to_read_out = []
    for curr_depth_mm in depths_to_read_out_mm
        # idx_to_read_out = findfirst(curr_depth_mm .<= y)
        append!(idx_to_read_out, findfirst(curr_depth_mm .<= lower_boundaries))
    end
    return idx_to_read_out
end

function get_auxiliary_variables(simulation::DiscretizedSPAC; days_to_read_out_d = nothing)
    solution = simulation.ODESolution
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"

    get_auxiliary_variables(solution; days_to_read_out_d = days_to_read_out_d)
end
# TODO(bernhard): get rid of all uses of get_auxiliary_variables(solution::ODESolution; days_to_read_out_d = nothing)
function get_auxiliary_variables(solution::ODESolution; days_to_read_out_d = nothing)
    p_soil = solution.prob.p.p_soil
    NLAYER = p_soil.NLAYER
    if isnothing(days_to_read_out_d)
        u_SWATI = reduce(hcat, [solution[t_idx].SWATI.mm  for t_idx = eachindex(solution)])
    else
        u_SWATI = reduce(hcat, [solution(t_days).SWATI.mm for t_days = days_to_read_out_d])
    end

    u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK =
        (fill(NaN, size(u_SWATI)) for i in 1:5)
    for t in 1:size(u_SWATI,2)
        (u_aux_WETNES[:, t], u_aux_PSIM[:, t], u_aux_PSITI[:, t], u_aux_θ[:, t], p_fu_KK[:, t]) =
            KPT.derive_auxiliary_SOILVAR(u_SWATI[:,t], p_soil)
    end
    # returns arrays of dimenstion (t,z) where t is number of timesteps and z number of computational layers
    return (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK)
end

"""
    get_θ(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)

Returns a 2D matrix of volumetric soil moisture values (m3/m3) with soil layers as rows and time steps as columns.
The user can define timesteps as `days_to_read_out_d` or specific depths as `depths_to_read_out_mm`,
that are both optionally provided as numeric vectors, e.g. `depths_to_read_out_mm = [100, 150]` or `days_to_read_out_d = 1:1.0:100`
"""
function get_θ(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)
    solution = simulation.ODESolution
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"

    # get auxiliary variables with requested time resolution (i.e. days_to_read_out_d)
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        get_auxiliary_variables(simulation; days_to_read_out_d = days_to_read_out_d)

    # return requested soil layers
    if isnothing(depths_to_read_out_mm)
        return u_aux_θ[:, :]
    else
        return u_aux_θ[find_soilDiscr_indices(simulation, depths_to_read_out_mm), :]
    end
end
"""
    get_ψ(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)

Returns a 2D matrix of soil matric potential (kPa) with soil layers as rows and time steps as columns.
The user can define timesteps as `days_to_read_out_d` or specific depths as `depths_to_read_out_mm`,
that are both optionally provided as numeric vectors, e.g. `depths_to_read_out_mm = [100, 150]` or `saveat = 1:1.0:100`
"""
function get_ψ(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)
    solution = simulation.ODESolution
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"

    # get auxiliary variables with requested time resolution (i.e. days_to_read_out_d)
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        get_auxiliary_variables(simulation; days_to_read_out_d = days_to_read_out_d)

    # return requested soil layers
    if isnothing(depths_to_read_out_mm)
        return u_aux_PSIM[:, :]
    else
        return u_aux_PSIM[find_soilDiscr_indices(simulation, depths_to_read_out_mm), :]
    end
end
"""
    get_WETNES(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)

Returns a 2D matrix of soil wetness (-) with soil layers as rows and time steps as columns.
The user can define timesteps as `days_to_read_out_d` or specific depths as `depths_to_read_out_mm`,
that are both optionally provided as numeric vectors, e.g. `depths_to_read_out_mm = [100, 150]` or `saveat = 1:1.0:100`
"""
function get_WETNES(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)
    solution = simulation.ODESolution
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"

    # get auxiliary variables with requested time resolution (i.e. days_to_read_out_d)
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        get_auxiliary_variables(simulation; days_to_read_out_d = days_to_read_out_d)

    # return requested soil layers
    if isnothing(depths_to_read_out_mm)
        return u_aux_WETNES[:, :]
    else
        return u_aux_WETNES[find_soilDiscr_indices(simulation, depths_to_read_out_mm), :]
    end
end
"""
    get_SWATI(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)

Returns a 2D matrix of soil water volumes contained in discretized layers (mm) with soil layers as rows and time steps as columns.
Note that the values depend on the thickness of the layers and thus on the discretization.
The user can define timesteps as `days_to_read_out_d` or specific depths as `depths_to_read_out_mm`,
that are both optionally provided as numeric vectors, e.g. `depths_to_read_out_mm = [100, 150]` or `saveat = 1:1.0:100`
"""
function get_SWATI(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)
    solution = simulation.ODESolution
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"

    # get auxiliary variables with requested time resolution (i.e. days_to_read_out_d)
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        get_auxiliary_variables(simulation; days_to_read_out_d = days_to_read_out_d)

    # return requested soil layers
    if isnothing(depths_to_read_out_mm)
        return u_SWATI[:, :]
    else
        return u_SWATI[find_soilDiscr_indices(simulation, depths_to_read_out_mm), :]
    end
end
"""
    get_K(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)

Returns a 2D matrix of soil hydraulic conductivities (mm/day) with soil layers as rows and time steps as columns.
The user can define timesteps as `days_to_read_out_d` or specific depths as `depths_to_read_out_mm`,
that are both optionally provided as numeric vectors, e.g. `depths_to_read_out_mm = [100, 150]` or `saveat = 1:1.0:100`
"""
function get_K(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)
    solution = simulation.ODESolution
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"

    # get auxiliary variables with requested time resolution (i.e. days_to_read_out_d)
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        get_auxiliary_variables(simulation; days_to_read_out_d = days_to_read_out_d)

    # return requested soil layers
    if isnothing(depths_to_read_out_mm)
        return p_fu_KK[:, :]
    else
        return p_fu_KK[find_soilDiscr_indices(simulation, depths_to_read_out_mm), :]
    end
end

"""
    get_δsoil(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)

Returns tuple of two 2D matrices of isotopic signatures of soil water (δ in permil) for d18O and d2H.
The 2D matrix with soil layers as rows and time steps as columns can be accessed with `.d18O` and `.d2H`, respectively.
The user can define timesteps as `days_to_read_out_d` or specific depths as `depths_to_read_out_mm`,
that are both optionally provided as numeric vectors, e.g. `depths_to_read_out_mm = [100, 150]` or `saveat = 1:1.0:100`
"""
function get_δsoil(simulation::DiscretizedSPAC; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing)
    solution = simulation.ODESolution
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"
    @assert solution.prob.p.simulate_isotopes "Provided solution did not simulate isotopes"

    # get auxiliary variables with requested time resolution (i.e. days_to_read_out_d)
    if isnothing(days_to_read_out_d)
        # vector quantities δ-values:
        rows_SWAT_d18O = reduce(hcat, [solution[t_idx].SWATI.d18O for t_idx = eachindex(solution)])
        rows_SWAT_d2H  = reduce(hcat, [solution[t_idx].SWATI.d2H  for t_idx = eachindex(solution)])
    else
        # vector quantities δ-values:
        rows_SWAT_d18O = reduce(hcat, [solution(t_days).SWATI.d18O for t_days = days_to_read_out_d])
        rows_SWAT_d2H  = reduce(hcat, [solution(t_days).SWATI.d2H  for t_days = days_to_read_out_d])
    end
    # return requested soil layers
    if isnothing(depths_to_read_out_mm)
        return (d18O = rows_SWAT_d18O[:, :],
                d2H  = rows_SWAT_d2H[ :, :])
    else
        idx_soil_layers = find_soilDiscr_indices(simulation, depths_to_read_out_mm)
        return (d18O = rows_SWAT_d18O[idx_soil_layers, :],
                d2H  = rows_SWAT_d2H[ idx_soil_layers, :])
    end
end

##########################
# Functions to get values linked to aboveground:
"""
    get_aboveground(simulation::DiscretizedSPAC; days_to_read_out_d = nothing)

Returns a DataFrame with state variables: GWAT, INTS, INTR, SNOW, SNOWLQ in mm and CC in MJm2.
By default, the values are returned for each simulation timestep.
The user can define timesteps as `days_to_read_out_d` by optionally providing a numeric vector, e.g. saveat = 1:1.0:100
"""
function get_aboveground(simulation::DiscretizedSPAC; days_to_read_out_d = nothing)
    solution = simulation.ODESolution
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"

    # compartments_to_extract = [:GWAT, :INTS, :INTR, :SNOW, :SNOWLQ]
    compartments_to_extract = [:GWAT, :INTS, :INTR, :SNOW, :CC, :SNOWLQ] # CC is in MJm2

    if isnothing(days_to_read_out_d)
        # u_aboveground = [solution[t_idx][compartment].mm  for t_idx = eachindex(solution), compartment in compartments_to_extract]
        u_aboveground = [solution[t_idx][compartment][1]  for t_idx = eachindex(solution), compartment in compartments_to_extract]  # CC is in MJm2, by using [1] we avoid specifiyng :mm or :MJm2
        # [solution[t_idx][compartment].d18O  for t_idx = eachindex(solution), compartment in compartments_to_extract]
        # [solution[t_idx][compartment].d2H  for t_idx = eachindex(solution), compartment in compartments_to_extract]
    else
        # u_aboveground = [solution(t_days)[compartment].mm  for t_days = days_to_read_out_d, compartment in compartments_to_extract]
        u_aboveground = [solution(t_days)[compartment][1]  for t_days = days_to_read_out_d, compartment in compartments_to_extract] # CC is in MJm2, by using [1] we avoid specifiyng :mm or :MJm2
        # [solution[t_idx][compartment].d18O  for t_idx = eachindex(solution), compartment in compartments_to_extract]
        # [solution[t_idx][compartment].d2H  for t_idx = eachindex(solution), compartment in compartments_to_extract]
    end

    # returns DataFrame of dimenstion (t,6) where t is number of time steps and 6 number of aboveground compartments
    return DataFrame(u_aboveground, string.(compartments_to_extract))
end

"""
    get_δ(simulation::DiscretizedSPAC; days_to_read_out_d = nothing)

Returns a DataFrame with the isotopoic compositions of the inputs and state variables: PREC, GWAT, INTS, INTR, SNOW, RWU, XYLEM.
By default, the values are returned for each simulation timestep.
The user can define timesteps as `days_to_read_out_d` by optionally providing a numeric vector, e.g. saveat = 1:1.0:100
"""
function get_δ(simulation::DiscretizedSPAC; days_to_read_out_d = nothing)
    # function get_δ(simulation::DiscretizedSPAC; days_to_read_out_d = nothing, compartment::Symbol = :all)
    solution = simulation.ODESolution
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"
    @assert solution.prob.p.simulate_isotopes "Provided solution did not simulate isotopes"

    compartments_to_extract = [:GWAT, :INTS, :INTR, :SNOW, :RWU, :XYLEM]

    if isnothing(days_to_read_out_d)
        # row_PREC_d18O = reshape(solution.prob.p.p_δ18O_PREC.(solution.t), 1, :)
        # row_PREC_d2H  = reshape(solution.prob.p.p_δ2H_PREC.(solution.t), 1, :)
        PREC_d18O = solution.prob.p.p_δ18O_PREC.(solution.t)
        PREC_d2H  = solution.prob.p.p_δ2H_PREC.( solution.t)

        # u_aboveground = [solution[t_idx][compartment].mm  for t_idx = eachindex(solution), compartment in compartments_to_extract]
        d18O_aboveground = [solution[t_idx][compartment].d18O  for t_idx = eachindex(solution), compartment in compartments_to_extract]
        d2H_aboveground  = [solution[t_idx][compartment].d2H   for t_idx = eachindex(solution), compartment in compartments_to_extract]
    else
        # row_PREC_d18O = reshape(solution.prob.p.p_δ18O_PREC.(days_to_read_out_d), 1, :)
        # row_PREC_d2H  = reshape(solution.prob.p.p_δ2H_PREC.(days_to_read_out_d), 1, :)
        PREC_d18O = solution.prob.p.p_δ18O_PREC.(days_to_read_out_d)
        PREC_d2H  = solution.prob.p.p_δ2H_PREC.( days_to_read_out_d)

        # u_aboveground = [solution(t_days)[compartment].mm  for t_days = days_to_read_out_d, compartment in compartments_to_extract]
        d18O_aboveground = [solution(t_days)[compartment].d18O  for t_days = days_to_read_out_d, compartment in compartments_to_extract]
        d2H_aboveground  = [solution(t_days)[compartment].d2H   for t_days = days_to_read_out_d, compartment in compartments_to_extract]
    end

    return DataFrame([PREC_d18O d18O_aboveground PREC_d2H d2H_aboveground],
                     ["PREC_d18O"; string.(compartments_to_extract).*"_d18O";
                      "PREC_d2H" ; string.(compartments_to_extract).*"_d2H"])

end

# also make non-unicode variants:
get_theta     = get_θ
get_psi       = get_ψ
get_deltasoil = get_δsoil

get_delta     = get_δ

############################################################################################
############################################################################################
############################################################################################