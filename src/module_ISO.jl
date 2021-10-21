@doc raw"""
# Isotopic mixing and fractionation

This module defines constants and functions associated with mixing and isotopic fractionation
processes simulated by LWFBrook90.jl.

States INTS, INTR, SNOW are updated once per day. The isotopic composition at the end of the
day is computed based on an operator splitting approach. First, non-fractionating in- and
outflows are considered for each compartment and a resulting isotopic composition is computed.
Second, the fractionation associated with the evaporation process is computed, assuming a
Craig-Gordon fractionation process applied to a reservoir with the evaporation process as
single sink. The reservoirs inital condition is the isotopic composition resulting from the
mixing computed in the first operator step.

Turbulence conditions are fixed to different values for each compartment (X_INTR, X_INTS, X_SNOW).

Evaporation from the soil is modelled from the uppermost layer only.TODO(bernhard): Should this be improved?
Turbulence conditions for evaporation from the soil layer are conditional on the soil moisture
status $θ$ following (Zhou-2021-Environ_Model_Softw (X is called n_k there)) as a weighted
average between $X_s = 1$ (molecular diffusion only) and $X_s = 0.5$ (both molecular and turbulent diffusion):
```math
X_{soil} = \frac{(θ-θ_{res})*X_a + (θ_{sat}-θ)*X_s}{(θ_{sat}-θ_{res})}
```

Enrichment of xylem water due to evaporation from the stomatas is not modelled.
"""
module ISO # ISOTOPIC MIXING AND FRACTIONATION

using RecipesBase
using ..LWFBrook90: RelativeDaysFloat2DateTime

export α¹⁸O_dif, α²H_dif, α¹⁸O_eq, α²H_eq, δ_CraigGordon

# 1a) Kinetic fractionation (Gonfiantini 2018):
α¹⁸O_dif = 1.0285 # -, D/D_i i.e. ¹H¹H¹⁶O/¹H¹H¹⁸O: Taken from Gonfiantini 2018, citing Merlivat 1978
α²H_dif  = 1.0251 # -, D/D_i i.e. ¹H¹H¹⁶O/¹H¹H¹⁸O: Taken from Gonfiantini 2018, citing Merlivat 1978

# 1b) Equilibrium fractionation
function α¹⁸O_eq(temp_celsius = 25.) # temp_celsius is temperature in Celsius
    T_Kel = temp_celsius + 273.15
    # Majoube
    A,B,C = 1137,-0.4156,-0.00207 # taken from Gonfiantini 1986/2018 citing Majoube
    return(exp(A*T_Kel^(-2)+B*T_Kel^(-1)+C))
    # # Horita&Weselowski
    # T_mK = T_Kel/10^3 # divided by 1000 [mK]
    # K0,Km1,Km2,Km3 = -7.685,6.7123,-1.6664,0.3504 # taken from Benettin 2018 citing Horita&Weselowski
    # return(exp(10^-3 * (K0 + Km1*T_mK^-1 + Km2*T_mK^-2 + Km3*T_mK^-3)))
end
function α²H_eq(temp_celsius = 25.) # temp_celsius is temperature in Celsius
    T_Kel = temp_celsius + 273.15
    # Majoube
    A,B,C = 24844,-76.248,0.05261 # taken from Gonfiantini 1986/2018 citing Majoube
    return(exp(A*T_Kel^(-2)+B*T_Kel^(-1)+C))
    # # Horita&Weselowski
    # T_mK = T_Kel/10^3 # divided by 1000 [mK]
    # K3,K2,K1,K0,Km3 = 1158.8,-1620.1,794.84,-161.04,2.9992 # taken from Benettin 2018 citing Horita&Weselowski
    # return(exp(10^-3 * (K3*T_mK^3 + K2*T_mK^2 + K1*T_mK + K0 + Km3*T_mK^-3)))
end

# 7a) source: equations from Gonfiantini 2018
function δ_CraigGordon(δ0, δΑ, f, h, α_eq, α_dif, γ, X)
    # A(h, γ, α_dif, X)       = -h / (α_dif^X * (γ-h))            # Eq 12
    # B(h, γ, α_eq, α_dif, X) = γ / (α_eq * α_dif^X * (γ-h)) - 1  # Eq 13
    # Rw(A,B,RA,Rw0,f)        = -A/B*RA + (Rw0 + A/B*RA)*f^B      # Eq 14

    A = -h / (α_dif^X * (γ-h))            # Eq 12
    B = γ / (α_eq * α_dif^X * (γ-h)) - 1  # Eq 13

    return (-1 -A/B*(δΑ + 1) + (δ0 + 1 + A/B*(δΑ+1)) * f^B) # (-)
end


@userplot PlotIsotopes

@recipe function f(h::PlotIsotopes)
    # 0) parse input arguments
    # if length(h.args) != 1 # || !(typeof(h.args[1]) <: SciMLBase.ODESolution)
    #     error("Plot Isotopes should be given two arguments of type ODESolution.  Got: $(typeof(h.args))")
    # end
    if length(h.args) == 1 # || !(typeof(h.args[1]) <: SciMLBase.ODESolution)
        # all okay, use default tick_function
        sol = h.args[1]
        tick_function = (x1,x2) -> [range(x1, x2; length=10)][1]
    elseif length(h.args) == 2
        # all okay, tick_function was provided
        # check that second argument is indeed a function
        if isempty(methods(h.args[2]))
            error("Plot Isotopes' second argument should be a function for the ticks.")
        end
        sol = h.args[1]
        tick_function = h.args[2]
    else
        error("Plot Isotopes should be given two arguments of type ODESolution.  Got: $(typeof(h.args))")
    end



    # 1) data to plot

    # 1a) extract data from solution object `sol`
    # idx_u_vector_amounts       = sol.prob.p[1][4][4]
    idx_u_scalar_isotopes_d18O = sol.prob.p[1][4][5]
    idx_u_vector_isotopes_d18O = sol.prob.p[1][4][6]
    idx_u_scalar_isotopes_d2H  = sol.prob.p[1][4][7]
    idx_u_vector_isotopes_d2H  = sol.prob.p[1][4][8]
    # idx_u_vector_accumulators  = sol.prob.p[1][4][9]

    t_ref = sol.prob.p[2][17]
    x = RelativeDaysFloat2DateTime.(sol.t, t_ref);
    y = cumsum(sol.prob.p[1][1].p_THICK) - sol.prob.p[1][1].p_THICK/2

    row_PREC_amt = sol.prob.p[2][8].(sol.t)
    rows_SWAT_amt = sol[7 .+ (0:sol.prob.p[1][1].NLAYER-1),
                                1,
                                :]./sol.prob.p[1][1].p_THICK;
    rows_SWAT_d18O = sol[idx_u_vector_isotopes_d18O,1,:];
    rows_SWAT_d2H  = sol[idx_u_vector_isotopes_d2H, 1,:];

    row_NaN       = fill(NaN, 1,length(x))
    row_PREC_d18O = reshape(sol.prob.p[2][15].(sol.t), 1, :)
    row_INTS_d18O = reshape(sol[idx_u_scalar_isotopes_d18O[2],1,:], 1, :)
    row_INTR_d18O = reshape(sol[idx_u_scalar_isotopes_d18O[3],1,:], 1, :)
    row_SNOW_d18O = reshape(sol[idx_u_scalar_isotopes_d18O[4],1,:], 1, :)
    row_GWAT_d18O = reshape(sol[idx_u_scalar_isotopes_d18O[1],1,:], 1, :)
    row_PREC_d2H  = reshape(sol.prob.p[2][16].(sol.t), 1, :)
    row_INTS_d2H  = reshape(sol[idx_u_scalar_isotopes_d2H[2],1,:], 1, :)
    row_INTR_d2H  = reshape(sol[idx_u_scalar_isotopes_d2H[3],1,:], 1, :)
    row_SNOW_d2H  = reshape(sol[idx_u_scalar_isotopes_d2H[4],1,:], 1, :)
    row_GWAT_d2H  = reshape(sol[idx_u_scalar_isotopes_d2H[1],1,:], 1, :)

    # 1b) define some plot arguments based on the extracted data
    # color scheme:
    all_d18O_values = [rows_SWAT_d18O; row_PREC_d18O; row_INTS_d18O; row_INTR_d18O; row_SNOW_d18O; row_GWAT_d18O][:]
    all_d2H_values  = [rows_SWAT_d2H ; row_PREC_d2H ; row_INTS_d2H ; row_INTR_d2H ; row_SNOW_d2H ; row_GWAT_d2H ][:]
    clims_d18O = (-16, -6)
    clims_d2H  = (-125, -40)
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

    true_to_check_colorbar = true; # set this flag to false for final plot, true for debugging.



    # 2) set up the common arguments for all plots below
    # link := :x #TODO(bernhard): link doesn't seem to work...
    # layout --> (4,1)
    # layout := @layout [precip        , -
    #                    isotop_heatmap, colorbarlegend]
    # @layout unsupported: https://github.com/JuliaPlots/RecipesBase.jl/issues/15
    # TODO(bernhard): find an easy way to do: l = @layout [grid(2, 1, heights=[0.2, 0.8]) a{0.055w}]
    # Possibly it needs to be provided when calling `plot_LWFBrook90_isotopes(... ; layout = ...)`
    # Then make sure to number the subplots correctly and include the colorbars from 3c)


    # 3) generate plots
    # NOTE: --> sets attributes only when they don't already exist
    # NOTE: :=  sets attributes even when they already exist

    # 3a) Precipitation
    # We reproduce the following plot:
    # Workaround for barplot # https://github.com/JuliaPlots/Plots.jl/issues/3880
    # ts_PREC_δ18O = plot(reshape(x,1,:), reshape(row_PREC_amt,1,:), t = :bar, line_z = reshape(row_PREC_d18O,1,:), fill_z = reshape(row_PREC_d18O,1,:),
    #                 clims=clims_d18O, colorbar = true_to_check_colorbar, legend = false,
    #                 ylabel = "PREC [mm]"); # TODO(bernhard): make this work with a barplot
    # ts_PREC_δ2H = plot(reshape(x,1,:), reshape(row_PREC_amt,1,:), t = :bar, line_z = reshape(row_PREC_d2H,1,:), fill_z = reshape(row_PREC_d2H,1,:),
    #                 clims=clims_d2H,  colorbar = true_to_check_colorbar, legend = false,
    #                 ylabel = "PREC [mm]"); # TODO(bernhard): make this work with a barplot
    @series begin # ts_PREC_δ18O
        title := "δ18O"
        seriestype := :bar
        # linecolor := :match
        line_z := reshape(row_PREC_d18O,1,:)
        fill_z := reshape(row_PREC_d18O,1,:)
        clims := clims_d18O
        colorbar_title := "δ18O [mUr]"
        colorbar := true_to_check_colorbar # TODO: define this once for all plots (except force it for colorbar plot representing common legend)
        yguide := "PREC [mm]"
        legend := false
        subplot := 1

        # and other arguments:
        reshape(x,1,:), reshape(row_PREC_amt,1,:)
    end
    @series begin # ts_PREC_δ2H
        title := "δ2H"
        seriestype := :bar
        # linecolor := :match
        line_z := reshape(row_PREC_d2H,1,:)
        fill_z := reshape(row_PREC_d2H,1,:)
        clims := clims_d2H
        colorbar_title := "δ2H [mUr]"
        colorbar := true_to_check_colorbar # TODO: define this once for all plots (except force it for colorbar plot representing common legend)
        yguide := "PREC [mm]"
        legend := false
        subplot := 3

        # and other arguments:
        reshape(x,1,:), reshape(row_PREC_amt,1,:)
    end

    # 3b) Heatmap (containing SWATI and other compartments)
    y_extended = [-500; -350; -300; -250; -200; -150; -100;   -50;         y;          (maximum(y) .+ [10, 200])]
    # y_labels   = ["INTS"; ""; "INTR"; ""; "SNOW"; ""; round.(y); "";             "GWAT"]
    # y_soil_ticks = optimize_ticks(extrema(y)...; k_min = 4)[1]
    # y_soil_ticks = optimize_ticks(0., round(maximum(cumsum(sol.prob.p[1][1].p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
    y_soil_ticks = tick_function(0., round(maximum(cumsum(sol.prob.p[1][1].p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
    y_ticks    = [-500;       -300;       -200;       -100;     y_soil_ticks;                   (maximum(y) .+ [   200])]
    y_labels   = ["PREC";   "INTS";     "INTR";     "SNOW";     round.(y_soil_ticks; digits=0); "GWAT"]
    z2_extended = [row_PREC_d18O; row_NaN; row_INTS_d18O; row_NaN; row_INTR_d18O; row_NaN; row_SNOW_d18O; row_NaN; rows_SWAT_d18O; row_NaN; row_GWAT_d18O]
    z3_extended = [row_PREC_d2H;  row_NaN; row_INTS_d2H;  row_NaN; row_INTR_d2H;  row_NaN; row_SNOW_d2H;  row_NaN; rows_SWAT_d2H;  row_NaN; row_GWAT_d2H]

    # We reproduce the following plots:
    # pl_δ18O = heatmap(x, y_extended, z2_extended;
    #     yflip = true,
    #     yticks = (y_ticks, y_labels), colorbar = true_to_check_colorbar,
    #     ylabel = "Depth [mm]",
    #     colorbar_title = "δ18O [mUr]");
    # pl_δ2H = heatmap(x, y_extended, z3_extended;
    #     yflip = true,
    #     yticks = (y_ticks, y_labels), colorbar = true_to_check_colorbar,
    #     ylabel = "Depth [mm]",
    #     colorbar_title = "δ2H [mUr]");

    @series begin # pl_δ18O
        seriestype := :heatmap
        yflip := true
        yticks := (y_ticks, y_labels)
        colorbar := true_to_check_colorbar
        yguide := "Depth [mm]"
        colorbar_title := "δ18O [mUr]"
        clims := clims_d18O
        colorbar := true_to_check_colorbar # TODO: define this once for all plots (except force it for colorbar plot representing common legend)
        subplot := 2

        # and other arguments:
        # x, y_extended, z2_extended
        x, y_extended, z2_extended;
    end

    @series begin # pl_δ2H
        seriestype := :heatmap
        yflip := true
        yticks := (y_ticks, y_labels)
        colorbar := true_to_check_colorbar
        yguide := "Depth [mm]"
        colorbar_title := "δ2H [mUr]"
        clims := clims_d2H
        colorbar := true_to_check_colorbar # TODO: define this once for all plots (except force it for colorbar plot representing common legend)
        subplot := 4

        # and other arguments:
        # x, y_extended, z2_extended
        x, y_extended, z3_extended;
    end

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

end
