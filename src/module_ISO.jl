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

using RecipesBase, PlotUtils, Measures
using ..LWFBrook90: RelativeDaysFloat2DateTime, DiscretizedSPAC

export α¹⁸O_dif, α²H_dif, α¹⁸O_eq, α²H_eq, δ_CraigGordon, update_δ_with_mixing_and_evaporation
export R_VSMOW¹⁸O, R_VSMOW²H, Mi_¹⁸O, Mi_²H, Mw
export δ_to_x, x_to_δ, dxdt_to_dδdt, δ_to_C, C_to_δ

# Isotopic ratios of standard ocean water VSMOW (reference for definition of δ)
R_VSMOW¹⁸O = 2005.2e-6 # (source: Baertschi-1976-Earth_Planet_Sci_Lett)
R_VSMOW²H  = 155.76e-6 # (source: Hagemann-1970-Tellus)
Mi_¹⁸O = 0.020 # Molar mass of ¹H¹H¹⁸O in kg
Mi_²H  = 0.019 # Molar mass of ¹H²H¹⁶O in kg
Mw     = 0.018 # Molar mass of ¹H¹H¹⁶O in kg

#TODO(bernhard): debug issues and switch this back on...
# δ to C (and back) implementation below is approximative (assuming Ni*Mi << Nw*Mw)
δ_to_C(δ_permil,R_std, Mi; ρw_kg_m3 = 1000) = R_std .* (δ_permil./1000 .+ 1) .* ρw_kg_m3 .* Mi./Mw # returns C in kg/m3, source Eq.5 Braud et al. 2005
C_to_δ(C,R_std, Mi; ρw_kg_m3 = 1000)        = 1000 .* (C ./ (R_std .* ρw_kg_m3 .* Mi./Mw) .- 1 )   # returns δ in permil, source Eq.5 Braud et al. 2005
dCdt_to_dδdt(dCdt, C, R_std, Mi; ρw_kg_m3 = 1000) = dCdt .* 1 ./ (R_std .* ρw_kg_m3 .* Mi./Mw) .* 1000
# δ_to_C(δ_permil,R_std, Mi; ρw_kg_m3 = 1000) = δ_permil # R_std .* (δ_permil./1000 .+ 1) .* ρw_kg_m3 .* Mi./Mw # returns C in kg/m3, source Eq.5 Braud et al. 2005
# C_to_δ(C,R_std, Mi; ρw_kg_m3 = 1000)        = C #1000 .* (C ./ (R_std .* ρw_kg_m3 .* Mi./Mw) .- 1 )   # returns δ in permil, source Eq.5 Braud et al. 2005
# δ_to_C(δ_permil,R_std, Mi; ρw_kg_m3 = 1000) = 1 ./ (1 .+ 1 ./ (R_std .* ( δ_permil./1000 .+ 1 )) )
# C_to_δ(x,R_std, Mi; ρw_kg_m3 = 1000)        = (1 ./ ( R_std .* (1 ./ x .- 1) ) .- 1) .* 1000 # in permil

# δ to x (and back) implementation below is exact
δ_to_x(δ_permil,R_std)       = 1 ./ (1 .+ 1 ./ (R_std .* ( δ_permil./1000 .+ 1 )) )
x_to_δ(x,R_std)              = (1 ./ ( R_std .* (1 ./ x .- 1) ) .- 1) .* 1000 # in permil
dxdt_to_dδdt(dxdt, x, R_std) = dxdt .* 1 ./ R_std .* 1 ./ (x .- 1).^2 .* 1000  # using dδ/dt = dδ/dx * dx/dt



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

function update_δ_with_mixing_and_evaporation(dt, u₀, δ₀, inflow, δin, outflow, R_std, E, δₐ, h, α_eq, α_dif, γ, X; do_fractionation=false)
    # dt       [day]      , time step
    # u₀       [mm]       , initial amount
    # δ₀       [‰]      , initial isotopic signature
    # inflow   [mm/day]   , tuple/vector of inflows
    # δin      [‰]      , tuple/vector of inflow isotopic signatures
    # outflow  [mm/day]   , outflow rates
    # R_std     [-]       , standard ratio between number of heavy and number of light molecules
    # E        [mm/day]   , evaporation rate
    # δₐ, h, γ, α_eq, α_dif, X (see below...) other evaporation specific parameters

    # inflow and δin as well as outflow can be vectors [] or tuples ()
    # NOTE: sum((1.0, 2.0).*(1.0, 2.0)) is the same as sum([1.0, 2.0].*[1.0, 2.0])

    # Operator splitting: first compute all in- and outflows and corresponding mixing to yield δ⁺
    #                     second compute the fractionation from δ⁺ to δ⁺⁺ considering ...
    #                     ... a reservoir with a single sink: δ₁ + 1 = (1 + δ₀)*(N₁/N₀)^ε = (1 + δ₀)*f^ε (with f fraction remaining)

    # Going from R₀,W₀ to (1) R⁺,W⁺ and then to (2) R⁺⁺,W⁺⁺

    # ⁺ First step: compute only in- and outflows and neglect fractionation due to evaportion
    # Mass balance: dW/dt = ∑ Win - ∑ Wout,                   and neglect in Wout the evaporation flux E
    # Isot.balance: d(R*W)/dt = ∑ Win * Rin - ∑ Wout * R(t),  and neglect in Wout the evaporation flux E
    #               when simpliying Rout: R(t) := R(t=0) as constant:
    #               R⁺*W⁺ = R₀*W₀ + ∑ Win*Rin - ∑ Wout*R₀
    #
    #        (1) ==> W⁺ = W₀ + dt * (∑ Win - ∑ Wout)
    #        (2) ==> R⁺ = R₀ * W₀/W⁺ + 1/W⁺ (∑ Rin * dt*Win - ∑ R₀ * dt*Wout)
    #
    # ⁺⁺ Second step: going from R⁺,W⁺ to R⁺⁺,W⁺⁺
    #    Using reservoir with one evaporating sink (Gonfiantini 2018)
    #       (3) ==> W⁺⁺ = W⁺ - dt*E
    #       (4) ==> R⁺⁺ = (R⁺ + Rₐ * A/B)*f^B - A/B*Rₐ, (Gonfiantini 2018: Craig-Gordon model)
    #           where A/B = h * α_eq/(γ - α_eq * α_dif^X *(γ-h))
    #           where B   = γ/(α_eq * α_dif^X *(γ-h)) - 1
    #           where f   = W⁺⁺ / W⁺ , fraction of remaining liquid
    #           where Rₐ is atmospheric isotopic composition
    #           where γ is thermodynamic activity coefficient of evaporating water
    #           where h is relative humidty of environmental atmosphere
    #           where α_eq is isotopic fractionation factor at equilibrium between liquid water and vapor, always > 1
    #           where α_dif is isotopic fractionation factor between vapor in the saturated equilibrium layer and vapor escaping by diffusion, always >= 1
    #           where X is teh turbulence index of the atmosphere above the evaporating water (0=fully turbulent, 1=still)

    # (2) and (4) can be approximated in δ notation instead of R notation as [using R = Rstd * (δ+1)]:
    # (2') δ⁺  = -1 + (δ₀ + 1) * W₀/W⁺ + dt/W⁺ (∑ (δin + 1) * Win - ∑ (δ₀ + 1) * Wout)
    # (4') δ⁺⁺ = [δ⁺ + 1 + A/B(δₐ + 1)]*f^B - [1 + A/B*(δₐ + 1)]


    # Step 0) make the balance exact by using atom fractions x
    # x = R/(1+R) = ((δ+1)*R_std) / (1 + ((δ+1)*R_std))
    x₀ = (δ₀/1000 + 1)*R_std / (1 + (δ₀/1000 + 1)*R_std)
    xin = (δin./1000 .+ 1).*R_std ./ (1 .+ (δin./1000 .+ 1).*R_std)

    #### TODO 16.03: reactivate @assert all(inflow .>= 0) "Inflows should not be negative"
    #########
    # Step 1)
    u⁺ = u₀ + dt * (sum(inflow) - sum(outflow)) # [mm]
    x⁺ = (x₀ * u₀ + dt * (sum(xin .* inflow) - sum(x₀ .* outflow)) ) / u⁺
    # old, unused: δ⁺ = -1 + (δ₀ / 1000. + 1) * u₀/u⁺ + dt/u⁺ * (sum( (δin ./ 1000. .+ 1).* inflow) - sum( (δ₀ / 1000. + 1) * outflow))

    # Edge cases:
    # a) u₀ = 0, then we had δ₀ = NA, and we want δ⁺ = δin (weighted by inflows...)
    if u⁺ == 0
        # δ⁺ = NaN
        x⁺ = NaN
    elseif u₀ == 0
        #### TODO 16.03: reactivate @assert sum(outflow) <= sum(inflow)
        x⁺ = sum(xin .* inflow) / sum(inflow)
    end

    #########
    # Step 2)
    #### TODO 16.03: reactivate @assert dt*E >= 0 "Evaporation rate must be positive"
    u⁺⁺ = u⁺ - dt*E
    # u⁺⁺ = max(0, u⁺ - dt*E) # TODO(bernhard): as below assert was triggered, this is a quick workaround
    #### TODO 16.03: reactivate @assert u⁺⁺ >= 0 "End amount must still be positive after evaporation"

    if do_fractionation
        #### TODO(bernhard): reactivation fractionation f     = u⁺⁺ / u⁺
        #### TODO(bernhard): reactivation fractionation AdivB = h * α_eq/(γ - α_eq * α_dif^X *(γ-h))
        #### TODO(bernhard): reactivation fractionation B     = γ/(α_eq * α_dif^X *(γ-h)) - 1
        #### TODO(bernhard): reactivation fractionation δ⁺⁺   = (δ⁺ + 1 + AdivB * (δₐ / 1000. + 1)) * f^B - [1 + AdivB*(δₐ / 1000. + 1)]

        #### TODO(bernhard): reactivation fractionation # Edge cases:
        #### TODO(bernhard): reactivation fractionation # b) u⁺ = 0, then we had δ⁺ = NA and thus no evaporation from u⁺ to u⁺⁺ --> δ⁺⁺ = NA
        #### TODO(bernhard): reactivation fractionation # c) u⁺⁺ = 0, then we have δ⁺⁺ = NA
        #### TODO(bernhard): reactivation fractionation if u⁺ == 0 | u⁺⁺ == 0
        #### TODO(bernhard): reactivation fractionation     u⁺⁺ = 0
        #### TODO(bernhard): reactivation fractionation     δ⁺⁺ = NaN
        #### TODO(bernhard): reactivation fractionation else
        #### TODO(bernhard): reactivation fractionation     @assert isnan(δ⁺⁺) # TODO(bernhard) include bang: @assert !isnan(δ⁺⁺)
        #### TODO(bernhard): reactivation fractionation end

        #### TODO(bernhard): reactivation fractionation # Step 3)
        #### TODO(bernhard): reactivation fractionation # consider some edge cases:
        #### TODO(bernhard): reactivation fractionation # 3a) u₀ = 0, then we had δ₀ = NA
        #### TODO(bernhard): reactivation fractionation # 3b) u⁺ = 0, then we had δ⁺ = NA and thus no evaporation from u⁺ to u⁺⁺ --> δ⁺⁺ = NA
        #### TODO(bernhard): reactivation fractionation # 3c) u⁺⁺ = 0, then we have δ⁺⁺ = NA
        #### TODO(bernhard): reactivation fractionation # 3d) else we have the formulas (3) and (4')
        #### TODO(bernhard): reactivation fractionation return (u⁺⁺, δ⁺⁺)
    else
        x⁺⁺ = x⁺
    end

    # TODO(bernhard): for debugging: don't use fractionation effect!

    # Go back to using δ in ‰
    # δ = R/R_std - 1 = x/(1-x)/R_std - 1
    # δ⁺ = 1000 * (x⁺/(1-x⁺)/R_std - 1)
    δ⁺⁺ = 1000 * (x⁺⁺/(1-x⁺⁺)/R_std - 1)

    return (u⁺⁺, δ⁺⁺)
end

# defines mutable struct `PlotIsotopes` and sets shorthands `plotisotopes` and `plotisotopes!`
@userplot PlotIsotopes
@recipe function f(pliso::PlotIsotopes)

    # 0) parse input arguments
    # @show typeof.(pliso.args)
    # @show length(pliso.args) != 1

    if length(pliso.args) != 1
        error("plotisotopes requires a single argument of type DiscretizedSPAC. Other arguments to plot() should be separated by `;`,")
    end

    simulation = pliso.args[1]
    if !(simulation isa DiscretizedSPAC)
        error("plotisotopes should be given a first argument of type DiscretizedSPAC. Got: $(typeof(simulation))")
    end
    if isnothing(simulation.ODESolution)
        error("plotisotopes requires a solved system. Please `simulate!()` the DiscretizedSPAC first.")
    end
    sol = simulation.ODESolution


    # 1) prepare data to plot

    # 1a) extract data from solution object `sol`
        # # Two ways to extract data from soil object: using `[]` or `()`
        # u_SWATI = reduce(hcat, [sol[t_idx].SWATI.mm  for t_idx = eachindex(sol)])
        # u_SWATI = reduce(hcat, [sol(t_days).SWATI.mm for t_days = saveat])

        # saveat decides which points to use:
        # saveat = nothing # read out all simulation steps
        saveat = :daily  # read out only daily values
        if isnothing(saveat)
            saveat = sol.t # warning this can
            error("Too many output times requested. This easily freezes the program...")
        elseif saveat == :daily
            saveat = unique(round.(sol.t))
        end
        t_ref = sol.prob.p.REFERENCE_DATE
        x = RelativeDaysFloat2DateTime.(saveat, t_ref);
        y_center = cumsum(sol.prob.p.p_soil.p_THICK) - sol.prob.p.p_soil.p_THICK/2

    simulate_isotopes = sol.prob.p.simulate_isotopes
    @assert simulate_isotopes "Provided DiscretizedSPAC() did not simulate isotopes"

    # Some hardcoded options:
    xlimits = RelativeDaysFloat2DateTime.(sol.prob.tspan, t_ref)

    clims_d18O = (-16, -6)
    clims_d2H  = (-125, -40)
    true_to_check_colorbar = true; # set this flag to false for final plot, true for debugging.
    tick_function = (x1, x2) -> PlotUtils.optimize_ticks(x1, x2; k_min = 4)

    # Results
    row_PREC_amt = sol.prob.p.p_PREC.(saveat)
    rows_SWAT_amt  = reduce(hcat, [sol[t].SWATI.mm   for t in eachindex(sol)]) ./ sol.prob.p.p_soil.p_THICK
    rows_SWAT_d18O = reduce(hcat, [sol[t].SWATI.d18O for t in eachindex(sol)])
    rows_SWAT_d2H  = reduce(hcat, [sol[t].SWATI.d2H  for t in eachindex(sol)])

    # rows_SWAT_amt = sol[sol.prob.p.row_idx_SWATI, 1, :]./sol.prob.p.p_soil.p_THICK;
    rows_SWAT_amt  = reduce(hcat, [sol(t).SWATI.mm   for t in saveat]) ./ sol.prob.p.p_soil.p_THICK
    rows_SWAT_d18O = reduce(hcat, [sol(t).SWATI.d18O for t in saveat])
    rows_SWAT_d2H  = reduce(hcat, [sol(t).SWATI.d2H  for t in saveat])
    row_NaN       = fill(NaN, 1,length(x))
    row_PREC_d18O = reshape(sol.prob.p.p_δ18O_PREC.(saveat), 1, :)
    row_INTS_d18O = reduce(hcat, [sol(t).INTS.d18O for t in saveat])
    row_INTR_d18O = reduce(hcat, [sol(t).INTR.d18O for t in saveat])
    row_SNOW_d18O = reduce(hcat, [sol(t).SNOW.d18O for t in saveat])
    row_GWAT_d18O = reduce(hcat, [sol(t).GWAT.d18O for t in saveat])
    row_RWU_d18O  = reduce(hcat, [sol(t).RWU.d18O for t in saveat])
    row_XYL_d18O  = reduce(hcat, [sol(t).XYLEM.d18O for t in saveat])
    row_PREC_d2H  = reshape(sol.prob.p.p_δ2H_PREC.(saveat), 1, :)
    row_INTS_d2H  = reduce(hcat, [sol(t).INTS.d2H for t in saveat])
    row_INTR_d2H  = reduce(hcat, [sol(t).INTR.d2H for t in saveat])
    row_SNOW_d2H  = reduce(hcat, [sol(t).SNOW.d2H for t in saveat])
    row_GWAT_d2H  = reduce(hcat, [sol(t).GWAT.d2H for t in saveat])
    row_RWU_d2H   = reduce(hcat, [sol(t).RWU.d2H for t in saveat])
    row_XYL_d2H   = reduce(hcat, [sol(t).XYLEM.d2H for t in saveat])

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
    # link := :x #TODO(bernhard): link doesn't seem to work...
    layout --> RecipesBase.grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4])
    #layout --> (4,1)
    size --> (1000,1400)
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
    @series begin # ts_PREC_δ18O
        title := "δ18O"
        seriestype := :bar
        # linecolor := :match
        line_z := reshape(row_PREC_d18O,1,:)
        fill_z := reshape(row_PREC_d18O,1,:)
        clims := clims_d18O
        colorbar_title := "δ18O [‰]"
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
        colorbar_title := "δ2H [‰]"
        colorbar := true_to_check_colorbar # TODO: define this once for all plots (except force it for colorbar plot representing common legend)
        yguide := "PREC [mm]"
        legend := false
        subplot := 3

        # and other arguments:
        reshape(x,1,:), reshape(row_PREC_amt,1,:)
    end

    # 3b) Heatmap (containing SWATI and other compartments)
    # y_labels   = ["INTS"; ""; "INTR"; ""; "SNOW"; ""; round.(y_center); "";             "GWAT"]
    # y_soil_ticks = optimize_ticks(extrema(y_center)...; k_min = 4)[1]
    # y_soil_ticks = optimize_ticks(0., round(maximum(cumsum(sol.prob.p.p_soil.p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()

    y_extended = [-500; -350; -300; -250; -200; -150; -100; -50;         y_center;             (maximum(cumsum(sol.prob.p.p_soil.p_THICK)) .+ [50; 100; 150; 250; 300;400])]
    y_soil_ticks = tick_function(0., round(maximum(cumsum(sol.prob.p.p_soil.p_THICK))))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
    y_ticks    = [-500;       -300;       -200;       -100;          y_soil_ticks;             (maximum(cumsum(sol.prob.p.p_soil.p_THICK)) .+ [    100;      250;     400])]
    y_labels   = ["PREC";   "INTS";     "INTR";     "SNOW";     round.(y_soil_ticks; digits=0);                                                "GWAT";    "RWU";     "XYLEM"]
    z2_extended = [row_PREC_d18O; row_NaN; row_INTS_d18O; row_NaN; row_INTR_d18O; row_NaN; row_SNOW_d18O; row_NaN; rows_SWAT_d18O; row_NaN; row_GWAT_d18O; row_NaN; row_RWU_d18O; row_NaN; row_XYL_d18O]
    z3_extended = [row_PREC_d2H;  row_NaN; row_INTS_d2H;  row_NaN; row_INTR_d2H;  row_NaN; row_SNOW_d2H;  row_NaN; rows_SWAT_d2H;  row_NaN; row_GWAT_d2H;  row_NaN; row_RWU_d2H;  row_NaN; row_XYL_d2H ]

    @series begin # pl_δ18O
        seriestype := :heatmap
        yflip := true
        yticks := (y_ticks, y_labels)
        colorbar := true_to_check_colorbar
        yguide := "Depth [mm]"
        colorbar_title := "δ18O [‰]"
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
        colorbar_title := "δ2H [‰]"
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

end
