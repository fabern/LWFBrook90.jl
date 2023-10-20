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

using ..PET: ESAT

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

"""
    compute_isotope_U_of_INTS_INTR_SNOW_and_SLFL(
        p_δ2H_PREC, p_δ18O_PREC, p_fT_TADTM, p_fT_VAPPRES,
        # for INTS (in: SINT; out: ISVP):
        u_INTS, aux_du_SINT, aux_du_ISVP, p_DTP, u_δ2H_INTS, u_δ18O_INTS,
        # for INTR (in: RINT; out: IRVP):
        u_INTR, aux_du_RINT, aux_du_IRVP, u_δ2H_INTR, u_δ18O_INTR,
        # for SNOW (in: STHR, RSNO (both δ_PREC); out: SMLT, SNVP (δ_SNOW and fractionated)):
        p_fu_STHR, aux_du_RSNO, aux_du_SMLT, aux_du_SNVP, u_δ2H_SNOW, u_δ18O_SNOW,
        # to compute isotopic signature of soil infiltration: SLFL
        p_fu_RNET)

Computes updated values of states INTS, INTR, and SNOW as well as their isotopic composition.
Compute mixing and evaporative fractionation, and also compute the isotopic composotion of
the resulting flux that infiltrates into the soil: p_fu_δ18O_SLFL, p_fu_δ2H_SLFL

The function is called in the daily callback.

The function returns: p_fu_δ18O_SLFL, p_fu_δ2H_SLFL, as well as:
    u_INTS
    u_δ18O_INTS
    u_δ2H_INTS
    u_INTR
    u_δ18O_INTR
    u_δ2H_INTR
    u_SNOW
    u_δ18O_SNOW
    u_δ2H_SNOW
"""
function compute_isotope_U_of_INTS_INTR_SNOW_and_SLFL(
    p_δ2H_PREC, p_δ18O_PREC, p_fT_TADTM, p_fT_VAPPRES,
    # for INTS (in: SINT; out: ISVP):
    u_INTS, aux_du_SINT, aux_du_ISVP, p_DTP, u_δ2H_INTS, u_δ18O_INTS,
    # for INTR (in: RINT; out: IRVP):
    u_INTR, aux_du_RINT, aux_du_IRVP, u_δ2H_INTR, u_δ18O_INTR,
    # for SNOW (in: STHR, RSNO (both δ_PREC); out: SMLT, SNVP (δ_SNOW and fractionated)):
    u_SNOW, p_fu_STHR, aux_du_RSNO, aux_du_SMLT, aux_du_SNVP, u_δ2H_SNOW, u_δ18O_SNOW,
    # to compute isotopic signature of soil infiltration: SLFL
    p_fu_RNET)
    # check chart "../docs/src/assets/b90flow.gif"

    ##################
    ##################
    ##################
    # Define conditions for isotope calculations:
    Tc = p_fT_TADTM  # °C, average daytime air temperature
    h = min(1.0, p_fT_VAPPRES / ESAT(Tc)[1]) # -, relative humidity of the atmosphere (vappress_atm/1 atm)
    γ = 1.0          # -, thermodynamic activity coefficient of evaporating water
    X_INTS = 0.5  # -, turbulence incex of the atmosphere above the evaporating water
    X_INTR = 0.5  # -, turbulence incex of the atmosphere above the evaporating water
    X_SNOW = 1.0  # -, turbulence incex of the atmosphere above the evaporating water

    # 1c) Atmospheric vapor composition assumed to be in equilibrium with precipitation
    #     Benettin-2018-Hydrol_Earh_Syst_Sci citing Gibson et al. 2008
    δₐ(δ_p, α_eq) = (δ_p - 1000 * (α_eq - 1))/α_eq
    δ²H_a  = δₐ(p_δ2H_PREC,  ISO.α²H_eq(Tc))
    δ¹⁸O_a = δₐ(p_δ18O_PREC, ISO.α¹⁸O_eq(Tc))

    # # For soil:
    # Xa = 0.5 # between molecular and turbulent
    # Xs = 1.0 # molecular only
    # X_SOIL = ((θ_s-θ_res)*Xa + (θ_sat-θ_s)*Xs)/(θ_sat-θ_res)
    # # Taken from Zhou-2021-Environ_Model_Softw (X is called n_k there)
    ##################
    ##################
    ##################


                # # 2a) INTS (in: SINT*δ_SINT; out: ISVP*δ_ISVP)
                # #          with δ_SINT = δ_PREC; δ_ISVP = f(f, α, ...)
                # # Operator step 1
                # @assert aux_du_SINT >= 0 "aux_du_SINT should not be negative"
                # if ((u_INTS == 0) & (aux_du_SINT == 0)) # initially no intercepted snow and no new is added
                #     u_δ18O_INTS_final = 0
                #     u_δ2H_INTS_final  = 0
                #     u_INTS_final      = 0
                # else
                #     if ((u_INTS == 0) & (aux_du_SINT > 0)) # initially no intercepted snow but some is added
                #         u_δ18O_INTS_first = p_δ18O_PREC
                #         u_δ2H_INTS_first  = p_δ2H_PREC
                #     else # u_INTS is not zero, and some/none is added (aux_du_SINT=0 or aux_du_SINT>0)
                #         u_δ18O_INTS_first = u_δ18O_INTS + aux_du_SINT*p_DTP/u_INTS * (p_δ18O_PREC-u_δ18O_INTS)
                #         u_δ2H_INTS_first  = u_δ2H_INTS  + aux_du_SINT*p_DTP/u_INTS * (p_δ2H_PREC -u_δ2H_INTS )
                #     end
                #     u_INTS_first = u_INTS + aux_du_SINT*p_DTP
                #     # Operator step 2: now taking care of aux_du_ISVP
                #     u_INTS_final = u_INTS + (aux_du_SINT - aux_du_ISVP) * p_DTP
                #     f_INTS = max(0, u_INTS_final/u_INTS_first) # fraction remaining after evaporation
                #     # ε_δ18O = 1/(ISO.α¹⁸O_eq(Tc) * ISO.α¹⁸O_dif^X_INTS) - 1
                #     # ε_δ2H  = 1/(ISO.α²H_eq(Tc)  * ISO.α²H_dif^X_INTS)  - 1
                #     # u_δ18O_INTS_final = 1000 * ( (1 + u_δ18O_INTS_first/1000)(f_INTS)^ε_δ18O - 1 )
                #     # u_δ2H_INTS_final  = 1000 * ( (1 + u_δ2H_INTS_first /1000)(f_INTS)^ε_δ2H  - 1 )
                #     u_δ18O_INTS_final = u_δ18O_INTS_first#1000 * ISO.δ_CraigGordon.(u_δ18O_INTS_first/1000, δ¹⁸O_a/1000, f_INTS, h, ISO.α¹⁸O_eq(Tc), ISO.α¹⁸O_dif, γ, X_INTS)
                #     u_δ2H_INTS_final  = u_δ2H_INTS_first#1000 * ISO.δ_CraigGordon.(u_δ2H_INTS_first /1000,  δ²H_a/1000, f_INTS, h, ISO.α²H_eq(Tc),  ISO.α²H_dif,  γ, X_INTS)
                # end

                # # 2b) INTR (in: RINT*δ_RINT; out: IRVP*δ_IRVP)
                # #          with δ_RINT = δ_PREC; δ_IRVP = f(f, α, ...)
                # # Operator step 1
                # @assert aux_du_RINT >= 0 "aux_du_RINT should not be negative"
                # if ((u_INTR == 0) & (aux_du_RINT == 0)) # initially no intercepted rain and no new is added
                #     u_δ18O_INTR_final = 0
                #     u_δ2H_INTR_final  = 0
                #     u_INTR_final      = 0
                # else
                #     if ((u_INTR == 0) & (aux_du_RINT > 0)) # initially no intercepted rain but some is added
                #         u_δ18O_INTR_first = p_δ18O_PREC
                #         u_δ2H_INTR_first  = p_δ2H_PREC
                #     else # u_INTR is not zero, and some/none is added (aux_du_RINT=0 or aux_du_RINT>0)
                #         u_δ18O_INTR_first = u_δ18O_INTR + aux_du_RINT*p_DTP/u_INTR * (p_δ18O_PREC-u_δ18O_INTR)
                #         u_δ2H_INTR_first  = u_δ2H_INTR  + aux_du_RINT*p_DTP/u_INTR * (p_δ2H_PREC -u_δ2H_INTR )
                #     end
                #     u_INTR_first = u_INTR + aux_du_RINT*p_DTP
                #     # Operator step 2: now taking care of aux_du_IRVP
                #     u_INTR_final = u_INTR + (aux_du_RINT - aux_du_IRVP) * p_DTP
                #     f_INTR = max(0, u_INTR_final/u_INTR_first) # fraction remaining after evaporation
                #     # ε_δ18O = 1/(ISO.α¹⁸O_eq(Tc) * ISO.α¹⁸O_dif^X_INTR) - 1
                #     # ε_δ2H  = 1/(ISO.α²H_eq(Tc)  * ISO.α²H_dif^X_INTR)  - 1
                #     # u_δ18O_INTR_final = 1000 * ( (1 + u_δ18O_INTR_first/1000)(f_INTR)^ε_δ18O - 1 )
                #     # u_δ2H_INTR_final  = 1000 * ( (1 + u_δ2H_INTR_first /1000)(f_INTR)^ε_δ2H  - 1 )
                #     u_δ18O_INTR_final = u_δ18O_INTR_first# 1000 * ISO.δ_CraigGordon.(u_δ18O_INTR_first/1000, δ¹⁸O_a/1000, f_INTR, h, ISO.α¹⁸O_eq(Tc), ISO.α¹⁸O_dif, γ, X_INTR)
                #     u_δ2H_INTR_final  = u_δ2H_INTR_first# 1000 * ISO.δ_CraigGordon.(u_δ2H_INTR_first /1000,  δ²H_a/1000, f_INTR, h, ISO.α²H_eq(Tc),  ISO.α²H_dif,  γ, X_INTR)
                # end

                # # 2c) SNOW (in: STHR*δ_STHR, RSNO*δ_RSNO; out: SMLT*δ_SMLT, SNVP*δ_SNVP)
                # #          with δ_STHR = δ_PREC, δ_RSNO = δ_PREC; δ_SMLT = δ_SNOW, δ_SNVP = f(f, α, ...)

                # # NOTE: for SNOW the isotope balance is greatly simplified. The most precise
                # #       approach would be to define mixing concentrattion within `SNOWPACK()` in `module_SNO.jl`
                # # Operator step 1
                # @assert p_fu_STHR + aux_du_RSNO >= 0 "p_fu_STHR + aux_du_RSNO should not be negative"
                # if ((u_SNOW == 0) & (p_fu_STHR + aux_du_RSNO == 0)) # initially no snowpack and no new is added
                #     u_δ18O_SNOW_final = 0
                #     u_δ2H_SNOW_final  = 0
                #     u_SNOW_final      = 0
                # else
                #     if ((u_SNOW == 0) & (p_fu_STHR + aux_du_RSNO > 0)) # initially no snowpack but some is added
                #         u_δ18O_SNOW_first = p_δ18O_PREC
                #         u_δ2H_SNOW_first  = p_δ2H_PREC
                #     else # u_SNOW is not zero, and some/none is added ((p_fu_STHR + aux_du_RSNO)=0 or (p_fu_STHR + aux_du_RSNO)>0)
                #         # p_fu_STHR, aux_du_RSNO, aux_du_SMLT, aux_du_SNVP
                #         u_δ18O_SNOW_first = u_δ18O_SNOW + (p_fu_STHR + aux_du_RSNO)*p_DTP/u_SNOW * (p_δ18O_PREC - u_δ18O_SNOW)
                #                                     # NOTE: because the outflow term for aux_du_SMLT has
                #                                     #       an isotope concentration of u_δ18O_SNOW and it is thus not needed:
                #                                     # - (aux_du_SMLT)*p_DTP/u_SNOW * (u_δ18O_SNOW-u_δ18O_SNOW)
                #         u_δ2H_SNOW_first = u_δ2H_SNOW + (p_fu_STHR + aux_du_RSNO)*p_DTP/u_SNOW * (p_δ2H_PREC - u_δ2H_SNOW)
                #                                     # NOTE: because the outflow term for aux_du_SMLT has
                #                                     #       an isotope concentration of u_δ2H_SNOW and it is thus not needed:
                #                                     # - (aux_du_SMLT)*p_DTP/u_SNOW * (u_δ2H_SNOW-u_δ2H_SNOW)
                #     end
                #     u_SNOW_first = u_SNOW + p_DTP * (p_fu_STHR + aux_du_RSNO - aux_du_SMLT)
                #     # Operator step 2: now taking care of aux_du_SNVP
                #     u_SNOW_final = u_SNOW + p_DTP * (p_fu_STHR + aux_du_RSNO - aux_du_SMLT - aux_du_SNVP)
                #     f_SNOW = max(0, u_SNOW_final/u_SNOW_first) # fraction remaining after evaporation
                #     # ε_δ18O = 1/(ISO.α¹⁸O_eq(Tc) * ISO.α¹⁸O_dif^X_SNOW) - 1
                #     # ε_δ2H  = 1/(ISO.α²H_eq(Tc)  * ISO.α²H_dif^X_SNOW)  - 1
                #     # u_δ18O_SNOW_final = 1000 * ( (1 + u_δ18O_SNOW_first/1000)(f_SNOW)^ε_δ18O - 1 )
                #     # u_δ2H_SNOW_final  = 1000 * ( (1 + u_δ2H_SNOW_first /1000)(f_SNOW)^ε_δ2H  - 1 )
                #     u_δ18O_SNOW_final = u_δ18O_SNOW_first#1000 * ISO.δ_CraigGordon.(u_δ18O_SNOW_first/1000, δ¹⁸O_a/1000, f_SNOW, h, ISO.α¹⁸O_eq(Tc), ISO.α¹⁸O_dif, γ, X_SNOW)
                #     u_δ2H_SNOW_final  = u_δ2H_SNOW_first#1000 * ISO.δ_CraigGordon.(u_δ2H_SNOW_first /1000,  δ²H_a/1000, f_SNOW, h, ISO.α²H_eq(Tc),  ISO.α²H_dif,  γ, X_SNOW)
                # end

    # 2a) u_INTS (in: aux_du_SINT*δ_SINT; out: aux_du_ISVP*δ_ISVP)
    #            with δ_SINT = δ_PREC; δ_ISVP = f(f, α, ...)
    # 2b) u_INTR (in: aux_du_RINT*δ_RINT; out: aux_du_IRVP*δ_IRVP)
    #            with δ_RINT = δ_PREC; δ_IRVP = f(f, α, ...)
    # 2c) SNOW (in: p_fu_STHR*δ_STHR, aux_du_RSNO*δ_RSNO; out: aux_du_SMLT*δ_SMLT, aux_du_SNVP*δ_SNVP)
    #          with δ_STHR = δ_PREC, δ_RSNO = δ_PREC; δ_SMLT = δ_SNOW, δ_SNVP = f(f, α, ...)





    # update_δ_with_mixing_and_evaporation(dt, u₀, δ₀, inflow, δin, outflow, E, δₐ, h, α_eq, α_dif, γ, X)
    u_INTS_final, u_δ2H_INTS_final  = ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_INTS, u_δ2H_INTS,  aux_du_SINT,               p_δ2H_PREC,                0,           ISO.R_VSMOW²H,  aux_du_ISVP, δ²H_a,  h, ISO.α²H_eq(Tc),  ISO.α²H_dif,  γ, X_INTS)
    _,            u_δ18O_INTS_final = ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_INTS, u_δ18O_INTS, aux_du_SINT,               p_δ18O_PREC,               0,           ISO.R_VSMOW¹⁸O, aux_du_ISVP, δ¹⁸O_a, h, ISO.α¹⁸O_eq(Tc), ISO.α¹⁸O_dif, γ, X_INTS)
    u_INTR_final, u_δ2H_INTR_final  = ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_INTR, u_δ2H_INTR,  aux_du_RINT,               p_δ2H_PREC,                0,           ISO.R_VSMOW²H,  aux_du_IRVP, δ²H_a,  h, ISO.α²H_eq(Tc),  ISO.α²H_dif,  γ, X_INTR)
    _,            u_δ18O_INTR_final = ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_INTR, u_δ18O_INTR, aux_du_RINT,               p_δ18O_PREC,               0,           ISO.R_VSMOW¹⁸O, aux_du_IRVP, δ¹⁸O_a, h, ISO.α¹⁸O_eq(Tc), ISO.α¹⁸O_dif, γ, X_INTR)
    # NOTE: for SNOW the isotope balance is greatly simplified. The most precise
    #       approach would be to define mixing concentrattion within `SNOWPACK()` in `module_SNO.jl`
    u_SNOW_final, u_δ2H_SNOW_final  = ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_SNOW, u_δ2H_SNOW,  (p_fu_STHR, aux_du_RSNO), (p_δ2H_PREC,  p_δ2H_PREC),  aux_du_SMLT, ISO.R_VSMOW²H,  aux_du_SNVP, δ²H_a,  h, ISO.α²H_eq(Tc),  ISO.α²H_dif,  γ, X_SNOW)
    _,            u_δ18O_SNOW_final = ISO.update_δ_with_mixing_and_evaporation(p_DTP, u_SNOW, u_δ18O_SNOW, (p_fu_STHR, aux_du_RSNO), (p_δ18O_PREC, p_δ18O_PREC), aux_du_SMLT, ISO.R_VSMOW¹⁸O, aux_du_SNVP, δ¹⁸O_a, h, ISO.α¹⁸O_eq(Tc), ISO.α¹⁸O_dif, γ, X_SNOW)


    δ18O_empty = NaN
    δ2H_empty  = NaN

    # 3) also compute δ_SLFL as mix of δ_SMLT with δ_RNET (i.e. water that infiltrates)
    if (aux_du_SMLT + p_fu_RNET == 0)
        p_fu_δ18O_SLFL = δ18O_empty
        p_fu_δ2H_SLFL  = δ2H_empty
    elseif (aux_du_SMLT == 0)
        p_fu_δ18O_SLFL = p_δ18O_PREC
        p_fu_δ2H_SLFL  = p_δ2H_PREC
    elseif (p_fu_RNET == 0)
        p_fu_δ18O_SLFL = u_δ18O_SNOW_final  # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
        p_fu_δ2H_SLFL  = u_δ2H_SNOW_final   # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
    else # both fluxes are non-null and we need to compute their mix
        p_fu_δ18O_SLFL = (u_δ18O_SNOW_final * aux_du_SMLT + p_δ18O_PREC * p_fu_RNET) / (aux_du_SMLT + p_fu_RNET)  # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
        p_fu_δ2H_SLFL  = (u_δ2H_SNOW_final * aux_du_SMLT  + p_δ2H_PREC * p_fu_RNET)  / (aux_du_SMLT + p_fu_RNET)  # TODO(bernhard): using final is effectively operator splitting, (the isotope mass balance is not exact)
    end
    # TODO(bernhard): deactivate the following workaround and activate above code
    # p_fu_δ18O_SLFL = p_δ18O_PREC
    # p_fu_δ2H_SLFL  = p_δ2H_PREC

    return (p_fu_δ18O_SLFL, p_fu_δ2H_SLFL,
        u_INTS_final, u_δ18O_INTS_final, u_δ2H_INTS_final,
        u_INTR_final, u_δ18O_INTR_final, u_δ2H_INTR_final,
        u_SNOW_final, u_δ18O_SNOW_final, u_δ2H_SNOW_final)

    # In-place modify
    # u_INTS      = u_INTS_final
    # u_δ18O_INTS = u_δ18O_INTS_final
    # u_δ2H_INTS  = u_δ2H_INTS_final
    # u_INTR      = u_INTR_final
    # u_δ18O_INTR = u_δ18O_INTR_final
    # u_δ2H_INTR  = u_δ2H_INTR_final
    # u_SNOW      = u_SNOW_final
    # u_δ18O_SNOW = u_δ18O_SNOW_final
    # u_δ2H_SNOW  = u_δ2H_SNOW_final
end

end
