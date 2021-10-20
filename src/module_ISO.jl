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


export α¹⁸O_dif, α²H_dif, α¹⁸O_eq, α²H_eq, δ_CraigGordon
export plot_LWFBrook90Iso

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


end
