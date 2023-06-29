@doc raw"""
# Soil water properties

Text copied from Ecoshift on module KPT:

"
This section describes the soil water properties in BROOK90, and the four routines used to
calculate their values. The section also discusses the concept of field capacity, and how
it is or is not used. Finally the section comments on the selection of soil water
parameters and the use of the menu item EditParameters-KPTGraph. This section uses both
the standard algebraic notation for soil water variables as in Clapp and Hornberger(1978)
and the BROOK90 variable names. The correspondence is
θ  = THETA   W  = WETNES  ψ  = PSIM    K  = KK    b = BEXP
θf = THETAF  Wf = WETF    ψf = PSIF    Kf = KF    n = CHN
θs = THSAT   Wi = WETINF  ψi = PSIINF  Ks = KSAT  m = CHM
T  = THICK   S  = STONEF

Ecoshift-KPT:
Functional relationships among soil-water content, θ, matric potential, ψ,
and hydraulic conductivity, K, are required for any simulation that moves water through
the soil. There is a rather vast literature on this subject, though much of it emphasizes
agricultural soils rather than natural soils, which are less disturbed and higher in
organic matter.

##   Clapp Hornberger (FLAG_MualVanGen=0) see: http://www.ecoshift.net/brook/kpt.html
BROOK90 uses a modification of the Campbell (1974) expressions with the near-saturation
interpolation of Clapp and Hornberger (1978).
1) W = (θ - θr) / θs - θr)
2) Ψ = Ψs*W^(-b)              with b=1/λ
3) K = Ks*W^(2b+3)            with b=1/λ
near saturated:
4) Ψ = -m*(W - n)*(1 - W)     = (-mW + mn) +m*W^2 -Wmn = m*W^2 -m*(1+n)*W + m*n
5)      m = -Ψi/(1 - Wi) * [1/(1-Wi) - b/Wi]
6)      n = 2 Wi - 1 + (b*Ψi / (m*Wi))  ### where Ψi is Ψ at Wi obtained from 2)
clearly unsaturated:
7) Ψ = Ψf*(W/Wf)^(-b)       = Ψf*Wf^(b)*W^(-b)
8) K = Kf*(W/Wf)^(2b+3)  ### where Wf = θf/θs is the wetness at θf (with θr=0)

   ==> dψdW = Ψf*Wf^(b)*    (-b*W^(-b-1))    # in clearly unsaturated range
   ==> dψdW = 2*m*W - m*(1+n)

##   Mualem van Genuchten (FLAG_MualVanGen=1)
1) W = 1/(1 + (αh)^n )^m                         (Radcliffe eq.2.47)
2) ψ = -1/α * [ (1-W^(1/m))*(W^(-1/m)) ]^(1/n)   (Radcliffe p.122)
...
8) K = Ks*W^l*[ 1 - (1-W^(1/m))^m ]^2

   ==> dψdW = .....
"

"""
module KPT # SOIL WATER PROPERTIES

export SOILPAR_CH, derive_auxiliary_SOILVAR
export KPT_SOILPAR_Mvg1d, KPT_SOILPAR_Ch1d, FWETNES
export AbstractSoilHydraulicParams, MualemVanGenuchtenSHP
using Roots: find_zero, Bisection # to find wetness for a given hydraulic conductivity
using ..CONSTANTS: p_RHOWG  # https://discourse.julialang.org/t/large-programs-structuring-modules-include-such-that-to-increase-performance-and-readability/29102/5
using DataFrames: DataFrame
using Printf: @sprintf
using UnPack: @unpack

### ### Parameters
### #   Clapp and Hornberger (FLAG_MualVanGen=0)
### #   Par_model0 = [THSAT,THETAF,KF,PSIF,WETF,KSAT,CHM,CHN,BEXP,WETINF]
### # from LWFBrook90R: !   real(kind=8) :: PSIF     ! matrix potential at field capacity, kPa
### # from LWFBrook90R:     real(kind=8) :: BEXP     ! exponent for psi-theta relation
### # from LWFBrook90R: !   real(kind=8) :: WETINF   ! wetness at dry end of near-saturation range
### # from LWFBrook90R:     real(kind=8) :: WETF     ! saturation fraction at field capacity
### # from LWFBrook90R: !   real(kind=8) :: CHM      ! Clapp and Hornberger m, kPa
### # from LWFBrook90R: !   real(kind=8) :: CHN      ! Clapp and Hornberger n
### # from LWFBrook90R:     real(kind=8) :: KF       ! hydraulic conductivity at field capacity, mm/d
### #   Mualem van Genuchten (FLAG_MualVanGen=1)
### #   Par_model1 = [θs  ,THETAF,?, PSIF,WETF,KSAT,MvGα,MvGn,MvGl   ,θr]
### # from LWFBrook90R: !   real(kind=8) :: θs      ! water content at saturation
### # from LWFBrook90R: !   real(kind=8) :: θr      ! residual water content
### # from LWFBrook90R:     real(kind=8) :: MvGα     ! parameter alpha, 1/m
### # from LWFBrook90R:     real(kind=8) :: MvGn     ! parameter n
### # from LWFBrook90R:     real(kind=8) :: KSAT     ! hydraulic conductivity at saturation, mm/d
### # from LWFBrook90R:     real(kind=8) :: MvGl     ! tortuosity parameter (default = 0.5)


### Define Julia types for soil parametrization:

# brief INTRO to custom types: https://syl1.gitbook.io/julia-language-a-concise-tutorial/language-core/custom-types
# medium INTRO: https://scls.gitbooks.io/ljthw/content/_chapters/06-ex3.html
# long INTRO: https://benlauwens.github.io/ThinkJulia.jl/latest/book.html#chap15
# DOC: https://docs.julialang.org/en/v1/manual/types/#Composite-Types
# Careful: https://docs.julialang.org/en/v1/manual/performance-tips/#The-dangers-of-abusing-multiple-dispatch-(aka,-more-on-types-with-values-as-parameters)
# Careful: Don't use abstract types for fields (p_THICK::AbstractVector). Use {T} https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-fields-with-abstract-type
#
# Because we want to derive additional fields at initialization, an inner construcor is needed.
# In order that it works with parametric types ( d{T<:Number} ), an outer constructer that derives
# the type from its elements is needed. Without this one would need to provide the type explicitly
# when instantiating e.g. d{Float64}(1.2, 2.0) instead of d(1.2, 2.0).
# struct d{T<:Number}
#     # Input fields
#     field1::T
#     field2::T
#     field3::T

#     function d{T}(field1, field2) where {T<:Number}
#         new(field1, field2, field1+field2)
#     end
# end
# d(field1::T, field2::T) where {T<:Number} = d{T}(field1,field2) # see https://docs.julialang.org/en/v1/manual/constructors/#Parametric-Constructors

"""
    AbstractSoilHydraulicParams

Represents an abstract parametrization of soil hydraulics.
Examples of soil hydraulic parametrizations are:
- Mualem-van Genuchten
- Clapp-Hornberger

A summary is available in
Shao, Y. and Irannejad, P.: On the Choice of Soil Hydraulic Models in Land-Surface Schemes, Boundary Layer Meterol., 90, 83–115, https://doi.org/10.1023/A:1001786023282, 1999.
"""
abstract type AbstractSoilHydraulicParams end # The concrete types are defined in module KPT

@doc raw"""
    MualemVanGenuchtenSHP

Mualem-van Genuchten parametrization of soil hydraulics

```math
\Theta    = \frac{\theta - \theta_r}{\theta_s - \theta_r} = \left( \frac{1}{1+(-\alpha \psi)^n} \right)^m \\
n         = \frac{1}{m-1} \\
K(\theta) = K_s \Theta^{1/2}\left[ 1 - (1 - \Theta^{1/m})^m \right ]^2 \\
K(\psi)   = K_s\frac{\left[ 1- (-\alpha\psi)^{n-1} (1 + (-\alpha \psi)^n)^{-m} \right ]^2}{\left[ 1 + (-\alpha\psi)^n \right ]^{m/2}}
```
source: Shao, Y. and Irannejad, P.: On the Choice of Soil Hydraulic Models in Land-Surface Schemes, Boundary Layer Meterol., 90, 83–115, https://doi.org/10.1023/A:1001786023282, 1999.
"""
mutable struct MualemVanGenuchtenSHP <: AbstractSoilHydraulicParams
    # Soil hydraulic parameters: Mualem-van Genuchten
    # Input fields
    "Saturation volumetric soil water content [m³ m⁻³]"
    p_THSAT::Real
    "Residual volumetric soil water content [m³ m⁻³]"
    p_θr::Real
    "Mualem-van Genuchten α [m⁻¹]"
    p_MvGα::Real
    "Mualem-van Genuchten n [-]"
    p_MvGn::Real
    "Mualem-van Genuchten n [-]"
    p_MvGm::Real
    "Saturated hydraulic conductivity [mm day⁻¹]"
    p_KSAT::Real
    "Tortuosity [-]"
    p_MvGl::Real
    "Gravel/stone fraction [m³ m⁻³]"
    p_STONEF::Real
end
function MualemVanGenuchtenSHP(;p_THSAT,p_θr,p_MvGα,p_MvGn,p_KSAT,p_MvGl,p_STONEF)
    MualemVanGenuchtenSHP(
        p_THSAT,
        p_θr,
        p_MvGα,
        p_MvGn,
        1 .- 1 ./ p_MvGn, # = p_MvGm
        p_KSAT,
        p_MvGl,
        p_STONEF)
end
function MualemVanGenuchtenSHP(df::DataFrame)
    [MualemVanGenuchtenSHP(
        p_THSAT  = dfrow.ths_volFrac,     p_θr     = dfrow.thr_volFrac,
        p_MvGα   = dfrow.alpha_perMeter,  p_MvGn   = dfrow.npar_,
        p_KSAT   = dfrow.ksat_mmDay,      p_MvGl   = dfrow.tort_,
        p_STONEF = dfrow.gravel_volFrac)
    for dfrow in eachrow(df)]
end

function Base.show(io::IO, shp::MualemVanGenuchtenSHP)
    # function Base.show(io::IO, mime::MIME"text/plain", shp::MualemVanGenuchtenSHP)
    print(io,
        @sprintf("(θ from θr=%.3f to θs=%.3f, Ks =% 8.1f, STONEF=%.1f, l=%.1f, n=%.1f, α=% 7.1f)",
                shp.p_θr, shp.p_THSAT,
                shp.p_KSAT,
                shp.p_STONEF,
                shp.p_MvGl,
                shp.p_MvGn,
                shp.p_MvGα))
end

@doc raw"""
    ClappHornbergerSHP

Clapp-Hornberger parametrization of soil hydraulics

```math
\frac{\psi}{\psi_s} = w^{-1/\lambda} \qquad \textup{if: }\psi \leq \psi_i \\
\psi                = -m(w-n)(w-1) \qquad \textup{if: }\psi_i \leq \psi \leq 0 \\
K(\theta)           = K_s \left( \frac{\theta}{\theta_s} \right)^{(3+2/\lambda)} \\
\\
\textup{Parameters:} \\
m = \frac{\psi_i}{(1-w_i)^2} - \frac{\psi_i}{w_i(1-w_i)\lambda} \\
n = 2w_i - \frac{\psi_i}{mw_i\lambda} -1
```
source: Shao, Y. and Irannejad, P.: On the Choice of Soil Hydraulic Models in Land-Surface Schemes, Boundary Layer Meterol., 90, 83–115, https://doi.org/10.1023/A:1001786023282, 1999.
"""
struct ClappHornbergerSHP <: AbstractSoilHydraulicParams
    # Input fields
    p_THSAT::Real
    p_THETAF::Real
    p_KF::Real
    p_PSIF::Real
    p_BEXP::Real
    p_WETINF::Real
    p_STONEF::Real
    # TODO: finalize by commenting (see also l.246 in "definition_p.jl")
end
# TODO(bernhard): Rewrite functions below for AbstractSoilHydraulicParams and specialize
# derive_auxiliary_SOILVAR(shp::AbstractSoilHydraulicParams)
# derive_auxiliary_SOILVAR(shp::MualemVanGenuchtenSHP)
# derive_auxiliary_SOILVAR(shp::ClappHornbergerSHP)
# FK(shp::AbstractSoilHydraulicParams, θ)
# FK(shp::MualemVanGenuchtenSHP, θ)
# FK(shp::ClappHornbergerSHP, θ)

abstract type AbstractKptSoilpar end
"""
Represents a discretized 1D column of soil with Clapp-Hornberger parametrization.

Input fields: p_THICK, p_STONEF, p_THSAT, p_PSIF, p_THETAF, p_KF, p_BEXP, p_WETINF
Derived fields: NLAYER, p_CHM, p_CHN, p_THETAF, p_PSIG, p_SWATMAX, p_WETF
"""
struct KPT_SOILPAR_Ch1d{T<:AbstractVector} <: AbstractKptSoilpar
    # Input fields
    p_THICK::T
    p_STONEF::T
    p_THSAT::T
    p_PSIF::T
    p_THETAF::T
    p_KF::T
    p_BEXP::T
    p_WETINF::T
    # Derived fields
    NLAYER::Int64
    p_CHM::T
    p_CHN::T
    p_PSIG::T
    p_SWATMAX::T
    p_WETF::T

    # # Inner constructor:
    function KPT_SOILPAR_Ch1d{T}(;p_THICK, p_STONEF, p_THSAT, p_PSIF, p_THETAF, p_KF, p_BEXP, p_WETINF) where {T<:AbstractVector}
        NLAYER = length(p_THICK)
        @assert size(p_THICK) == size(p_STONEF) == size(p_THSAT) == size(p_PSIF) ==
            size(p_THETAF) == size(p_KF) == size(p_BEXP) == size(p_WETINF)

        # Derive fields
        # Variant A: used to be SOILPAR()
        # # removed...

        # Variant B: rewritten SOILPAR() here
        p_WETF   = p_THETAF ./ p_THSAT # wetness at field capacity, dimensionless
        p_KSAT   = p_KF .* (1 ./ p_WETF) .^ (2 .* p_BEXP .+ 3) # saturated hydraulic conductivity, mm/d
        p_SWATMAX = p_THICK .* p_THSAT .* (1.0 .- p_STONEF) # maximum water storage for layer, mm
        p_PSIG   = fill(NaN, NLAYER) # gravity potential negative down from surface, kPa
        for i = 1:NLAYER
            if i == 1
                p_PSIG[1] = -p_RHOWG * p_THICK[1] / 2
            else
                p_PSIG[i] = p_PSIG[i - 1] - p_RHOWG * ((p_THICK[i - 1] + p_THICK[i]) / 2)
            end
        end
        p_CHM    = fill(NaN, NLAYER) # Clapp and Hornberger m, kPa
        p_CHN    = fill(NaN, NLAYER) # Clapp and Hornberger n, dimensionless
        for i = 1:NLAYER
            # p_ψ∞: potential at dry end of near saturation range (kPa)
            p_ψ∞     = p_PSIF[i] * (p_WETINF[i] / p_WETF[i]) ^ -p_BEXP[i]
            p_CHM[i] = (-p_ψ∞ / (1 - p_WETINF[i]) ^ 2) - p_BEXP[i] * (-p_ψ∞) / (p_WETINF[i] * (1 - p_WETINF[i]))
            p_CHN[i] = 2 * p_WETINF[i] - 1 - (-p_ψ∞ * p_BEXP[i] / (p_CHM[i] * p_WETINF[i]))
        end
        # unused p_WETc = p_WETF .* (1000 .* p_PSICR ./ p_PSIF) .^ (-1 ./ p_BEXP) # wetness at p_PSICR, dimensionless

        # NOTE(bernhard): removed further tasks that were in SOILPAR() in LWFBrook90R
        #                  - computation of p_WETc = p_WETF .* (1000 .* p_PSICR ./ p_PSIF) .^ (-1 ./ p_BEXP) # wetness at p_PSICR, dimensionless
        #                  - test of validity of PSIM (assert that PSIM <= 0, i.e. error("matrix psi must be negative or zero"))
        #                  - computation of u_aux_WETNES:
        #                       if: PSIM==0 -> u_aux_WETNES = 1
        #                       else: u_aux_WETNES = p_WETF .* (u0_aux_PSIM ./ p_PSIF).^(1/.BEXP)
        #                             if u_aux_WETNES[i]>p_WETINF[i]: u_aux_WETNES[i] = (1 + p_CHM[i])/2 + 0.5*sqrt(p_CHN[i]^2 - 2*p_CHN[i] + 1 + 4 * u0_aux_PSIM[i]/p_CHM[i])
        #                  - computation of u_aux_SWATI = u_aux_WETNES .* p_SWATMAX
        #                  - computation of p_PsiCrit = FPSIM_CH(p_ThCrit./p_THSAT, p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN)

        # Instantiate
        new(p_THICK, p_STONEF, p_THSAT, p_PSIF, p_THETAF, p_KF, p_BEXP, p_WETINF,
            NLAYER, p_CHM, p_CHN, p_PSIG, p_SWATMAX, p_WETF)
    end
end
KPT_SOILPAR_Ch1d(;p_THICK::T, p_STONEF::T, p_THSAT::T, p_PSIF::T, p_THETAF::T, p_KF::T, p_BEXP::T, p_WETINF::T) where {T<:AbstractVector} =
    KPT_SOILPAR_Ch1d{T}(;p_THICK, p_STONEF, p_THSAT, p_PSIF, p_THETAF, p_KF, p_BEXP, p_WETINF)
    # for explanation see https://docs.julialang.org/en/v1/manual/constructors/#Parametric-Constructors

"""
Represents a discretized 1D column of soil with Mualem-van Genuchten parametrization.

Input fields: p_THICK, p_STONEF, p_THSAT, p_Kθfc, p_KSAT, p_MvGα, p_MvGn, p_MvGl, p_θr
Derived fields: p_PSIF, p_THETAF, p_PSIG, p_SWATMAX, p_WETF
"""
struct KPT_SOILPAR_Mvg1d{T<:AbstractVector} <: AbstractKptSoilpar
    # Input fields
    p_THICK::T
    p_STONEF::T
    p_THSAT::T
    p_Kθfc::T
    p_KSAT::T
    p_MvGα::T
    p_MvGn::T
    p_MvGm::T
    p_MvGl::T
    p_θr::T
    # Derived fields
    NLAYER::Int64
    p_PSIF::T    # matric potential at field capacity, kPa
    p_THETAF::T  # soil moisture θ at field capacity, m3/m3
    p_PSIG::T    # gravity potential negative down from surface, kPa
    p_SWATMAX::T  # maximum water storage for layer, mm
    p_WETF::T    # wetness at field capacity, dimensionless


    # # Inner constructor:
    function KPT_SOILPAR_Mvg1d{T}(;p_THICK, p_STONEF, p_THSAT, p_Kθfc, p_KSAT, p_MvGα, p_MvGn, p_MvGm, p_MvGl, p_θr) where {T<:AbstractVector}
        NLAYER = length(p_THICK)
        @assert size(p_THICK) == size(p_STONEF) == size(p_THSAT) == size(p_Kθfc) ==
                size(p_KSAT) == size(p_MvGα) == size(p_MvGn) == size(p_MvGm) ==
                size(p_MvGl) == size(p_θr)
        @assert !any(isnan.(p_Kθfc))

        # # Variant A: used to be SOILPAR()
        # # removed...

        # Variant B: rewritten SOILPAR() here
        p_SWATMAX = p_THICK .* p_THSAT .* (1.0 .- p_STONEF) # maximum water storage for layer, mm
        p_PSIG   = fill(NaN, NLAYER) # gravity potential negative down from surface, kPa
        for i = 1:NLAYER
            if i == 1
                p_PSIG[1] = -p_RHOWG * p_THICK[1] / 2
            else
                p_PSIG[i] = p_PSIG[i - 1] - p_RHOWG * ((p_THICK[i - 1] + p_THICK[i]) / 2)
            end
        end
        p_WETF   = fill(NaN, NLAYER) # wetness at field capacity, dimensionless
        for i = 1:NLAYER
            # Find wetness at field capacity based on hydraulic conductivity
            # using a nonlinear solver for: find θ such that: K(θ) - p_Kθfc = 0
            # Actually not finding θ but wetness = θ/porosity.
            # plot(-0.15:0.01:1.0, FK_MvG.(-0.15:0.01:1.0, p_KSAT[i], p_MvGl[i], p_MvGn[i]))
            p_WETF[i] = try
                f_root(x) = - p_Kθfc[i] + FK_MvG(x, p_KSAT[i], p_MvGl[i], p_MvGn[i], p_MvGm[i])
                find_zero(f_root, (0.0, 1.0), Bisection())  # find_zero((x) -> x, (0.1, 1.0), Bisection())
                # In LWFBrook90R: this was done using FWETK()
                # In LWFBrook90R: if p_WETF[i] == -99999. error("Computed invalid p_WETF by FWETK()") end
            catch e
                if e isa ArgumentError
                    function θ(PSIM_toFind) p_θr[i] + (p_THSAT[i] - p_θr[i])*(1+(-p_MvGα[i]*PSIM_toFind/9.81)^p_MvGn[i])^(-(1-1/p_MvGn[i])) end
                    # function ψ(TH_toFind) 9.81 .* (-1 ./ p_MvGα[i]) .* (max.((TH_toFind - p_θr[i])/(p_THSAT[i] .- p_θr[i]),  1.e-6) .^ (-1 ./ (1 .- 1 ./ p_MvGn[i])) .-1) .^ (1 ./ p_MvGn[i]) end
                    # function θ2(PSIM_toFind) p_θr[i] + (p_THSAT[i] - p_θr[i]) * LWFBrook90.KPT.FWETNES(u_aux_PSIM, p_soil) end
                    # function ψ2(TH_toFind) LWFBrook90.KPT.FPSIM_MvG((TH_toFind - p_θr[i])/(p_THSAT[i] .- p_θr[i]), p_MvGα[i], p_MvGn[i]) end
                    # plot(    -1:-1:-100,      θ.(-1:-1:-100), xlabel = "PSIM (kPa)", ylabel = "θ (-)", label = "function θ")
                    # plot!(ψ.(0.48:0.01:0.60), 0.48:0.01:0.60, xlabel = "PSIM (kPa)", ylabel = "θ (-)", label = "function ψ")
                    # plot!(ψ2.(0.48:0.01:0.60), 0.48:0.01:0.60, xlabel = "PSIM (kPa)", ylabel = "θ (-)", label = "function ψ2")
                    ψ_fc = -33
                    @warn "When setting up soil hydr. parameters: $e"*
                        "\n\n"*"Ignoring hardcoded p_Kθfc ($(p_Kθfc[i]) mm/day) and"*
                        " using ψ = $ψ_fc kPa for field capacity\n\n "*
                        " Layer:$(i), $(p_Kθfc[i]), $(p_KSAT[i]), $(p_MvGl[i]), $(p_MvGn[i])"
                    θ(ψ_fc)/(p_THSAT[i] .- p_θr[i])
                    # find_zero((θ_toFind) -> LWFBrook90.KPT.FPSIM_MvG((θ_toFind - p_θr[i])./(p_THSAT .- p_θr), p_MvGα[i], p_MvGn[i]), (0.0, 1.0), Bisection())
                else
                    error("When setting up soil hydr. parameters: Computed invalid p_WETF in KPT_SOILPAR_Mvg1d")
                end
            end
        end

        p_PSIF   = FPSIM_MvG(p_WETF, p_MvGα, p_MvGn, p_MvGm)  # matric potential at field capacity, kPa
        p_THETAF = FTheta_MvG(p_WETF, p_THSAT, p_θr)  # soil moisture θ at field capacity, m3/m3
        p_Kθfc .= NaN                                 # p_Kθfc is only used as input parameter for setup, not for calculation

        # NOTE(bernhard): removed further tasks that were in SOILPAR() in LWFBrook90R
        #                  - computation of p_WETc = FWETNES_MvG.(1000 .* p_PSICR, p_MvGα, p_MvGn) # wetness at p_PSICR, dimensionless
        #                  - test of validity of PSIM (assert that PSIM <= 0, i.e. error("matrix psi must be negative or zero"))
        #                  - derivation of state dependent parameters u_aux_WETNES and u_aux_SWATI
        #                       - u_aux_WETNES[i] = FWETNES_MvG(u_aux_PSIM[i], p_MvGα[i], p_MvGn[i])
        #                       - u_aux_SWATI[i] = FTheta_MvG(u_aux_WETNES[i], p_THSAT[i], p_θr[i]) * p_SWATMAX[i]/p_THSAT[i]
        #                  - computation of p_PsiCrit = FPSIM_CH(p_ThCrit./p_THSAT, p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN)

        # Instantiate
        new(p_THICK,p_STONEF,p_THSAT,p_Kθfc,p_KSAT,p_MvGα,p_MvGn,p_MvGm,p_MvGl,p_θr,
            NLAYER, p_PSIF, p_THETAF,p_PSIG,p_SWATMAX,p_WETF)
    end
end
KPT_SOILPAR_Mvg1d(;p_THICK::T, p_STONEF::T, p_THSAT::T, p_Kθfc::T, p_KSAT::T, p_MvGα::T, p_MvGn::T, p_MvGm::T, p_MvGl::T, p_θr::T) where {T<:AbstractVector} =
    KPT_SOILPAR_Mvg1d{T}(;p_THICK, p_STONEF, p_THSAT, p_Kθfc, p_KSAT, p_MvGα, p_MvGn, p_MvGm, p_MvGl, p_θr)
    # for explanation see https://docs.julialang.org/en/v1/manual/constructors/#Parametric-Constructors

# TODO: at some point in time define iterators or getindex for my type
    # function Base.getindex(p::KPT_SOILPAR_Mvg1d{Vector{Float64}}, idx::Int)
#     1 <= idx <= p.NLAYER || throw(BoundsError(p, idx))
#     p.p_THSAT[idx]
# end

"""
    derive_auxiliary_SOILVAR(u_SWATI, p_soil)

Derive alternative representations of soil water status.

Based on the state `u_SWATI` it returns (`u_aux_WETNES`, `u_aux_PSIM`, `u_aux_PSITI`, `p_fu_KK`)
- `u_aux_WETNES`: wetness, fraction of saturation
- `u_aux_PSIM`:   matric soil water potential for layer, kPa
- `p_PSIG`:       gravity potential, kPa
- `u_aux_PSITI`:  total potential ψt = ψm + ψg (sum of matrix potential and gravity potential)
- `u_aux_θ`:      volumetric soil water content, m3/m3
- `p_fu_KK`:      unsaturated hydraulic conductivity: K(Se) a.k.a. K(W)
"""
function derive_auxiliary_SOILVAR(u_SWATI,  p::KPT_SOILPAR_Mvg1d)

    u_aux_WETNES = (p.p_THSAT .* u_SWATI ./ p.p_SWATMAX .- p.p_θr) ./ (p.p_THSAT .- p.p_θr)
    u_aux_WETNES .= min.(1, u_aux_WETNES)

    u_aux_PSIM   = fill(NaN, length(p.p_SWATMAX)); FPSIM!( u_aux_PSIM, u_aux_WETNES, p)
    u_aux_θ      = fill(NaN, length(p.p_SWATMAX)); FTheta!(u_aux_θ,    u_aux_WETNES, p)
    p_fu_KK      = fill(NaN, length(p.p_SWATMAX)); FK_MvG!(p_fu_KK,    u_aux_WETNES, p.p_KSAT, p.p_MvGl, p.p_MvGn, p.p_MvGm)
    u_aux_PSITI  = u_aux_PSIM .+ p.p_PSIG

    return (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK)
end

function derive_auxiliary_SOILVAR!(u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK,
                                   u_SWATI,  p::KPT_SOILPAR_Mvg1d)

    u_aux_WETNES .= (p.p_THSAT .* u_SWATI ./ p.p_SWATMAX .- p.p_θr) ./ (p.p_THSAT .- p.p_θr)
    u_aux_WETNES .= min.(1, u_aux_WETNES)

    FPSIM!( u_aux_PSIM, u_aux_WETNES, p)
    FTheta!(u_aux_θ,    u_aux_WETNES, p)
    FK_MvG!(p_fu_KK,    u_aux_WETNES, p.p_KSAT, p.p_MvGl, p.p_MvGn, p.p_MvGm)

    u_aux_PSITI  .= u_aux_PSIM .+ p.p_PSIG

    return nothing
end

function derive_auxiliary_SOILVAR(u_SWATI,  p::KPT_SOILPAR_Ch1d)
    NLAYER   = size(p.p_KSAT,1)

    u_aux_WETNES = u_SWATI./p.p_SWATMAX
    u_aux_PSIM   = fill(NaN, length(p.p_SWATMAX)); FPSIM!( u_aux_PSIM, u_aux_WETNES, p)
    u_aux_θ      = fill(NaN, length(p.p_SWATMAX)); FTheta!(u_aux_θ,    u_aux_WETNES, p)

    u_aux_PSITI = fill(NaN, NLAYER)
    p_fu_KK     = fill(NaN, NLAYER)
    for i = 1:NLAYER
        u_aux_PSITI[i] = u_aux_PSIM[i] + p.p_PSIG[i]
        if u_aux_WETNES[i] > 0.00010
            p_fu_KK[i] = p.p_KF[i] * (u_aux_WETNES[i] / p.p_WETF[i]) ^ (2 * p.p_BEXP[i] + 3)
        else # extremely dry
            p_fu_KK[i] = 1E-10
        end
    end
    return (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK)
end

"""

    FK_MvG(WETNES, KSAT, MvGl, MvGn)

Compute hydraulic conductivity from wetness for the Mualem van Genuchten parametrization.

Compute unsaturated hydraulic conductivity: K(Se) a.k.a. K(W) using MvG equation 8)
K = Ks*W^l*[ 1 - (1-W^(1/m))^m ]^2 using m = 1-1/n yields: K = Ks*W^l*[ 1 - (1-W^(n/(n-1)))^(1-1/n) ]^2
"""
# function FK_MvG(WETNES, KSAT, MvGl, MvGn)
#     eps = 1.e-6
#     # return: hydraulic conductivity, mm/d
#     return KSAT .* max.(WETNES, eps) .^ MvGl .*
#         (1 .- (1 .- max.(WETNES, eps) .^ (MvGn ./ (MvGn .- 1)) ).^(1 .- 1 ./ MvGn) ) .^ 2
# end
function FK_MvG(WETNES, KSAT, MvGl, MvGn, MvGm)
    # return: hydraulic conductivity, mm/d
    eps = 1.e-6
    KSAT .* max.(WETNES, eps) .^ MvGl .*
        (1 .- (1 .- max.(WETNES, eps) .^ (1 ./ MvGm) ).^(MvGm) ) .^ 2
end
function FK_MvG!(result, WETNES, KSAT, MvGl, MvGn, MvGm)
    eps = 1.e-6
    # return: hydraulic conductivity, mm/d
    result .= KSAT .* max.(WETNES, eps) .^ MvGl .*
        (1 .- (1 .- max.(WETNES, eps) .^ (1 ./ MvGm) ).^(MvGm) ) .^ 2
end

"""
    FPSIM(u_aux_WETNES, p_soil)

Compute ψ(Se) = h(Se) a.k.a ψ(W) = h(W).
"""
function FPSIM!(ψM, u_aux_WETNES, p::KPT_SOILPAR_Ch1d)
    # FPSIM obtains Ψi from Wi for one layer using equation (7) in the clearly
    # unsaturated region and equation (4) in the near-saturation region.
    NLAYER = length(p.p_SWATMAX)
    for i in 1:NLAYER
        if u_aux_WETNES[i] <= 0.
            ψM[i] = -10000000000
        elseif u_aux_WETNES[i] < p.p_WETINF[i]
            # in clearly unsaturated range (eq. 7)
            ψM[i] = p.p_PSIF[i] * (u_aux_WETNES[i] / p.p_WETF[i]) ^ (-p.p_BEXP[i])
        elseif u_aux_WETNES[i] < 1.0
            # in near-saturated range (eq. 4)
            ψM[i] = p.p_CHM[i]* (u_aux_WETNES[i] - p.p_CHN[i]) * (u_aux_WETNES[i] - 1)
        else
            # saturated
            ψM[i] = 0.0
        end
    end
    return nothing
end
function FPSIM(u_aux_WETNES, p::KPT_SOILPAR_Ch1d)
    # FPSIM obtains Ψi from Wi for one layer using equation (7) in the clearly
    # unsaturated region and equation (4) in the near-saturation region.
    NLAYER = length(p.p_SWATMAX)
    ψM     = fill(NaN, NLAYER)
    for i in 1:NLAYER
        if u_aux_WETNES[i] <= 0.
            ψM[i] = -10000000000
        elseif u_aux_WETNES[i] < p.p_WETINF[i]
            # in clearly unsaturated range (eq. 7)
            ψM[i] = p.p_PSIF[i] * (u_aux_WETNES[i] / p.p_WETF[i]) ^ (-p.p_BEXP[i])
        elseif u_aux_WETNES[i] < 1.0
            # in near-saturated range (eq. 4)
            ψM[i] = p.p_CHM[i]* (u_aux_WETNES[i] - p.p_CHN[i]) * (u_aux_WETNES[i] - 1)
        else
            # saturated
            ψM[i] = 0.0
        end
    end
    return ψM
end
function FPSIM_CH(u_aux_WETNES,p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN)
    # needed only for construction of KPT_SOILPAR_Ch1d

    # FPSIM obtains Ψi from Wi for one layer using equation (7) in the clearly
    # unsaturated region and equation (4) in the near-saturation region.
    NLAYER = length(u_aux_WETNES)
    ψM     = fill(NaN, NLAYER)
    for i in 1:NLAYER
        if u_aux_WETNES[i] <= 0.
            ψM[i] = -10000000000
        elseif u_aux_WETNES[i] < p_WETINF[i]
            # in clearly unsaturated range (eq. 7)
            ψM[i] = p_PSIF[i] * (u_aux_WETNES[i] / p_WETF[i]) ^ (-p_BEXP[i])
        elseif u_aux_WETNES[i] < 1.0
            # in near-saturated range (eq. 4)
            ψM[i] = p_CHM[i]* (u_aux_WETNES[i] - p_CHN[i]) * (u_aux_WETNES[i] - 1)
        else
            # saturated
            ψM[i] = 0.0
        end
    end
    return ψM
end
function FPSIM!(result, u_aux_WETNES, p::KPT_SOILPAR_Mvg1d)
    FPSIM_MvG!(result, u_aux_WETNES, p.p_MvGα, p.p_MvGn, p.p_MvGm)
end
function FPSIM_MvG!(ψM, u_aux_WETNES, p_MvGα, p_MvGn, p_MvGm)
    eps = 1.e-6
    # MvGm = 1-1/MvGn
    # AWET = max.(u_aux_WETNES, eps)
    # ψM = 9.81 .* (-1 ./ p_MvGα) .* (AWET .^ (-1 ./ (1 .- 1 ./ p_MvGn)) .-1) .^ (1 ./ p_MvGn)
    # ψM = 9.81 .* (-1 ./ p_MvGα) .* (max.(u_aux_WETNES, eps) .^ (-1 ./ (1 .- 1 ./ p_MvGn)) .-1) .^ (1 ./ p_MvGn)
    ψM .= 9.81 .* (-1 ./ p_MvGα) .* (max.(u_aux_WETNES, eps) .^ (-1 ./ (p_MvGm)) .-1) .^ (1 ./ p_MvGn)

    # 9.81 conversion from m to kPa #TODO define and use const
end
function FPSIM(u_aux_WETNES, p::KPT_SOILPAR_Mvg1d)
    FPSIM_MvG(u_aux_WETNES, p.p_MvGα, p.p_MvGn, p.p_MvGm)
end
function FPSIM_MvG(u_aux_WETNES, p_MvGα, p_MvGn, p_MvGm)
    eps = 1.e-6
    # MvGm = 1-1/MvGn
    # AWET = max.(u_aux_WETNES, eps)
    # ψM = 9.81 .* (-1 ./ p_MvGα) .* (AWET .^ (-1 ./ (1 .- 1 ./ p_MvGn)) .-1) .^ (1 ./ p_MvGn)
    # ψM = 9.81 .* (-1 ./ p_MvGα) .* (max.(u_aux_WETNES, eps) .^ (-1 ./ (1 .- 1 ./ p_MvGn)) .-1) .^ (1 ./ p_MvGn)
    ψM = 9.81 .* (-1 ./ p_MvGα) .* (max.(u_aux_WETNES, eps) .^ (-1 ./ (p_MvGm)) .-1) .^ (1 ./ p_MvGn)

    # 9.81 conversion from m to kPa #TODO define and use const
end

"""
    FDPSIDWF!(dψδW, u_aux_WETNES, p_soil)

Compute derivative dψi/dWi for each layer i.

Ecoshift:
FDPSIDW returns dψi/dWi for one layer, which is needed for the selection of iteration
time-step.

For Clapp-Hornberger:
Differentiation of (7) and (4) leads to
dψi / dWi = ( -b ψf / Wf ) ( Wi / Wf )-b-1
in the unsaturated range,
dψi / dWi= m ( 2 Wi - n - 1 )
in the near saturation range, and
dψi / dWi = 0
when the soil is saturated (Wi = 1).

For Mualem-van Genuchten:
dψi / dWi = TODO.... write documentation here
"""
function FDPSIDWF!(dψδW, u_aux_WETNES, p::KPT_SOILPAR_Ch1d)
    # FDPSIDW returns dΨi/dWi for one layer, which is needed for the selection of iteration time-step.
    @unpack p_WETINF, p_BEXP, p_PSIF, p_WETF, p_CHM, p_CHN = p

    # dψδW  = zero(u_aux_WETNES)
    for i = 1:length(u_aux_WETNES)
        if u_aux_WETNES[i] < p_WETINF[i]
            # in clearly unsaturated range (eq. 7)
            dψδW[i] = (-p_BEXP[i] * p_PSIF[i] / p_WETF[i]) * (u_aux_WETNES[i] / p_WETF[i]) ^ (-p_BEXP[i] - 1)
            # TODO(bernhard): slightly faster:
            # dψδW[i] = p_PSIF[i]*p_WETF[i]^p_BEXP[i] * -p_BEXP[i]*u_aux_WETNES[i]^(-p_BEXP[i]-1)
        elseif u_aux_WETNES[i] < 1.0
            # in near-saturated range (eq. 4)
            dψδW[i] = p_CHM[i] * (2 * u_aux_WETNES[i] - p_CHN[i] - 1)
        else
            # saturated
            dψδW[i] = 0.0
        end
    end
    return dψδW # d PSI/d WETNES, kPa
end
function FDPSIDWF!(dψδW, u_aux_WETNES, p::KPT_SOILPAR_Mvg1d)
    α = p.p_MvGα
    n = p.p_MvGn

    eps = 1.e-6
    m = p.p_MvGm
    # dψδW = zeros(size(u_aux_WETNES))
    @inbounds for i = 1:length(u_aux_WETNES)
        # 9.81 conversion from m to kPa #TODO define and use const
        if (u_aux_WETNES[i] <= eps)
            dψδW[i] = 9.81 * (-1/α[i])*(1/n[i])*(            eps^(-1/(m[i]))-1)^(1/n[i]-1)*(-1/(m[i]))*            eps^(-1/(m[i])-1)
        end
        if ((u_aux_WETNES[i] > eps) && (u_aux_WETNES[i] < 1.0))
            dψδW[i] = 9.81 * (-1/α[i])*(1/n[i])*(u_aux_WETNES[i]^(-1/(m[i]))-1)^(1/n[i]-1)*(-1/(m[i]))*u_aux_WETNES[i]^(-1/(m[i])-1)
        end
        if (u_aux_WETNES[i] > 1.0)
            dψδW[i] = 0.0
        end
    end

    return dψδW # d PSI/d WETNES, kPa
end

"""
    FTheta(u_aux_WETNES, p_soil)

Compute θ based on Se.
"""
function FTheta(u_aux_WETNES, p::KPT_SOILPAR_Mvg1d)
    return FTheta_MvG(u_aux_WETNES, p.p_THSAT, p.p_θr)
end
function FTheta_MvG(u_aux_WETNES, p_THSAT, p_θr)
    # Computes θ(Se) = Se*(θs-θr) + θr
    return u_aux_WETNES .* (p_THSAT .- p_θr) .+ p_θr
end
function FTheta(u_aux_WETNES, p::KPT_SOILPAR_Ch1d)
    # Computes θ(Se) = Se*(θs-θr) + θr

    # variant 1:
    # p_θr = 0
    # return u_aux_WETNES*(p_THSAT-p_θr)+p_θr
    # variant 2:
    return u_aux_WETNES .* p.p_THSAT
end
function FTheta!(result, u_aux_WETNES, p::KPT_SOILPAR_Mvg1d)
    FTheta_MvG!(result, u_aux_WETNES, p.p_THSAT, p.p_θr)
end
function FTheta_MvG!(result, u_aux_WETNES, p_THSAT, p_θr)
    # Computes θ(Se) = Se*(θs-θr) + θr
    result .= u_aux_WETNES .* (p_THSAT .- p_θr) .+ p_θr
end
function FTheta!(result, u_aux_WETNES, p::KPT_SOILPAR_Ch1d)
    # Computes θ(Se) = Se*(θs-θr) + θr

    # variant 1:
    # p_θr = 0
    # u_aux_WETNES*(p_THSAT-p_θr)+p_θr
    # variant 2:
    result .= u_aux_WETNES .* p.p_THSAT
end
"""
    FWETNES(u_aux_PSIM, p_soil)

Computes θ(ψ) = θ(h) by computing first Se(ψ)=Se(h) a.k.a  W(ψ)=W(h)
"""
function FWETNES(u_aux_PSIM, p::KPT_SOILPAR_Mvg1d)
    NLAYER = length(u_aux_PSIM)
    WETNES = fill(NaN, NLAYER)
    for i = 1:NLAYER
        if u_aux_PSIM[i] <= 0
            # 1) W = (1 + (αh)^n )^(-m)
            α = p.p_MvGα[i]
            n = p.p_MvGn[i]
            m = p.p_MvGm[i]
            # MvGm = 1-1/MvGn
            WETNES[i] = (1+(-α*u_aux_PSIM[i]/9.81)^n)^(-(m)) # 9.81 conversion from kPa to m #TODO define and use const
        else
            WETNES[i] = 1.0
        end
    end
    return WETNES
end
function FWETNES(u_aux_PSIM, p::KPT_SOILPAR_Ch1d)
    error("FWETNES not implemented for Clapp+Hornberger. Use FSWATIPSIMF!")

    p_WETF = p.p_WETF
    p_WETINF = p.p_WETINF
    p_BEXP = p.p_BEXP
    p_PSIF = p.p_PSIF
    p_CHM = p.p_CHM
    p_CHN = p.p_CHN

    NLAYER = length(u_aux_PSIM)
    WETNES = fill(NaN, NLAYER)
    for i = 1:NLAYER
        if u_aux_PSIM[i] > 0
            error("STOP: Received positive ψ (u_aux_PSIM).")
        elseif (u_aux_PSIM[i] == 0.)
            WETNES[i] = 1.0
        else
            # in clearly unsaturated range (eq. 7)
            WETNES[i] = p_WETF[i] * (u_aux_PSIM[i] / p_PSIF[i]) ^ (-1 / p_BEXP[i])
            if (WETNES[i] > p_WETINF[i])
                # in near-saturated range (eq. 4)
                WETNES[i] = 0.5*( (1 + p_CHN[i]) + ((p_CHN[i]-1)^2              + 4*u_aux_PSIM[i]/p_CHM[i])^0.5)
                #WETNES[i] = 0.5*( (1 + p_CHN[i]) + (p_CHN[i]^2 - 2*p_CHN[i] + 1 + 4*u_aux_PSIM[i]/p_CHM[i])^0.5)
                #WETNES[i] = (1. + p_CHN[i])*0.5 + 0.5*(p_CHN[i]^2.0 - 2.0*p_CHN[i] + 1. + 4.0*u_aux_PSIM[i] / p_CHM[i])^(0.5)
            end
        end
    end
    return WETNES
end



"""
    SWCHEK(u_SWATI)

Correct u_SWATI and throw error if too far away.
"""
function SWCHEK!(u_SWATI, p_SWATMAX, t)
    # test for u_SWATI < 0 or u_SWATI > SWATMAX(I)

    # TODO(bernhard): instead of this manual correction, do it the DiffEq.jl way by using
    #                 a callback imposing the boundaries.
    #                 QUESTION: will the callback approach reduce the time step if boundaries are hit?
    #                           will this ensure zero mass balance error?

    for i in 1:length(u_SWATI)
        if u_SWATI[i] <= 0.0
            # error("You lose! u_SWATI = $(u_SWATI[i]) for layer $i at time $t."*
            # " Examine input and parameters to try to determine the cause")
            u_SWATI[i] = 0.0 # TODO(bernhard): remove correction and reactivate error

            # TODO(bernhard): This manual correction goes against the mass balance

        elseif u_SWATI[i] > p_SWATMAX[i]

            if u_SWATI[i] > p_SWATMAX[i] + 0.00001
                # TODO(bernhard): reactivate
                # error("u_SWATI > SWATMAX ($(u_SWATI[i]) > $(p_SWATMAX[i])) for layer $iat time $t."*
                # " Examine input and parameters to try to determine the cause")
                u_SWATI[i] = p_SWATMAX[i] # TODO(bernhard): remove correction and reactivate error
            else
                # rounding error only
                u_SWATI[i] = p_SWATMAX[i]
            end
        end
    end
end

end