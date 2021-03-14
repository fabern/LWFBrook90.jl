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

##   Clapp Hornberger (IMODEL=0) see: http://www.ecoshift.net/brook/kpt.html
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

##   Mualem van Genuchten (iModel=1)
1) W = 1/(1 + (αh)^n )^m                         (Radcliffe eq.2.47)
2) ψ = -1/α * [ (1-W^(1/m))*(W^(-1/m)) ]^(1/n)   (Radcliffe p.122)
...
8) K = Ks*W^l*[ 1 - (1-W^(1/m))^m ]^2

   ==> dψdW = .....
"

"""
module KPT # SOIL WATER PROPERTIES

export SOILPAR_CH, derive_auxiliary_SOILVAR
export KPT_SOILPAR_Mvg1d, KPT_SOILPAR_Ch1d
using Roots: find_zero, Bisection # to find wetness for a given hydraulic conductivity
using ..CONSTANTS: p_ThCrit,p_RHOWG # https://discourse.julialang.org/t/large-programs-structuring-modules-include-such-that-to-increase-performance-and-readability/29102/5

### ### Parameters
### #   Clapp and Hornberger (iModel=0)
### #   Par_model0 = [THSAT,THETAF,KF,PSIF,WETF,KSAT,CHM,CHN,BEXP,WETINF]
### # from LWFBrook90R: !   real(kind=8) :: PSIF     ! matrix potential at field capacity, kPa
### # from LWFBrook90R:     real(kind=8) :: BEXP     ! exponent for psi-theta relation
### # from LWFBrook90R: !   real(kind=8) :: WETINF   ! wetness at dry end of near-saturation range
### # from LWFBrook90R:     real(kind=8) :: WETF     ! saturation fraction at field capacity
### # from LWFBrook90R: !   real(kind=8) :: CHM      ! Clapp and Hornberger m, kPa
### # from LWFBrook90R: !   real(kind=8) :: CHN      ! Clapp and Hornberger n
### # from LWFBrook90R:     real(kind=8) :: KF       ! hydraulic conductivity at field capacity, mm/d
### #   Mualem van Genuchten (iModel=1)
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
abstract type AbstractKptSoilpar end
"""
Represents a discretized 1D column of soil with Clapp-Hornberger parametrization.

Input fields: p_THICK, p_STONEF, p_THSAT, p_PSIF, p_THETAF, p_KF, p_BEXP, p_WETINF
Derived fields: p_CHM, p_CHN, p_THETAF, p_PSIG, p_SWATMX, p_WETF, p_PsiCrit
"""
struct KPT_SOILPAR_Ch1d <: AbstractKptSoilpar
    # Input fields
    p_THICK::AbstractVector
    p_STONEF::AbstractVector
    p_THSAT::AbstractVector
    p_PSIF::AbstractVector
    p_THETAF::AbstractVector
    p_KF::AbstractVector
    p_BEXP::AbstractVector
    p_WETINF::AbstractVector
    # Derived fields
    p_CHM::AbstractVector
    p_CHN::AbstractVector
    p_PSIG::AbstractVector
    p_SWATMX::AbstractVector
    p_WETF::AbstractVector
    p_PsiCrit::AbstractVector
    # p_PsiCrit is the ψ value that corresponds to the constant, critical θ value p_ThCrit
    # Note that p_PSICR is different!


    # # Inner constructor:
    function KPT_SOILPAR_Ch1d(;p_THICK, p_STONEF, p_THSAT, p_PSIF, p_THETAF, p_KF, p_BEXP, p_WETINF)
        NLAYER = length(p_THICK)
        @assert size(p_THICK) == size(p_STONEF) == size(p_THSAT) == size(p_PSIF) ==
            size(p_THETAF) == size(p_KF) == size(p_BEXP) == size(p_WETINF)

        # Derive fields
        # Derive 2:
        # Variant 2a: used to be SOILPAR()
        # p_Kθfc = p_KSAT = p_MvGl = p_MvGn = p_MvGα = p_θr = fill(NaN, NLAYER)
        # (p_PSIG, p_SWATMX, p_WETF, p_CHM, p_CHN, p_KSAT, p_PSIF, p_THETAF ) =
        # SOILPAR_CH(p_RHOWG,
        #                             p_THICK,p_THETAF,p_THSAT,p_STONEF,p_BEXP,
        #                             p_KF,p_PSIF,p_WETINF,p_Kθfc,p_PSICR,p_KSAT,
        #                             p_MvGl, p_MvGn, p_MvGα, p_θr,
        #                             NLAYER)

        # Variant 2b:
        p_WETF   = p_THETAF ./ p_THSAT # wetness at field capacity, dimensionless
        p_KSAT   = p_KF .* (1 ./ p_WETF) .^ (2 .* p_BEXP .+ 3) # saturated hydraulic conductivity, mm/d
        p_SWATMX = p_THICK .* p_THSAT .* (1.0 .- p_STONEF) # maximum water storage for layer, mm
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

        # Derive 1:
        # # Variant 1a:
        # p_PsiCrit = FPSIMF_CH.(p_ThCrit./p_THSAT,
        #                        p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN)
        # Variant 1b:
        WetnesCrit = p_ThCrit ./ p_THSAT
        # Obtain Ψi from Wi for one layer using equation (7) in the clearly
        # unsaturated region and equation (4) in the near-saturation region.
        p_PsiCrit = fill(NaN, NLAYER)
        for i = 1:NLAYER
            if WetnesCrit[i] <= 0.
                p_PsiCrit[i] = -10000000000
            elseif WetnesCrit[i] < p_WETINF[i]
                # in clearly unsaturated range (eq. 7)
                p_PsiCrit[i] = p_PSIF[i] * (WetnesCrit[i] / p_WETF[i]) ^ (-p_BEXP[i])
            elseif WetnesCrit[i] < 1.0
                # in near-saturated range (eq. 4)
                p_PsiCrit[i] = p_CHM[i] * (WetnesCrit[i] - p_CHN[i]) * (WetnesCrit[i] - 1)
            else
                # saturated
                p_PsiCrit[i] = 0.0
            end
        end

        # NOTE(bernhard): removed further tasks that were in SOILPAR() in LWFBrook90R
        #                  - computation of p_WETc = p_WETF .* (1000 .* p_PSICR ./ p_PSIF) .^ (-1 ./ p_BEXP) # wetness at p_PSICR, dimensionless
        #                  - test of validity of PSIM (assert that PSIM <= 0, i.e. error("matrix psi must be negative or zero"))
        #                  - computation of u_aux_WETNES:
        #                       if: PSIM==0 -> u_aux_WETNES = 1
        #                       else: u_aux_WETNES = p_WETF .* (u0_aux_PSIM ./ p_PSIF).^(1/.BEXP)
        #                             if u_aux_WETNES[i]>p_WETINF[i]: u_aux_WETNES[i] = (1 + p_CHM[i])/2 + 0.5*sqrt(p_CHN[i]^2 - 2*p_CHN[i] + 1 + 4 * u0_aux_PSIM[i]/p_CHM[i])
        #                  - computation of u_aux_SWATI = u_aux_WETNES .* p_SWATMX

        # Instantiate
        new(p_THICK, p_STONEF, p_THSAT, p_PSIF, p_THETAF, p_KF, p_BEXP, p_WETINF,
            p_CHM, p_CHN, p_PSIG, p_SWATMX, p_WETF, p_PsiCrit)
    end
end

"""
Represents a discretized 1D column of soil with Mualem-van Genuchten parametrization.

Input fields: p_THICK, p_STONEF, p_THSAT, p_Kθfc, p_KSAT, p_MvGα, p_MvGn, p_MvGl, p_θr
Derived fields: p_PSIF, p_THETAF, p_PSIG, p_SWATMX, p_WETF, p_PsiCrit
"""
struct KPT_SOILPAR_Mvg1d <: AbstractKptSoilpar
    # Input fields
    p_THICK::AbstractVector
    p_STONEF::AbstractVector
    p_THSAT::AbstractVector
    p_Kθfc::AbstractVector
    p_KSAT::AbstractVector
    p_MvGα::AbstractVector
    p_MvGn::AbstractVector
    p_MvGl::AbstractVector
    p_θr::AbstractVector
    # Derived fields
    p_PSIF::AbstractVector    # matric potential at field capacity, kPa
    p_THETAF::AbstractVector  # soil moisture θ at field capacity, m3/m3
    p_PSIG::AbstractVector    # gravity potential negative down from surface, kPa
    p_SWATMX::AbstractVector  # maximum water storage for layer, mm
    p_WETF::AbstractVector    # wetness at field capacity, dimensionless
    p_PsiCrit::AbstractVector
    # p_PsiCrit is the ψ value that corresponds to the constant, critical θ value p_ThCrit
    # Note that p_PSICR is different!


    # # Inner constructor:
    function KPT_SOILPAR_Mvg1d(;p_THICK, p_STONEF, p_THSAT, p_Kθfc, p_KSAT, p_MvGα, p_MvGn, p_MvGl, p_θr)
        NLAYER = length(p_THICK)
        @assert size(p_THICK) == size(p_STONEF) == size(p_THSAT) == size(p_Kθfc) ==
                size(p_KSAT) == size(p_MvGα) == size(p_MvGn) == size(p_MvGl) == size(p_θr)

        # Derive fields
        # Derive 1:
        # Variant 1a:
        # p_PsiCrit =
        #     FPSIM_MvG.(p_ThCrit./(p_THSAT .- p_θr),
        #                                 p_MvGα, p_MvGn)
        # Variant 1b:
        eps       = 1.e-6
        # MvGm      = 1-1/MvGn
        AWET      = max.(eps, p_ThCrit./(p_THSAT .- p_θr))
        p_PsiCrit = 9.81 * (-1 ./ p_MvGα).*(AWET.^(-1 ./ (1 .- 1 ./ p_MvGn)).-1).^(1 ./ p_MvGn) # 9.81 conversion from m to kPa

        # Derive 2:
        # # Variant 2a: used to be SOILPAR()
        # p_THETAF = p_PSIF = p_BEXP = p_KF = p_WETINF = fill(NaN, NLAYER)
        # (p_PSIG, p_SWATMX, p_WETF, p_KSAT, p_PSIF, p_THETAF) =
        #     SOILPAR_MvG(p_RHOWG,
        #                                 p_THICK,p_THETAF,p_THSAT,p_STONEF,p_BEXP,
        #                                 p_KF,p_PSIF,p_WETINF,p_Kθfc,p_PSICR, p_KSAT,
        #                                 p_MvGl, p_MvGn, p_MvGα, p_θr,
        #                                 NLAYER)

        # Variant 2b:
        p_SWATMX = p_THICK .* p_THSAT .* (1.0 .- p_STONEF) # maximum water storage for layer, mm
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
            # plot(-0.15:0.01:1.0, FK_MvG.(-0.15:0.01:1.0, p_KSAT[i], p_MvGl[i], p_MvGn[i]))
            f(x) = - p_Kθfc[i] + FK_MvG(x, p_KSAT[i], p_MvGl[i], p_MvGn[i])
            p_WETF[i] = find_zero(f, (0.0, 1.0), Bisection())
            # In LWFBrook90R: this was done using FWETK()
            # In LWFBrook90R: if p_WETF[i] == -99999. error("Computed invalid p_WETF by FWETK()") end
        end
        p_PSIF   = FPSIM_MvG.(p_WETF, p_MvGα, p_MvGn)  # matric potential at field capacity, kPa
        p_THETAF = FTheta_MvG.(p_WETF, p_THSAT, p_θr)  # soil moisture θ at field capacity, m3/m3

        # NOTE(bernhard): removed further tasks that were in SOILPAR() in LWFBrook90R
        #                  - computation of p_WETc = FWETNES_MvG.(1000 .* p_PSICR, p_MvGα, p_MvGn) # wetness at p_PSICR, dimensionless
        #                  - test of validity of PSIM (assert that PSIM <= 0, i.e. error("matrix psi must be negative or zero"))
        #                  - derivation of state dependent parameters u_aux_WETNES and u_aux_SWATI
        #                       - u_aux_WETNES[i] = FWETNES_MvG(u_aux_PSIM[i], p_MvGα[i], p_MvGn[i])
        #                       - u_aux_SWATI[i] = FTheta_MvG(u_aux_WETNES[i], p_THSAT[i], p_θr[i]) * p_SWATMX[i]/p_THSAT[i]

        # Instantiate
        new(p_THICK,p_STONEF,p_THSAT,p_Kθfc,p_KSAT,p_MvGα,p_MvGn,p_MvGl,p_θr,
            p_PSIF, p_THETAF,p_PSIG,p_SWATMX,p_WETF,p_PsiCrit)
    end
end


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
    NLAYER   = size(p.p_KSAT,1)

    u_aux_WETNES = fill(NaN, NLAYER)
    u_aux_PSIM   = fill(NaN, NLAYER)
    u_aux_θ      = fill(NaN, NLAYER)
    u_aux_PSITI  = fill(NaN, NLAYER)
    p_fu_KK      = fill(NaN, NLAYER)
    for i = 1:NLAYER
        u_aux_WETNES[i] = (p.p_THSAT[i] * u_SWATI[i] / p.p_SWATMX[i] -p.p_θr[i]) / (p.p_THSAT[i] - p.p_θr[i])
        u_aux_WETNES[i] = min(1, u_aux_WETNES[i])

        u_aux_PSIM[i]   = FPSIM_MvG(u_aux_WETNES[i], p.p_MvGα[i], p.p_MvGn[i])
        u_aux_θ[i]      = FTheta_MvG(u_aux_WETNES[i], p.p_THSAT[i], p.p_θr[i])
        p_fu_KK[i]      = FK_MvG(u_aux_WETNES[i], p.p_KSAT[i], p.p_MvGl[i], p.p_MvGn[i])

        u_aux_PSITI[i]  = u_aux_PSIM[i] + p.p_PSIG[i]
    end
    return (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK)
end
function derive_auxiliary_SOILVAR(u_SWATI,  p::KPT_SOILPAR_Ch1d)
    NLAYER   = size(p.p_KSAT,1)

    u_aux_WETNES = u_SWATI./p.p_SWATMX

    u_aux_PSIM  = fill(NaN, NLAYER)
    u_aux_θ     = fill(NaN, NLAYER)
    u_aux_PSITI = fill(NaN, NLAYER)
    p_fu_KK     = fill(NaN, NLAYER)

    for i = 1:NLAYER
        u_aux_PSIM[i]  = FPSIMF_CH(
            u_aux_WETNES[i],p.p_PSIF[i], p.p_BEXP[i], p.p_WETINF[i], p.p_WETF[i], p.p_CHM[i], p.p_CHN[i])
        u_aux_θ[i]     = FTheta_CH(u_aux_WETNES[i], p.p_THSAT[i])
        u_aux_PSITI[i] = u_aux_PSIM[i] + p.p_PSIG[i]
        if u_aux_WETNES[i] > 0.00010
            p_fu_KK[i] = p.p_KF[i] * (u_aux_WETNES[i] / p.p_WETF[i]) ^ (2 * p.p_BEXP[i] + 3)
        else # extremely dry
            p_fu_KK[i] = 1E-10
        end
    end
    return (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK)
end


# TODO(bernhard): remove remark that I renamed LWFBrook90_MvG_FK() to FK_MvG()
"""

    FK_MvG(WETNES, KSAT, MvGl, MvGn)

Compute hydraulic conductivity from wetness for the Mualem van Genuchten parametrization.

Compute unsaturated hydraulic conductivity: K(Se) a.k.a. K(W) using MvG equation 8)
K = Ks*W^l*[ 1 - (1-W^(1/m))^m ]^2 using m = 1-1/n yields: K = Ks*W^l*[ 1 - (1-W^(n/(n-1)))^(1-1/n) ]^2
"""
function FK_MvG(WETNES, KSAT, MvGl, MvGn)
    eps = 1.e-6
    AWET = max.(WETNES, eps)

    # return: hydraulic conductivity, mm/d
    return KSAT*AWET^MvGl*(1 - (1-AWET^(MvGn/(MvGn-1)) )^(1-1/MvGn) )^2
end

"""
    FPSIMF_CH(u_aux_WETNES,p_PSIF, p_BEXP, p_WET∞, p_WETF, p_CHM, p_CHN)

Compute ψ(Se) = h(Se) a.k.a ψ(W) = h(W).
"""

"""
    FPSIM_MvG(u_aux_WETNES, p_MvGα, p_MvGn)

Compute ψ(Se) = h(Se) a.k.a ψ(W) = h(W).
"""
function FPSIMF_CH(u_aux_WETNES,p_PSIF, p_BEXP, p_WET∞, p_WETF, p_CHM, p_CHN)
    # Computes ψ(Se) = h(Se) a.k.a ψ(W) = h(W)
    #
    # iModel == 0 (Clapp + Hornberger)
    # FPSIM obtains Ψi from Wi for one layer using equation (7) in the clearly
    # unsaturated region and equation (4) in the near-saturation region.
    if u_aux_WETNES <= 0.
        ψM = -10000000000
    elseif u_aux_WETNES < p_WET∞
        # in clearly unsaturated range (eq. 7)
        ψM = p_PSIF * (u_aux_WETNES / p_WETF) ^ (-p_BEXP)
    elseif u_aux_WETNES < 1.0
        # in near-saturated range (eq. 4)
        ψM = p_CHM* (u_aux_WETNES - p_CHN) * (u_aux_WETNES - 1)
    else
        # saturated
        ψM = 0.0
    end
    return ψM
end
function FPSIM_MvG(u_aux_WETNES, p_MvGα, p_MvGn)
    # iModel == 1 (Mualem - van Genuchten)
    # ...
    eps = 1.e-6
    # MvGm = 1-1/MvGn
    AWET = max.(u_aux_WETNES, eps)
    ψM = (-1/p_MvGα)*(AWET^(-1/(1-1/p_MvGn))-1)^(1/p_MvGn)
    ψM = ψM * 9.81 # 9.81 conversion from m to kPa #TODO define and use const

    return ψM
end


"""
Ecoshift:
FDPSIDW returns dψi/dWi for one layer, which is needed for the selection of iteration
time-step. Differentiation of (7) and (4) leads to
dψi / dWi = ( -b ψf / Wf ) ( Wi / Wf )-b-1
in the unsaturated range,
dψi / dWi= m ( 2 Wi - n - 1 )
in the near saturation range, and
dψi / dWi = 0
when the soil is saturated (Wi = 1).
"""

function FDPSIDWF_CH(u_aux_WETNES, p_WET∞, p_BEXP, p_PSIF, p_WETF, p_CHM, p_CHN)
    # FDPSIDW returns dΨi/dWi for one layer, which is needed for the selection of iteration time-step.

    dψδW = zeros(size(u_aux_WETNES))
    for i = 1:length(u_aux_WETNES)
        if u_aux_WETNES[i] < p_WET∞[i]
            # in clearly unsaturated range (eq. 7)
            dψδW[i] = (-p_BEXP[i] * p_PSIF[i] / p_WETF[i]) * (u_aux_WETNES[i] / p_WETF[i]) ^ (-p_BEXP[i] - 1)
            # TODO(bernhard): slightly faster: dψδW[i] = p_PSIF[i]*p_WETF[i]^p_BEXP[i] * -p_BEXP[i]*u_aux_WETNES[i]^(-p_BEXP[i]-1)
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

function FDPSIDWF_MvG(u_aux_WETNES, α, n)
    eps = 1.e-6
    m = 1 .- 1 ./ n
    dψδW = zeros(size(u_aux_WETNES))
    for i = 1:length(u_aux_WETNES)
        if (u_aux_WETNES[i] <= eps)
            dψδW[i] = (-1/α[i])*(1/n[i])*(            eps^(-1/(m[i]))-1)^(1/n[i]-1)*(-1/(m[i]))*            eps^(-1/(m[i])-1)
        end
        if ((u_aux_WETNES[i] > eps) && (u_aux_WETNES[i] < 1.0))
            dψδW[i] = (-1/α[i])*(1/n[i])*(u_aux_WETNES[i]^(-1/(m[i]))-1)^(1/n[i]-1)*(-1/(m[i]))*u_aux_WETNES[i]^(-1/(m[i])-1)
        end
        if (u_aux_WETNES[i] > 1.0)
            dψδW[i] = 0.0
        end
    end

    return dψδW = dψδW * 9.81 # # 9.81 conversion from m to kPa #TODO define and use const
                              # d PSI/d WETNES, kPa
end

"""
    FTheta_MvG(u_aux_WETNES, p_θs, p_θr, iModel)

Compute θ based on Se.
"""
function FTheta_MvG(u_aux_WETNES, p_θs, p_θr)
    # Computes θ(Se) = Se*(θs-θr) + θr
    return u_aux_WETNES*(p_θs-p_θr)+p_θr
end
function FTheta_CH(u_aux_WETNES, p_θs)
    # Computes θ(Se) = Se*(θs-θr) + θr

    # variant 1:
    # p_θr = 0
    # return u_aux_WETNES*(p_θs-p_θr)+p_θr
    # variant 2:
    return u_aux_WETNES*p_θs
end

# TODO(bernhard): remove remark that I renamed LWFBrook90_MvG_FWETNES() to FWETNES_MvG()
"""
    FWETNES_MvG(u_aux_PSIM, p_MvGα, p_MvGn)

Computes θ(ψ) = θ(h) by computing first Se(ψ)=Se(h) a.k.a  W(ψ)=W(h)
"""
function FWETNES_MvG(u_aux_PSIM, p_MvGα, p_MvGn)
    # Computes θ(ψ) = θ(h) by computing first Se(ψ)=Se(h) a.k.a  W(ψ)=W(h)
    #
    # FWETNES obtains wetness from matrix potential PSIM

    # WETNEs = NaN # fill(NaN, length(u_aux_PSIM)) # NOTE(bernhard): works only for a single argument
    #for i = 1:length(u_aux_PSIM)
    if u_aux_PSIM <= 0
        # 1) W = (1 + (αh)^n )^(-m)
        # MvGm = 1-1/MvGn
        WETNEs = (1+(-p_MvGα*u_aux_PSIM/9.81)^p_MvGn)^(-(1-1/p_MvGn)) # 9.81 conversion from kPa to m #TODO define and use const
    else
        WETNEs = 1.0
    end
    #end

    return WETNEs
end

"""
    FWETNES_CH(u_aux_PSIM,p_WETF, p_WET∞, p_BEXP, p_PSIF, p_CHM, p_CHN)

Computes θ(ψ) = θ(h) by computing first Se(ψ)=Se(h) a.k.a  W(ψ)=W(h)
"""
function FWETNES_CH(u_aux_PSIM,p_WETF, p_WET∞, p_BEXP, p_PSIF, p_CHM, p_CHN)
    # Computes θ(ψ) = θ(h) by computing first Se(ψ)=Se(h) a.k.a  W(ψ)=W(h)
    #
    # FWETNES obtains wetness from matrix potential PSIM

    error("FWETNES not implemented for Clapp+Hornberger. Use FSWATIPSIMF!")

    # WETNEs = NaN # fill(NaN, length(u_aux_PSIM)) # NOTE(bernhard): works only for a single argument
    #for i = 1:length(u_aux_PSIM)
    if u_aux_PSIM > 0
        error("STOP: Received positive ψ (u_aux_PSIM).")
        # TODO(bernhard): or should it set WETNEs = u_aux_WETNES
    elseif (u_aux_PSIM == 0.)
        WETNEs = 1.0
    else
        # in clearly unsaturated range (eq. 7)
        WETNEs= p_WETF * (u_aux_PSIM / p_PSIF) ^ (-1 / p_BEXP)
        if (WETNEs > p_WET∞)
            # in near-saturated range (eq. 4)
            WETNEs = 0.5*( (1 + p_CHN) + ((p_CHN-1)^2              + 4*u_aux_PSIM/p_CHM)^0.5)
            #WETNEs[i] = 0.5*( (1 + p_CHN[i]) + (p_CHN[i]^2 - 2*p_CHN[i] + 1 + 4*u_aux_PSIM[i]/p_CHM[i])^0.5)
            #WETNEs[i] = (1. + p_CHN[i])*0.5 + 0.5*(p_CHN[i]^2.0 - 2.0*p_CHN[i] + 1. + 4.0*u_aux_PSIM[i] / p_CHM[i])^(0.5)
        end
    end
    #end
    return WETNEs
end

end