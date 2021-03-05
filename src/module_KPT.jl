# fabian.bernhard@wsl.ch, 2021-02-03

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

export SOILPAR, derive_auxiliary_SOILVAR

using Roots: find_zero, Bisection # to find wetness for a given hydraulic conductivity


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


"""SOILPAR()\n Function that defines constant parameters from input_siteparam. TO BE REDEFINED
"""
function SOILPAR(p_RHOWG,
                 p_THICK, # layer thicknesses (mm)
                 p_THETAF,# θ at field capacity (-)
                 p_THSAT, # θ at saturation == matrix porosity (-)

                 p_STONEF,# stone volume fraction, unitless
                 p_BEXP,  # exponent for psi-theta relation
                 p_KF,    # hydraulic conductivity at field capacity (mm/d)
                 p_PSIF,  # ψ at field capacity (kPa)
                 p_WET∞,  # wetness at dry end of near-saturation range

                 p_Kθfc,
                 p_PSICR, # minimum plant leaf water potential (MPa)
                 p_Ksat, p_MvGl, p_MvGn, p_MvGα, p_θr,

                 NLAYER, IMODEL)

    if IMODEL == 0
        (p_ψg, p_SWATMX, p_WETfc, p_CHm, p_CHn, p_Ksat, p_PSIF, p_THETAF) =
            SOILPAR_CH(p_RHOWG,p_THICK,p_THETAF,p_THSAT,p_STONEF,p_BEXP,
                                          p_KF,p_PSIF,p_WET∞,p_Kθfc,p_PSICR,p_Ksat,
                                          p_MvGl, p_MvGn, p_MvGα, p_θr,
                                          NLAYER)
    elseif IMODEL == 1
        (p_ψg, p_SWATMX, p_WETfc, p_CHm, p_CHn, p_Ksat, p_PSIF, p_THETAF) =
            SOILPAR_MvG(p_RHOWG,p_THICK,p_THETAF,p_THSAT,p_STONEF,p_BEXP,
                                          p_KF,p_PSIF,p_WET∞,p_Kθfc,p_PSICR,p_Ksat,
                                          p_MvGl, p_MvGn, p_MvGα, p_θr,
                                          NLAYER)
    else
        error("Error in SOILPAR(), unexpected input IMODEL: $IMODEL. Valid values ar 0 or 1.")
    end

    return (p_ψg, p_SWATMX, p_WETfc, p_CHm, p_CHn, p_Ksat, p_PSIF, p_THETAF)
end

function SOILPAR_CH(p_RHOWG,p_THICK,p_THETAF,p_THSAT,p_STONEF,p_BEXP,
                                          p_KF,p_PSIF,p_WET∞,
                                          p_Kθfc,p_PSICR,p_Ksat,
                                          p_MvGl, p_MvGn, p_MvGα, p_θr,
                                          NLAYER)
    # Clapp and Hornberger:
    # Compute the quantities
    p_ψg     = fill(NaN, NLAYER) # gravity potential, kPa
    p_SWATMX = fill(NaN, NLAYER) # maximum water storage for layer, mm
    p_WETfc  = fill(NaN, NLAYER) # wetness at field capacity, dimensionless
    p_WETc   = fill(NaN, NLAYER) # wetness at PSICR, dimensionless
    p_CHm    = fill(NaN, NLAYER) # Clapp and Hornberger m, kPa
    p_CHn    = fill(NaN, NLAYER) # Clapp and Hornberger n, dimensionless
    p_Ksat   = fill(NaN, NLAYER) # saturated hydraulic conductivity, mm/d
    #u_aux_WETNES = fill(NaN, NLAYER) # wetness at current state u_aux_PSIM, dimensionelss
    #u_aux_SWATI  = fill(NaN, NLAYER) # water volume in layer, mm

    # using as inputs: p_RHOWG, p_THICK, p_THETAF, p_THSAT
    # using as inputs: p_THSAT, p_STONEF, p_BEXP,
    #                  p_THETAF, p_KF, p_PSIF
    #                  p_WET∞,

    # state variables:
    # u_aux_WETNES    # wetness, fraction of saturation
    # u_aux_SWATI     # water volume in layer, mm
    # u_aux_PSIM      # matric soil water potential for layer, kPa

    # Fill layer by layer
    for i = 1:NLAYER
        # compute p_ψg (gravity potential) as negative down from surface
        if i == 1
            p_ψg[1] = -p_RHOWG * p_THICK[1] / 2
        else
            p_ψg[i] = p_ψg[i - 1] - p_RHOWG * ((p_THICK[i - 1] + p_THICK[i]) / 2)
        end

        ### constant parameters:
        p_SWATMX[i] = p_THICK[i] * p_THSAT[i] * (1 - p_STONEF[i])
        p_WETfc[i]  = p_THETAF[i] / p_THSAT[i]
        p_Ksat[i]   = p_KF[i] * (1 / p_WETfc[i]) ^ (2 * p_BEXP[i] + 3)

        # compute p_CHm and p_CHn using intermediate p_ψ∞
        p_ψ∞ = p_PSIF[i] * (p_WET∞[i] / p_WETfc[i]) ^ -p_BEXP[i] # p_ψ∞: potential at dry end of near saturation range (kPa)
        p_CHm[i] = (-p_ψ∞ / (1 - p_WET∞[i]) ^ 2) - p_BEXP[i] * (-p_ψ∞) / (p_WET∞[i] * (1 - p_WET∞[i]))
        p_CHn[i] = 2 * p_WET∞[i] - 1 - (-p_ψ∞ * p_BEXP[i] / (p_CHm[i] * p_WET∞[i]))

        # NOTE(bernhard): removed:
        # NOTE(bernhard): in LWFBrook90R there was addtionally:
        #                  - test of validity of PSIM (assert that PSIM <= 0, i.e. error("matrix psi must be negative or zero"))
        #                  - computation of u_aux_WETNES:
        #                       if: PSIM==0 -> u_aux_WETNES = 1
        #                       else: u_aux_WETNES = p_WETfc[i] * (u0_aux_PSIM[i]/p_PSIF[i])^(1/BEXP)
        #                             if u_aux_WETNES>p_WET∞: u_aux_WETNES[i] = (1 + p_CHm[i])/2 + 0.5*sqrt(p_CHn[i]^2 - 2*p_CHn[i] + 1 + 4 * u0_aux_PSIM[i]/p_CHm[i])
        #                  - computation of u_aux_SWATI[i] = u_aux_WETNES[i] * p_SWATMX[i]
        #                  - computation of p_WETc:
        #                       p_WETc[i] = p_WETfc[i] * (1000 * p_PSICR / p_PSIF[i]) ^ (-1 / p_BEXP[i])
    end

    return (p_ψg, p_SWATMX, p_WETfc, p_CHm, p_CHn, p_Ksat, p_PSIF, p_THETAF)#, u_aux_WETNES, u_aux_SWATI)
end

function SOILPAR_MvG(p_RHOWG,p_THICK,p_THETAF,p_THSAT,p_STONEF,p_BEXP,
                                          p_KF,p_PSIF,p_WET∞,p_Kθfc,p_PSICR,p_Ksat,
                                          p_MvGl, p_MvGn, p_MvGα, p_θr,
                                          NLAYER)
    # For Mualem-van Genuchten
    # compute
    p_ψg     = fill(NaN, NLAYER) # gravity potential, kPa
    p_SWATMX = fill(NaN, NLAYER) # maximum water storage for layer, mm
    p_WETfc  = fill(NaN, NLAYER) # wetness at field capacity, dimensionless
    p_WETc   = fill(NaN, NLAYER) # wetness at PSICR, dimensionless
    p_PSIF       = fill(NaN, NLAYER) # matric potential at field capacity (kPa)
    p_THETAF     = fill(NaN, NLAYER) # θ at field capacity (-)
    # u_aux_WETNES = fill(NaN, NLAYER) # wetness at current state: u_aux_PSIM (dimensionelss)
    # u_aux_SWATI  = fill(NaN, NLAYER) # water volume in layer (mm)

    # using as inputs: p_RHOWG, p_THICK, p_THETAF, p_THSAT
    # using as inputs: p_Kθfc, p_THSAT, p_THICK, p_STONEF, u_aux_PSIM, p_PSICR, p_Ksat, p_MvGl, p_MvGn, p_MvGα

    # Fill layer by layer
    for i = 1:NLAYER
        # compute p_ψg (gravity potential) as negative down from surface
        if i == 1
            p_ψg[1] = -p_RHOWG * p_THICK[1] / 2
        else
            p_ψg[i] = p_ψg[i - 1] - p_RHOWG * ((p_THICK[i - 1] + p_THICK[i]) / 2)
        end

        ### constant parameters:

        # Find wetness at field capacity based on hydraulic conductivity
        # using a nonlinear solver for: find θ such that: K(θ) - p_Kθfc = 0
        # plot(-0.15:0.01:1.0, FK_MvG.(-0.15:0.01:1.0, p_Ksat[i], p_MvGl[i], p_MvGn[i]))
        f(x) = FK_MvG(x, p_Ksat[i], p_MvGl[i], p_MvGn[i]) - p_Kθfc[i];
        p_WETfc[i] = find_zero(f, (0.0, 1.0), Bisection())
        # p_WETfc[i] = FWETK(p_Kθfc[i], p_THSAT[i]) #FB: wetness at field capacity
        # if p_WETfc[i] == -99999.
        #     error("Computed invalid p_WETfc by FWETK()")
        # end

        p_SWATMX[i] = p_THICK[i] * p_THSAT[i] * (1.0 - p_STONEF[i])
        p_PSIF[i]   = FPSIM_MvG(p_WETfc[i], p_MvGα[i], p_MvGn[i]) #FB: matrix potential at field capacity, kPa
        p_THETAF[i] = FTheta_MvG(p_WETfc[i], p_THSAT[i], p_θr[i])  #FB: soil moisture at field capacity (mm)
        p_WETc[i]   = FWETNES_MvG(1000 * p_PSICR, p_MvGα[i], p_MvGn[i])
    end
    # for i = 1:NLAYER
    #     # state variables:
    #     # u_aux_WETNES(*) # wetness, fraction of saturation
    #     # u_aux_SWATI(*)  ! water volume in layer, mm
    #     # PSIM(*)   ! matric soil water potential for layer, kPa

    #     ### state dependent parameters:
    #     # NOTE(bernhard): removed:
    #     # NOTE(bernhard): in LWFBrook90R there was addtionally:
    #     # if u_aux_PSIM[i] > 0
    #     #     error("matrix psi must be negative or zero")
    #     # else
    #     #     u_aux_WETNES[i] = FWETNES_MvG(u_aux_PSIM[i], p_MvGα[i], p_MvGn[i])
    #     # end
    #     # u_aux_SWATI[i] = FTheta_MvG(u_aux_WETNES[i], p_THSAT[i], p_θr[i]) * p_SWATMX[i]/p_THSAT[i]
    # end
    p_CHm, p_CHn = (NaN, NaN)
    return (p_ψg, p_SWATMX, p_WETfc, p_CHm, p_CHn, p_Ksat, p_PSIF, p_THETAF)#, u_aux_WETNES, u_aux_SWATI)
end

""" derive_auxiliary_SOILVAR(u_SWATI, p_SWATMX, p_THSAT,
                 p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_KF,
                 p_θr, p_MvGα, p_MvGn, p_MvGl, p_Ksat,
                 p_PSIG, NLAYER, IMODEL)\n
Derives alternative representations of soil water status.
I.e. based on the state u_SWATI it returns (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, p_fu_KK)"""
function derive_auxiliary_SOILVAR(u_SWATI, p_SWATMX, p_THSAT,
                 p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_KF,
                 p_θr, p_MvGα, p_MvGn, p_MvGl, p_Ksat,
                 p_PSIG, NLAYER, IMODEL)
        if IMODEL == 0
            u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK =
            SOILVAR_CH(u_SWATI, p_SWATMX, p_THSAT,
                            p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_KF,
                            p_PSIG, NLAYER)
        else # IMODEL == 1
            u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK =
            SOILVAR_MvG(u_SWATI, p_SWATMX, p_THSAT,
                            p_θr, p_MvGα, p_MvGn, p_MvGl, p_Ksat,
                            p_PSIG, NLAYER)
        end
    return (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK)
end
function SOILVAR_MvG(u_SWATI,  p_SWATMX,
                                   p_THSAT, p_θr, p_MvGα, p_MvGn, p_MvGl, p_Ksat,
                                   p_PSIG, NLAYER)
    # Case where IMODEL == 1
    u_aux_WETNES = fill(NaN, NLAYER)
    u_aux_PSIM   = fill(NaN, NLAYER)
    u_aux_θ      = fill(NaN, NLAYER)
    u_aux_PSITI  = fill(NaN, NLAYER)
    p_fu_KK      = fill(NaN, NLAYER)

    for i = 1:NLAYER
        u_aux_WETNES[i] = (p_THSAT[i] * u_SWATI[i] / p_SWATMX[i] -p_θr[i]) / (p_THSAT[i] - p_θr[i])
        u_aux_WETNES[i] = min(1, u_aux_WETNES[i])

        u_aux_PSIM[i]   = FPSIM_MvG(u_aux_WETNES[i], p_MvGα[i], p_MvGn[i])
        u_aux_θ[i]      = FTheta_MvG(u_aux_WETNES[i], p_THSAT[i], p_θr[i])
        p_fu_KK[i]      = FK_MvG(u_aux_WETNES[i], p_Ksat[i], p_MvGl[i], p_MvGn[i])

        u_aux_PSITI[i]  = u_aux_PSIM[i] + p_PSIG[i]
    end

    # u_aux_WETNES: wetness, fraction of saturation
    # u_aux_PSIM:   matric soil water potential for layer, kPa
    # p_PSIG:       gravity potential, kPa
    # PSITI: total potential ψt = ψm + ψg (sum of matrix potential and gravity potential)
    # KK:    unsaturated hydraulic conductivity: K(Se) a.k.a. K(W)
    return (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK)
end
function SOILVAR_CH(u_SWATI,  p_SWATMX,
                                  p_THSAT, p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_KF,
                                  p_PSIG, NLAYER)
    # p_KF,         # hydraulic conductivity at field capacity, mm/d
    # p_WETF,       # wetness at field capacity, dimensionless
    # p_BEXP.       # exponent for psi-theta relation

    # Case where IMODEL == 0
    u_aux_WETNES = u_SWATI./p_SWATMX

    u_aux_PSIM = fill(NaN, NLAYER)
    u_aux_θ = fill(NaN, NLAYER)

    u_aux_PSITI=fill(NaN, NLAYER)
    p_fu_KK   =fill(NaN, NLAYER)

    for i = 1:NLAYER
        u_aux_PSIM[i] = FPSIMF_CH(u_aux_WETNES[i],
                                       p_PSIF[i], p_BEXP[i], p_WETINF[i], p_WETF[i], p_CHM[i], p_CHN[i])
        u_aux_θ[i]    = FTheta_CH(u_aux_WETNES[i], p_THSAT[i])

        u_aux_PSITI[i] = u_aux_PSIM[i] + p_PSIG[i]
        if u_aux_WETNES[i] > 0.00010
            p_fu_KK[i] = p_KF[i] * (u_aux_WETNES[i] / p_WETF[i]) ^ (2 * p_BEXP[i] + 3)
        else
            # extremely dry
            p_fu_KK[i] = 1E-10
        end
    end

    # u_aux_WETNES, # wetness, fraction of saturation
    # u_aux_PSIM,   # matric soil water potential for layer, kPa
    # p_PSIG,       # gravity potential, kPa
    return (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK)
end


# TODO(bernhard): remove remark that I renamed LWFBrook90_MvG_FK() to FK_MvG()
""" FK_MvG(WETNES, KSAT, MvGl, MvGn)\n Computes hydraulic conductivity from
wetness for the Mualem van Genuchten parametrization.\n\n

Computes unsaturated hydraulic conductivity: K(Se) a.k.a. K(W) using MvG equation 8)
K = Ks*W^l*[ 1 - (1-W^(1/m))^m ]^2 using m = 1-1/n yields: K = Ks*W^l*[ 1 - (1-W^(n/(n-1)))^(1-1/n) ]^2
"""
function FK_MvG(WETNES, KSAT, MvGl, MvGn)
    eps = 1.e-6
    AWET = max.(WETNES, eps)

    # return: hydraulic conductivity, mm/d
    return KSAT*AWET^MvGl*(1 - (1-AWET^(MvGn/(MvGn-1)) )^(1-1/MvGn) )^2
end

""" FPSIMF_CH(u_aux_WETNES,p_PSIF, p_BEXP, p_WET∞, p_WETF, p_CHM, p_CHN)\n Computes ψ(Se) = h(Se) a.k.a ψ(W) = h(W)
"""

""" FPSIM_MvG(u_aux_WETNES, p_MvGα, p_MvGn)\n Computes ψ(Se) = h(Se) a.k.a ψ(W) = h(W)
"""
#TODO(bernhard): check correct usage in SOILVAR() i.e. split in FPSIM_MvG() and FPSIMF_CH
# old definition: TODO: remove function LWFBrook90_MvG_FPSIMF(u_aux_WETNES,
# old definition: TODO: remove                                p_PSIF, p_BEXP, p_WET∞, p_WETF, p_CHM, p_CHN,
# old definition: TODO: remove                                p_MvGα, p_MvGn,
# old definition: TODO: remove                                iModel)
# old definition: TODO: replace by (for iModel = 0): FPSIMF_CH(u_aux_WETNES,p_PSIF, p_BEXP, p_WET∞, p_WETF, p_CHM, p_CHN)
# old definition: TODO: replace by (for iModel = 1): FPSIM_MvG(u_aux_WETNES, p_MvGα, p_MvGn)
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

# TODO(bernhard): remove remark that I renamed LWFBrook90_MvG_FTheta() to FTheta_MvG()
"""FTheta_MvG(u_aux_WETNES, p_θs, p_θr, iModel)
Computes θ based on Se.
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
""" FWETNES_MvG(u_aux_PSIM, p_MvGα, p_MvGn)\n
Computes θ(ψ) = θ(h) by computing first Se(ψ)=Se(h) a.k.a  W(ψ)=W(h)
"""

""" FWETNES_CH(u_aux_PSIM,p_WETF, p_WET∞, p_BEXP, p_PSIF, p_CHM, p_CHN)\n
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