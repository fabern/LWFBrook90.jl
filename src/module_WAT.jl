# fabian.bernhard@wsl.ch, 2021-02-03

module WAT # WATER MOVEMENT IN SOIL

export INFPAR, LWFRootGrowth, ITER

# Ecoshift-WAT:
# BROOK90 does not try to account for all possible paths of water movement through soil,
# which is a very complicated subject (McDonnell 1990). As a lumped parameter model it does
# not move water laterally between subareas. However, it does simulate several different
# pathways of water movement over or through the soil to streamflow and groundwater (see
# Flow Chart):

# snowmelt and throughfall on impervious areas goes directly to streamflow (SRFL) snowmelt
# and precipitation on variable saturated source areas goes directly to streamflow (SRFL)
# remaining snowmelt and streamflow (SLFL) infiltrates either all to the surface layer or to
# several layers via vertical pipes or macropores (INFL) some infiltrated water can go
# directly to streamflow as "bypass" flow in pipes (BYFL) classic vertical matric flow
# between soil layers (VRFLI) lateral or downslope movement of matric water to streamflow
# (DSFL) vertical drainage of matric water to groundwater (VRFLN) discharge of groundwater
# to streamflow (GWFL) deep seepage loss from groundwater (SEEP)

# SLFL, INFL, VRFLI, and VRFLN are internal flows. SRFL and BYFL produce streamflow only on
# the day of precipitation, simulating a streamflow response of less than 24 hours duration.
# Users generally should choose one or the other of these flows. DSFL produces a response
# over several days only if VRFL is limited by a low conductivity layer within the profile
# or VRFLN is limited by setting DRAIN < < 1. VRFL from the bottom of the profile produces a
# response of several days only if there is no groundwater. GWFL reponse can vary from
# several to many days. SEEP produces no streamflow at all.

# If the surface horizon becomes saturated, infiltration-excess overland flow is simulated
# as BYFL from the top layer. The input rate (SLFL) is constant over the precipitation
# interval, which may be a full day for most users. So with saturated hydraulic
# conductivities usually 200 mm/d or more, such overland flow will be rare.

# Vertical water movement out of a layer is a combination of matrix flow, VRFL(I), and the
# macropore infiltration, SLFL(I). SLFL(I) either becomes BYFL from deeper layers or becomes
# soil water in deeper layers. VRFL(I) will generally increase with depth as SLFL(I)
# decreases. Lysimeter users may want to assume that suction lysimeters collect VRFL(I)
# whereas zero-tension lysimeters collect both VRFL(I) and SLFLI(I).

# In this section a subscript i is used to indicate individual layers.

"""INFPAR(p_INFEXP, p_ILAYER, p_THICK, NLAYER)\n Computes fraction of infiltration to each
soil layer.\n\n
Arguments:
- p_INFEXP: infiltration exponent: 0 all to top, 1 uniform with depth, >1.0=more at bottom than at top
- p_ILAYER: number of layers over which infiltration is distributed
- p_THICK
- NLAYER
"""
function INFPAR(p_INFEXP, p_ILAYER, p_THICK, NLAYER)
    p_INFRAC = fill(NaN, NLAYER) # fraction of infiltration to each layer
    if p_INFEXP <= 0
        p_INFRAC[1]         = 1
        p_INFRAC[2:NLAYER] .= 0
    else
        THICKT =    sum(p_THICK[1:p_ILAYER])
        THICKA = cumsum(p_THICK[1:p_ILAYER])

        p_INFRAC[1]         = (THICKA[1]/THICKT)^p_INFEXP - (0/THICKT)^p_INFEXP
        for i=2:NLAYER
            if (i <= p_ILAYER)
                p_INFRAC[i] = (THICKA[i]/THICKT)^p_INFEXP - (THICKA[i-1]/THICKT)^p_INFEXP
            else
                p_INFRAC[i] = 0
            end
        end
    end

    return p_INFRAC
end




"""LWFRootGrowth(frelden, tini, age, rgroper, inirdep, inirlen, NLAYER)\n Computes root growth
according to LWF root growth model, (Hammel and Kennel 2000).\n\n
Arguments:
    - frelden[] :  final relative values of root length per unit volume
    - tini[]    :  initial time for root growth in layer
    - age       :  age of vegetation
    - rgroper   :  period of root growth in layer, a
    - inirdep   :  intial root depth, m
    - inirlen   :  intial total root length, m m-2
    - NLAYER    :  number of soil layers
Returns:
    - RELDEN[]  : current, age-dependent relative values of root length per unit volume
"""
function LWFRootGrowth(frelden, tini, age, rgroper, inirdep, inirlen, NLAYER)

    p_fT_RELDEN = fill(NaN, NLAYER)
    if rgroper > zero(rgroper)
        for i = 1:NLAYER
            if age < tini[i]
                p_fT_RELDEN[i]=0.0
            elseif age >= tini[i] && age <= tini[i] + rgroper
                # rl0:  constant intial root length density, m m-3
                rl0=inirlen/inirdep
                p_fT_RELDEN[i]=rl0*(frelden[i]/rl0)^((age-tini[i])/rgroper)
            elseif age > tini[i] + rgroper
                p_fT_RELDEN[i]=frelden[i]
            else
                error("In RootGrowth() unexpected error occurred.")
            end
        end
    else
        for i = 1:NLAYER
            p_fT_RELDEN[i]=frelden[i]
        end
    end

    return p_fT_RELDEN
end


"""
BYFLFR(NLAYER, p_BYPAR, p_QFPAR, p_QFFC, u_aux_WETNES, p_WETF)\n
TODO(bernhard):

Ecoshift:
Bypass flow (BYFL) and surface or source area flow (SRFL) are the two stormflow or quickflow
generating mechanisms in BROOK90. The conceptual difference is that SRFL is "new" water that
has not infiltrated but has moved across the surface to a channel, whereas BYFL is "new"
water that has moved to a channel below the surface via macropores or pipes. In BROOK90 the
amount of BYFL from each layer depends on the wetness of that particular layer and on the
amount of infiltration to it, which is controlled by INFEXP. The amount of SRFL depends on
the total wetness of all soil layers down to and including input parameter QDEPTH . In
general users should not try to model BYFL and SRFL simultaneously because trying to fit
parameters for both at the same time would be too complicated. The parameter BYPAR is set to
1 to allow BYFL and zero to prevent it. SRFL is prevented by setting both QDEPTH and IMPERV
to zero. The same parameters, QFFC and QFPAR, are used for SRFL and for BYFL from all
layers.

When BYPAR = 0, there is no bypass flow from deeper layers, but bypass flow is still
generated from layer 1 when the layer would otherwise become oversaturated. When INFEXP = 0
or IDEPTH = 0, BYFL can only be generated from layer 1 because it is the only layer
receiving infiltrated water.

In each iteration loop, subroutine BYFLFR calculates the fraction of the water infiltrating
to each layer that becomes bypass flow (BYFRACi) as

BYFRACi = QFFC ^ [ 1 - (1 / QFPAR) * (WETNESi - WETFi) / (1 - WETFi ) ]

where WETNESi is the layer wetness, WETFi is the layer wetness at field capacity, and QFFC
and QFPAR are parameters. QFFC is the bypass fraction at field capacity (WETNESi = WETFi).
QFPAR represents the fraction of the water content between field capacity and saturation at
which BYFRAC reaches 1 (Fig. WAT-2). BYFRAC increases exponentially with layer wetness and
is prevented from exceeding 1. Raising QFFC will raise BYFL proportionally at all water
contents. Raising QFPAR will increase BYFL for soil drier than field capacity and decrease
it for soil above field capacity (Fig. WAT-2).

With BYPAR = 1 and QFPAR = 0, BYFRAC is 0 below field capacity and 1 above it, and QFFC is
ignored. With a single soil layer and no SRFL, this produces a classic "bucket" model, with
the leakiness of the bucket determined by DRAIN. With multiple layers and BYPAR = 1, the
bucket model does not work because a layer at field capacity diverts all excess water to
BYFL, preventing wetting of deeper layers. See SRFLFR for a modified bucket model with
multiple layers.

When QFPAR = 1, BYFRAC reaches 1 when the layer is saturated (Fig WAT-2). A very large QFPAR
produces a constant BYFRAC of QFFC.

Note that BYFRAC is calculated from soil water prior to the input of water for the time
step.
"""
function BYFLFR(NLAYER, p_BYPAR, p_QFPAR, p_QFFC, u_aux_WETNES, p_WETF)
    # TODO(bernhard): could be optimized by not allocating each time new memory (versus in-place)
    p_fu_BYFRAC = fill(NaN, NLAYER)
    for i = 1:NLAYER
        if (isone(p_BYPAR))
            ####
            # variant Brook90: if (p_QFPAR > 0.01)
            # variant Brook90:     p_fu_BYFRAC[i] = p_QFFC ^ (1 - (1 / p_QFPAR) * (u_aux_WETNES[i] - p_WETF[i]) / (1 - p_WETF[i]))
            # variant Brook90:     if (p_fu_BYFRAC[i] > 1)
            # variant Brook90:         p_fu_BYFRAC[i] = 1
            # variant Brook90:     end
            # variant Brook90: else
            # variant Brook90:     # bucket for the layer
            # variant Brook90:     if (u_aux_WETNES[i] >= p_WETF[i])
            # variant Brook90:         p_fu_BYFRAC[i] = 1
            # variant Brook90:     else
            # variant Brook90:         p_fu_BYFRAC[i] = 0
            # variant Brook90:     end
            # variant Brook90: end
            ####

            ####
            # variant LWFBrook90R:
            p_fu_BYFRAC[i] = p_QFFC ^ (1 - (1 / p_QFPAR) * (u_aux_WETNES[i] - p_WETF[i]) / (1 - p_WETF[i]))
            if (p_fu_BYFRAC[i] > 1)
                p_fu_BYFRAC[i] = 1
            end

            # generate bypass flow to avoid saturation
            if (u_aux_WETNES[i] > 0.99)
                p_fu_BYFRAC[i] = 1
            end
            ####

        else
            p_fu_BYFRAC[i] = 0
        end
    end
    return p_fu_BYFRAC
end


""" DSLOP() downslope flow rate from layer
"""
function DSLOP(u_aux_PSIM_i, p_fu_KK_i, p_DSLOPE, p_LENGTH, p_THICK_i, p_STONEF_i, p_RHOWG)
    if( p_LENGTH == 0 || p_DSLOPE == 0)
            # added in Version 4
            return 0
    else
        LL = 1000 * p_LENGTH
        GRAD = p_RHOWG * sin(p_DSLOPE) + (2 * u_aux_PSIM_i / LL) * cos(p_DSLOPE)
        ARATIO = p_THICK_i * (1 - p_STONEF_i) * cos(p_DSLOPE) / LL

        aux_du_DSFLI_i = p_fu_KK_i * ARATIO * GRAD / p_RHOWG

        # no water uptake into dry soil because no free water at outflow face
        if (aux_du_DSFLI_i < 0)
            aux_du_DSFLI_i = 0
        end
        return aux_du_DSFLI_i
    end
end

""" VRFLI() vertical flow rate
"""
function VRFLI(i, NLAYER, u_aux_PSITI, p_fu_KK, p_KSAT, p_THICK, p_STONEF, p_RHOWG, p_DRAIN, p_DPSIMX)
    # out_of_place: allocates a new aux_du_VRFLI_i

    if (i < NLAYER)
        if (abs(u_aux_PSITI[i] - u_aux_PSITI[i+1]) < p_DPSIMX)
            # TODO(bernhard): put this in domain control. However, in the current position it affects
            aux_du_VRFLI_i = 0
        else
            aux_du_VRFLI_i =
                VERT(p_fu_KK[i],    p_fu_KK[i+1],
                    p_KSAT[i],      p_KSAT[i+1],
                    p_THICK[i],     p_THICK[i+1],
                    u_aux_PSITI[i], u_aux_PSITI[i+1],
                    p_STONEF[i],    p_STONEF[i+1],
                    p_RHOWG)
        end
    else
    # bottom layer i == NLAYER
        if (p_DRAIN > 0.00001)
        # gravity drainage only
            aux_du_VRFLI_i = p_DRAIN * p_fu_KK[NLAYER] * (1 - p_STONEF[NLAYER])
        else
        # bottom of profile sealed
            aux_du_VRFLI_i = 0
        end
    end

    return aux_du_VRFLI_i
end

""" VERT() vertical flow rate
"""
function VERT(KK_i, KK_iplus1,
              KSAT_i, KSAT_iplus1,
              THICK_i, THICK_iplus1,
              PSITI_i, PSITI_iplus1,
              STONEF_i, STONEF_iplus1,
              p_RHOWG)

    # NOTE(bernhard): different averaging for KKMEAN and GRAD exist in different Brook90 versions:
    # see: http://www.ecoshift.net/brook/update.htm (Version 4.4 - October 4, 2001)
    # TODO(bernhard):

    # # BROOK90: for Version 4.4 to 4.8 was
    # KKMEAN = exp((log(KK_i) + log(KK_iplus1)) / 2)
    # # BROOK90: for Version 4.2 and 4.3 was
    # KKMEAN = exp((THICK_i * log(KK_i) + THICK1 * log(KK_iplus1)) / (THICK_i + THICK_iplus1))
    # # BROOK90: for Version 4.1 and 3.25a and earlier was
    # KKMEAN = exp( (THICK_iplus1 * log(KK_i) + THICK_i * log(KK_iplus1)) /
    #               (THICK_i + THICK_iplus1) )
    #LWFBROOK90R (2021-02-22):
    KKMEAN = exp( (THICK_iplus1 * log(KK_i) + THICK_i * log(KK_iplus1)) /
                  (THICK_i + THICK_iplus1) )

    # limit KKMEAN to lesser saturated conductivity
    KKMEAN = min(KKMEAN, KSAT_i, KSAT_iplus1)

    # # BROOK90: through Version 4.3a was
    # GRAD = (PSITI_i - PSITI_iplus1) / ((THICK_i + THICK_iplus1) / 2)
    # # BROOK90: for Version 4.4 to 4.8 was
    # GRAD = (PSITI_i - PSITI_iplus1) / min(THICK_i, THICK_iplus1)
    #LWFBROOK90R (2021-02-22):
    GRAD = (PSITI_i - PSITI_iplus1) / ((THICK_i + THICK_iplus1) / 2)

    VRFLI = (GRAD * KKMEAN / p_RHOWG) * (1 - (STONEF_i + STONEF_iplus1) / 2)

    return(VRFLI)
end




"""net inflow to soil layer

Ecoshift:
In this routine, infiltrating water (SLFL) is allocated to soil water in each layer (INFLIi
) and to bypass flow from each layer (BYFLIi ). The fraction of SLFL going to each layer
(INFRACi ) is constant and is obtained in subroutine INFPAR. This fraction is separated into
water to bypass flow (BYFLIi ) and water to the soil matrix (INFLIi ) by the bypass flow
fraction (BYFRACi ) from subroutine BYFLFR. The routine then calculates net inflow to each
layer, including withdrawal by transpiration and soil evaporation. INFLOW is called once
each iteration time step and then is called once again if subroutine ITER produces a new,
shorter, iteration time step.

The INFLOW routine is passed through once for each layer, working from the bottom up,
preventing oversaturation of any layer. The total water allocated to layer i is

INFIL = INFRACi * SLFL.

Bypass flow rate is

BYFLIi = BYFRACi * INFIL

and the infiltration to soil matrix water in the layer is

INFLIi = INFIL - BYFLIi.

INFLOW next determines the maximum inflow rate (MAXINi ), to the layer that can be allowed
in the iteration time-step. The vertical flow to the next layer down (VRFLIi ) (which may be
negative), the transpiration withdrawal (TRANIi ), and the downslope flow (DSFLIi ) are
fixed. So

MAXINi = (SWATMXi - SWATIi) / DTI + VRFLIi + DSFLIi + TRANIi

where DTI is the iteration time step.

If VRFLIi-1 + INFLIi exceeds MAXIN, then oversaturation will occur. If BYFRACi > 0, INFLIi
is first reduced toward zero, then, if necessary, VRFLIi-1 is reduced, or even made negative
if VRFLIi is negative. BYFLIi is increased by the amount that INFLIi is reduced. If BYFRACi
= 0, VRFLIi-1 is reduced or even made negative. INFLOW finally calculates the net water flux
rate, NTFLIi into each soil layer.

NTFLI(I%) = VV(I% - 1) + INFLI(I%) - VV(I%) - DSFLI(I%) - TRANI(I%)

where VV is the final value of VRFLI.

In the top layer, soil evaporation withdrawal is also added to MAXIN. Because there is no
VRFLI(0) to reduce, excess water becomes negative INFLI(1) and increases BYFLI(1).

The modified values of VRFLIi are output from the INFLOW routine as variable VV because the
original VRFLIi are needed again if the iteration time step (DTI) is reduced.
"""
function INFLOW(NLAYER, DTI, p_INFRAC, p_fu_BYFRAC, p_fu_SLFL,
                aux_du_DSFLI, aux_du_TRANI, aux_du_SLVP, p_SWATMX, u_SWATI, VRFLI_prior)
                # This function a) computes all the fluxes involved in the
                # balance of a single soil layer and b) corrects the fluxes of
                # VRFLI, INFLI and BYFLI.
                #
                # Prior guess of VRFLI was computed beforehand.
                # Prior guess of INFLI and BYFLI are computed by this function.
                #
                #
                # It does so by allocating available infilterted water (p_fu_SLFL)
                # to either water supplying a layer (INFLI) or bypass flow (BYFLI)
                # which leaves the system

                ###= for first layer i=1 =======================================
                ###                                       |
                ###                               SLFL*p_INFRAC(1)
                ###                                       |
                ###                                       v
                ###  <-SLVP----- .----------. <-INFLI(1)--+--BYFLI(1)----------->
                ###              |          |             |
                ###  <-TRANI(i)- | SWATI(i) |             ┴ SLFLI(1) (unused)
                ###              |          |
                ###              .----------. --DSFLI(i)------------------------>
                ###                 |
                ###            VRFLI(i)
                ###                 |
                ###                 v
                ###=============================================================
                ###= for any layer i ===========================================
                ###                 |                     |
                ###            VRFLI(i-1)         SLFL*p_INFRAC(i)
                ###                 |                     |
                ###                 v                     v
                ###              .----------. <-INFLI(i)--+--BYFLI(i)----------->
                ###              |          |             |
                ###  <-TRANI(i)- | SWATI(i) |             ┴ SLFLI(1) (unused)
                ###              |          |
                ###              .----------. --DSFLI(i)------------------------>
                ###                 |
                ###            VRFLI(i)
                ###                 |
                ###                 v
                ###=============================================================

    # NLAYER    - number of soil layers being used, max of 20
    # DTI       - time step for iteration interval, d
    # INFRAC(*) - fraction of infiltration to each layer
    # BYFRAC(*) - fraction of layer infiltration to bypass flow
    # SLFL      - input rate to soil surface, mm/d
    # DSFLI(*)  - downslope flow rate from layer, mm/d
    # TRANI(*)  - transpiration rate from layer, mm/d
    # SLVP      - evaporation rate from soil, mm/d
    # SWATMX(*) - maximum water storage for layer, mm
    # SWATI(*)  - water volume in layer, mm
    # VRFLI(*)  - vertical drainage rate from layer, mm/d
    if DTI==0
        #Infiltrator.@infiltrate
        @error ("INFLOW() called with DTI == 0")
    end

    VRFLI_posterior=fill(NaN, NLAYER) #TODO(bernhard): assess whether in-place INFLOW!() is faster
    INFLI=fill(NaN, NLAYER) #TODO(bernhard): assess whether in-place INFLOW!() is faster
    BYFLI=fill(NaN, NLAYER) #TODO(bernhard): assess whether in-place INFLOW!() is faster
    NTFLI=fill(NaN, NLAYER) #TODO(bernhard): assess whether in-place INFLOW!() is faster

    # A) Compute prior estimates for BYFLI and INFLI
    for i = NLAYER:-1:1
        # Allocate available water (SLFL*INFRAC) to either pass by soil (BYFLI)
        # or to infiltrate into the soil layer (INFLI). Or then
        BYFLI[i] = p_fu_SLFL * p_INFRAC[i] * p_fu_BYFRAC[i]
        INFLI[i] = p_fu_SLFL * p_INFRAC[i] * (1.0- p_fu_BYFRAC[i])
    end

    # B) Check if with prior estimates layers become oversaturated and modify
    #    prior estimates accordingly.
    #    For layer i, relevant incoming fluxes are: VRFLi[i-1], BYFLI[i], and INFLI[i]
    #    Therefore VRFLI exiting last layer (VRFLi[NLAYER]) is only outgoing flux
    #    and cannot sature any other soil layer. It is therefore accepted as-is:
    VRFLI_posterior[NLAYER] = VRFLI_prior[NLAYER]

    # Now going through all layers bottom up and prevent oversaturation.
    # (To obey causailty we need to assume flow goes from top to bottom):
    for i = NLAYER:-1:1
        # Compute maximum possible inflow during time interval DTI:
        # maximum allowed rate of input of water to layer, mm/d:
        MAXIN = (p_SWATMX[i] - u_SWATI[i]) / DTI + VRFLI_posterior[i] + aux_du_DSFLI[i] + aux_du_TRANI[i]
        # inflow is constitued by INFLI(i), VRFLI(i-1) for any layer i
        # inflow is constitued by INFLI(1)             for the first layer 1
        if (i == 1)
            # In first layer there is additionally soil evaporation SLVP as an
            # outflow. A fact which increases the maximum possible inflow.
            MAXIN = MAXIN + aux_du_SLVP
            # If inflow is too large
            if (INFLI[1] > MAXIN)
                # Decrease INFLI, and increase BYFLI:
                # variant 1:
                tooMuchInflow = INFLI[1] - MAXIN
                INFLI[1] = INFLI[1] - tooMuchInflow
                BYFLI[1] = BYFLI[1] + tooMuchInflow
                # variant 2:
                #BYFLI[1] = BYFLI[1] + INFLI[1] - MAXIN
                #INFLI[1] = MAXIN
            else
                # Prior estimates of INFLI and BYFLI are accepted for that layer.
                # There is no incoming VRFLI to accept for first layer.
            end
        else
            # i = 2:NLAYER
            if (VRFLI_prior[i - 1] + INFLI[i] > MAXIN)
                # If inflow (VRFLI[i-1] + INFLI[i]) is too large:
                # oversaturation occurs
                # Therefore modify estimation of inflowing VRFLI
                if (p_fu_BYFRAC[i] > 0)
                    # If byflow can occur:
                    if (VRFLI_prior[i - 1] < MAXIN)
                        # If VRFLI on its own is small enough:
                        # a) keep VRFLI as-is,
                        VRFLI_posterior[i-1] = VRFLI_prior[i-1]
                        # b) and balance between BYFLI and INFLI (like in case i=1)
                        # variant 1:
                        tooMuchInflow = INFLI[i] + VRFLI_prior[i - 1] - MAXIN
                        INFLI[i] = INFLI[i] - tooMuchInflow
                        BYFLI[i] = BYFLI[i] + tooMuchInflow
                        # variant 2:
                        #BYFLI[i] = BYFLI[i] + INFLI[i] - (MAXIN - VRFLI_prior[i - 1])
                        #INFLI[i] = MAXIN - VRFLI_prior[i - 1]
                    else
                        # If VRFLI on its own is still too large:
                        # a) reduce VRFLi to max value and b) reroute all of INFLI to BYFLI
                        VRFLI_posterior[i-1] = MAXIN
                        BYFLI[i] = BYFLI[i] + INFLI[i]
                        INFLI[i] = 0
                    end
                else
                    # If no byflow can occur, decrease VRFLI (even making it
                    # negative) so that inflow just achieves complete saturation
                    VRFLI_posterior[i-1] = MAXIN - INFLI[i]
                    # With that VRFLI + INFLI == MAXIN
                    # By modifying VRFLI_posterior[i-1], this will also have an
                    # effect on the next loop iteration. Increasing MAXIN for
                    # the layer above.
                end
            else
                # If inflow is below MAXIN, accept prior estimate of incoming VRFLI for that layer:
                VRFLI_posterior[i-1] = VRFLI_prior[i-1]
                # as well as accept prior estimates of INFLI and BYFLI for that layer.
            end
        end

        if (i > 1)
            NTFLI[i] = VRFLI_posterior[i-1] + INFLI[i] - VRFLI_posterior[i] - aux_du_DSFLI[i] - aux_du_TRANI[i]
        else
            NTFLI[i] = INFLI[i] - VRFLI_posterior[i] - aux_du_DSFLI[i] - aux_du_TRANI[i] - aux_du_SLVP
        end

    end

    # VV(*)     - modified VRFLI, mm/d
    # BYFLI(*)  - bypass flow rate from layer, mm/d
    # INFLI(*)  - infiltration rate into layer, mm/d
    # NTFLI(*)  - net flow rate into layer, mm/d

    return (VRFLI_posterior, INFLI, BYFLI, NTFLI)
end



function SRFLFR(p_QLAYER, u_SWATI, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC)
    SUM = sum(u_SWATI[1:p_QLAYER])

    # Brook90 Variant: if p_QFPAR > 0.01
    # Brook90 Variant:     SAFRAC = min(
    # Brook90 Variant:         1.0,
    # Brook90 Variant:         p_QFFC ^ (1.0 - (1.0 / p_QFPAR) * (SUM - p_SWATQF) / (p_SWATQX - p_SWATQF)))
    # Brook90 Variant: else
    # Brook90 Variant:     # bucket over QLAYERs
    # Brook90 Variant:     if SUM >= p_SWATQF
    # Brook90 Variant:             SAFRAC = 1.0
    # Brook90 Variant:     else
    # Brook90 Variant:             SAFRAC = 0.0
    # Brook90 Variant:     end
    # Brook90 Variant: end

    # LWFBrook90 variant:
    SAFRAC = min(
            1.0,
            p_QFFC ^ (1.0 - (1.0 / p_QFPAR) * (SUM - p_SWATQF) / (p_SWATQX - p_SWATQF)))

    return SAFRAC
end

function ITER(IMODEL, NLAYER, DTI, DTIMIN,
    du_NTFLI, u_aux_PSITI, u_aux_θ,
    u_aux_WETNES, fdpsidwf_ch, fdpsidwf_mvg,
    p_WETINF, p_BEXP, p_PSIF, p_WETF, p_CHM, p_CHN,
    p_MvGα, p_MvGn,
    p_DSWMAX, p_DPSIMX, p_THICK, p_STONEF, p_THSAT, p_θr)
    # ITER() is a step size limiter

    # DTI       ! time step for iteration interval, d
    # DTIMIN    ! minimum time step for iteration interval, d
    # DPSIDW(*) ! rate of change of total potential with water content, kPa/mm
    # NTFLI(*)  ! net flow rate into layer, mm/d
    # PSITI(*)  ! total potential, kPa
    # u_aux_θ   ! volumetric soil moisture content, m3/m3
    # p_DSWMAX  ! maximum change allowed in SWATI, percent of SWATMX(i)
    # p_DPSIMX  ! maximum potential difference considered "equal", kPa
    # p_THICK   ! soil layer thicknesses, mm
    # p_THSAT   ! θ at saturation == matrix porosity (-)
    # unused: p_SWATMX   maximum water storage for layer, mm

    if IMODEL == 0
        DPSIDW = fdpsidwf_ch(u_aux_WETNES, p_WETINF, p_BEXP, p_PSIF, p_WETF, p_CHM, p_CHN)
    else # IMODEL == 1
        DPSIDW = fdpsidwf_mvg(u_aux_WETNES, p_MvGα, p_MvGn)
    end

    A = zeros(NLAYER)
    temp = zeros(NLAYER)
    # first approximation to new total potential
    if (IMODEL == 0)
        for i = 1:NLAYER
            # A[i]    = du_NTFLI[i] * DPSIDW[i] / p_SWATMX[i]
            # NOTE(bernhard): as p_SWATMX[i] = p_THICK[i] * p_THSAT[i] * (1 - p_STONEF[i])
            A[i]    = du_NTFLI[i]/p_THICK[i] * DPSIDW[i] / (p_THSAT[i] -     0. )/ (1 - p_STONEF[i])
            temp[i] = u_aux_PSITI[i] + A[i] * DTI
        end
    elseif (IMODEL == 1)
        for i = 1:NLAYER
            A[i]    = du_NTFLI[i]/p_THICK[i] * DPSIDW[i] / (p_THSAT[i] - p_θr[i])
            # TODO(bernhard): is ther no STONEF in IMODEL==1. Bug?
            temp[i] = u_aux_PSITI[i] + A[i] * DTI
        end
    else
        error("IMODEL unknown.")
    end

    # test to see if DTI should be reduced
    DTINEW = DTI

    for i = 1:NLAYER
        # 1) prevent too large a change in water content, reduce DTI to keep change below p_DSWMAX
        # DTINEW = min(DTINEW, 0.01 * p_DSWMAX * p_SWATMX[i] / max(0.000001, abs(du_NTFLI[i])))
        DTINEW = min(DTINEW,   0.01 * p_DSWMAX * p_THICK[i] * p_THSAT[i] * (1 - p_STONEF[i]) / max(0.000001, abs(du_NTFLI[i])))
        # TODO(bernhard): shouldn't here also be checked that DTINEW=max(DTIMIN, DTINEW) ?

        # 2) If water is flowing out of cell (du_NTFLI < 0), prevent
        #    a change in water content larger than total available water
        if (du_NTFLI[i] < 0)
            if (IMODEL == 0)
                available_water = (0.0     - u_aux_θ[i]) * p_THICK[i] # TODO(bernhard): is ther no STONEF in IMODEL==0. Bug?
            else # IMODEL == 1
                available_water = (p_θr[i] - u_aux_θ[i]) * p_THICK[i] # TODO(bernhard): is ther no STONEF in IMODEL==1. Bug?
            end

            DTINEW = min(DTINEW, available_water/du_NTFLI[i]/1.30)
            if (DTINEW < DTIMIN)
                # Bernhard: if debug: LWFBrook90R printed here full state vector
                DTINEW = DTIMIN

                # NOTE(Bernhard): if step is smaller than DTMIN, reduce TRANI and SLVP, but keep at least DTMIN
                # error("DTINEW is smaller than DTMIN. LWFBrook90R reduce in this case TRANI and SLVP. This is not implemented in LWFBrook90Jullia.")
                # TRANI[i] = 0 # TODO(Bernhard): shold this change leak out into main program? (side effect)
                # if (i == 1)
                #     SLVP=0   # TODO(Benrhard): should this change leak out into main program? (side effect)
                # end
                # NOTE: This original TRANI and SLVP correction violates the mass balance.
                # @warn "Reduced DTI was lower than DTIMIN. DTI was increased to DTIMIN. Warning: original Brook set TRANI and SLVP to zero in these cases. This is not done anymore."
            end
        end

        # 3) prevent oscillation of gradient in soil water potential
        if (IMODEL == 0)
            if (i < NLAYER)
                # total potential difference at beginning of iteration
                PP = u_aux_PSITI[i] - u_aux_PSITI[i + 1]
                # first approx to total potential difference at end of iteration
                TT = temp[i] - temp[i + 1]
                if ((abs(TT) > p_DPSIMX) && (abs(PP) > p_DPSIMX) && (sign(TT) != sign(PP)))
                    DTINEW = min(DTINEW, -PP / (A[i] - A[i + 1]))
                    DTINEW = max(DTINEW, DTIMIN) # Only in LWFBrook90R
                end
            end
        end
    end

    # if (DTINEW < DTI)
    #     @info("DTI reduced from $DTI to $DTINEW.")
    # end

    return DTINEW # return second estimate of DTI
end

"""calculates groundwater flow and seepage loss
"""
function GWATER(u_GWAT, p_GSC, p_GSP, p_DT, aux_du_VRFLI_N)
    if (p_GSC < 1.0e-8)
        # no groundwater
        du_SEEP = p_GSP * aux_du_VRFLI_N
        du_GWFL = aux_du_VRFLI_N - du_SEEP
    else
        du_SEEP = u_GWAT * p_GSC * p_GSP
        du_GWFL = u_GWAT * p_GSC * (1. - p_GSP)

        # prevent negative GWAT
        # TODO(bernhard): put this in domain control
        #                 Actually no, as this has no feedback on any other variable
        if (u_GWAT - (du_GWFL + du_SEEP)*p_DT < 0)
            # if groundwater would be more than emptied
            # reduce fluxes so that it is just emptied completely during p_DT
            du_SEEP = u_GWAT / p_DT * p_GSP
            du_GWFL = u_GWAT / p_DT * (1. - p_GSP)
        end
    end
    return (du_GWFL, du_SEEP)
end


end