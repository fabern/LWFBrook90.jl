# fabian.bernhard@wsl.ch, 2021-02-03

module WAT # WATER MOVEMENT IN SOIL

export INFPAR, LWFRootGrowth

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


end