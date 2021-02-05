# fabian.bernhard@wsl.ch, 2021-02-03

module WAT # WATER MOVEMENT IN SOIL

export INFPAR

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

end