# fabian.bernhard@wsl.ch, 2021-02-03

module SUN # RADIATION

export EQUIVSLP

# Ecoshift-SUN:
# Daily solar radiation on a horizontal surface (SOLRAD in MJ/m2) is an input variable to
# BROOK90. This value is sometimes called global radiation to emphasize that it includes
# both direct or beam radiation from the sun and diffuse radiation from the sky hemisphere.
# It is directly measured by various types of pyranometers at most research sites and some
# locations of the National Weather Service. SOLRAD is internally prevented from being
# larger than 0.99 times the potential insolation (in the absence of an atmosphere), I0HDAY.

# If no data are available and input SOLRAD is zero, BROOK90 uses a fixed fraction of the
# potential solar radiation on a horizontal surface. This fraction had been set at 0.55, but
# as of Version 4.8 it can be changed on the BROOK90 main window. Use of this generalized
# fraction loses the effects of day-to-day variation of solar radiation in the model.

# Alternatively, there are various simple to complicated methods for estimating SOLRAD from
# cloud cover measurements, which are more widely made (Brutsaert 1982, Dingman 1994).

# Another approach to estimating daily solar radiation uses the fact that the daily
# temperature range is usually larger on clear days than on overcast ones. The equation of
# Bristow and Campbell (1984) has been used worldwide in various modifications. For Hubbard
# Brook I use

# R / I0 = a (1 - exp(-b Δ ^ c)) + d

# where R is the daily solar radiation, I0 is the potential insolation (I0HDAY), and Δ is
# the daily temperature range (Tmax - Tmin), and a, b, c, and d are constants that depend
# more or less on location. I added the d value so that R is greater than zero when Δ is
# zero. For Hubbard Brook I use a = 0.78, c = 1.5, d = 0.05, and b = 0.245 + 0.0035 sin(2π
# (doy + 80) / 366) to account for some residual seasonality, where doy is day of the year.

# In BROOK90 subroutine EQUIVSLP is called once to obtain parameters that depend on
# latitude, slope, and aspect of the surface being simulated. Subroutine SUNDS uses
# functions FUNC3 and HAFDAY once a day to calculate daylength, potential radiation on a
# horizontal surface, and the ratio of potential radiation on the slope to that on a
# horizontal surface. AVAILEN calculates net radiation minus soil heat flux at the top and
# the bottom of the canopy, separately for daytime and for nighttime.

"""
Correction of solar radiation for slope and aspect requires calculation of an "equivalent
slope" in subroutine EQUIVSLP. The equivalent slope is defined as the location on the
earth's surface where a horizontal surface is parallel to the given sloping surface.
Following Swift (1976), L1 is the latitude of this "equivalent slope" (which is actually a
horizontal surface), and L2 is the difference in hour angle (longitude) between the two
locations. For any given slope and aspect, L1 and L2 need be found only once; they do not
change over time. So they are calculated in EQUIVSLP at the beginning of B90, as:

L1 = ASIN[ COS(SLOPE) * SIN(LAT) + SIN(SLOPE) * COS(LAT)* COS(ASPECT) ]
and
L2 = ATAN { SIN(SLOPE) * SIN(ASPECT) /[ COS(SLOPE) * COS(LAT) - SIN(SLOPE) * SIN(LAT) * COS(ASPECT) ] }

with fixes for negative or zero denominator in L2, where SLOPE (ESLOPE), ASPECT, and
latitude (LAT) are input parameters describing the location. All angles in the subroutine
are in radians.
"""
function EQUIVSLP(p_LAT, SLOPE, p_ASPECT)
    #Swift#s p_L1 and p_L2, Lee (3.31, 3.32)
    p_L1 = asin(cos(SLOPE) * sin(p_LAT) + sin(SLOPE) * cos(p_LAT) * cos(p_ASPECT))
    D1 = cos(SLOPE) * cos(p_LAT) - sin(SLOPE) * sin(p_LAT) * cos(p_ASPECT)
    if (D1 == 0)
        D1 = 0.0000000001
    end
    p_L2 = atan(sin(SLOPE) * sin(p_ASPECT) / D1)
    if (D1 < 0)
        p_L2 = p_L2 + LWFBrook90Julia.CONSTANTS.p_PI
    end
    return (p_L1, p_L2)
end

end