using Interpolations: interpolate, extrapolate, NoInterp, Gridded, Constant, Next, Previous, Flat, Throw, scale, BSpline
"""
    define_diff_eq_parameters()

Generate vector p needed for ODE() problem in DiffEq.jl package.

# Arguments
- `continuous_SPAC::SPAC`: Instance of a SPAC-model (soil-plant-atmosphere continuum)
- `compute_intermediate_quantities::...`: TODO argument description.
- `simulate_isotopes::...`: TODO argument description.
- `soil_output_depths`: vector of depths at which state variables should be extractable (negative numeric values [in meter])
"""
function define_LWFB90_p(continuous_SPAC::SPAC, soil_discr)
    ########
    ## Interpolate discretized root distribution
     p_RELDEN = LWFBrook90.interpolate_spaceDiscretized_root_density(;
        timepoints = continuous_SPAC.meteo_forcing.p_days,
        p_AGE      = continuous_SPAC.canopy_evolution.p_AGE,
        p_INITRDEP = continuous_SPAC.params[:INITRDEP],
        p_INITRLEN = continuous_SPAC.params[:INITRLEN],
        p_RGROPER  = continuous_SPAC.params[:RGROPER],
        p_tini     = soil_discr["tini"],
        p_frelden  = soil_discr["frelden"]);
    # TODO(bernhard): document input parameters: INITRDEP, INITRLEN, RGROPER, tini, frelden, MAXLAI, HEIGHT_baseline_m
    # TOOD(bernhard): remove from params: IDEPTH_m, QDEPTH_m, INITRDEP, RGRORATE, INITRDEP, INITRLEN, RGROPER
    ########

    ########
    ## Solver algorithm options
    p_DPSIMAX = continuous_SPAC.solver_options[:DPSIMAX] # maximum potential difference considered "equal", kPa               (BROOK90: DPSIMAX is fixed at 0.01 kPa)
    p_DSWMAX  = continuous_SPAC.solver_options[:DSWMAX] # maximum change allowed in soil wetness SWATI, percent of SWATMAX(i) (BROOK90: DSWMAX is fixed at 2 %)
    p_DTIMAX  = continuous_SPAC.solver_options[:DTIMAX] # maximum iteration time step, d (BROOK90: DTIMAX is fixed at 0.5 d)
    # Documentation from ecoshift:
    # DPSIMAX (Fixed parameter) - maximum potential difference considered equal during soil water integration, kPa. There is no vertical flow between layers whose potentials differ by less than DPSIMAX. This reduces oscillation initiated by flows that are the product of large conductivities and large time steps, but small gradients. The number of iterations used is not at all linearly related to the three iteration parameters DPSIMAX, DSWMAX, and DTIMAX. Selection of values depends on whether the user only wants monthly or daily totals, or is concerned with behaviour at shorter time steps. Generally, faster runs are obtained by using fewer thicker soil layers rather than by using large values of DSWMAX and DPSIMAX. DPSIMAX is fixed at 0.01 kPa. [see WAT-ITER]
    # DSWMAX (Fixed parameter) - maximum change allowed in soil wetness for any layer during an iteration, percent. DSWMAX sets the maximum change in soil wetness or saturation fraction (SWATI / SWATMAX(i)) allowed for any layer in an iteration. See also DPSIMAX. DSWMAX is fixed at 2 %. [see WAT-ITER]
    # DTIMAX (Fixed parameter) - maximum iteration time step, d. DTIMAX is fixed at 0.5 d, which forces at least two iterations per day. This is the largest value that should be used. Much smaller values of DTIMAX, between 0.01 and 0.001 d, will force many iterations per day and thus smooth integration. However, a run with such a small DTIMAX takes a long time. See also DPSIMAX. [see WAT-ITER]


    ## Heat flow (unimplemented)
    # p_HEAT    = continuous_SPAC.params[:HEAT]
    # p_TopInfT = soil_discr["TopInfT"]
    # unused p_HeatCapOld = soil_discr["HeatCapOld"]

    # Isotope transport parameters
    p_VXYLEM       = continuous_SPAC.params[:VXYLEM_mm] # mm, storage volume of well mixed xylem storage per ground area # TODO(bernhard): possibly link this to SAI...
    p_DISPERSIVITY = continuous_SPAC.params[:DISPERSIVITY_mm]/1000  # m, dispersivity length (0.04m is the average fitted dispersivity to the lysimters of Stumpp et al. 2012)

    ## Location / Meteo
    p_NPINT  = 1 # Hardcoded. If p_NPINT>1, then multiple precipitation intervals would need
                 #            to be defined in an additional input data set PRECDAT.
    if p_NPINT == 1
        p_DTP = 1 / p_NPINT
    else
        error("Case with multiple precipitation intervals (using PRECDAT and precip_interval != 1) is not implemented.")
    end
    p_DURATN = continuous_SPAC.storm_durations[1:12,"storm_durations_h"]# average storm duration for each month, hr
    p_LAT    = continuous_SPAC.params[:LAT_DEG]   /57.296  # (Location parameter), latitude, radians
    p_ASPECT = continuous_SPAC.params[:ASPECT_DEG]/57.296  # (Location parameter), aspect, radians through east from north
    p_ESLOPE = continuous_SPAC.params[:ESLOPE_DEG]/57.296  # (Location parameter), slope for evapotranspiration and snowmelt, radians
    p_MELFAC = continuous_SPAC.params[:MELFAC]    # (Location parameter), degree day melt factor for open land, MJ m-2 d-1 K-1
    p_RSTEMP = continuous_SPAC.params[:RSTEMP]    # (Location parameter), base temperature for snow-rain transition, °C
    # removed for LWFBrook90: RELHT
    # removed for LWFBrook90: RELLAI
    # Documentation from ecoshift:
    # DURATN(1 To 12) (Location parameter) - average duration of daily precipitation by month, hr. When a precipitation input file is not used, DURATN is used in subroutine INTER24 to estimate interception processes hourly. It is the average number of hours per day with precipitation. Do not use a DURATN of 1, or no interception will be produced. Fractional values are truncated to the next lower even integer. Provision of 12 monthly values for DURATN allows seasonal variation of storm type from low intensity, long DURATN, to high intensity, short DURATN. In reality, though, this seasonal difference is not large. The number of hours a day that precipitation exceeds 0.5 mm varies only from 2 to 6 for most seasons and locations in the United States (see EVP-INTER24. DURATN = 4 for all months is a satisfactory approximation. DURATN is ignored when a precipitation interval file is provided; then precipitation is assumed constant over the precipitation interval, and subroutine INTER is used. [see EVP-INTER24]
    # LAT and LAT_DEG (Location parameter)       - latitude, degrees.                                                          In code LAT_DEG    is in degrees, LAT    is in radians. LAT is the location latitude in degrees north for solar radiation calculations Negative values are used in the southern hemisphere. This is the only parameter that is specific to geographic location. [see SUN-EQUIVSLP]
    # ASPECT and ASPECT_DEG (Location parameter) - aspect, degrees through east from north. Used in radiation calculations.    In code ASPECT_DEG is in degrees, ASPECT is in radians. If ESLOPE = 0, ASPECT is ignored. [see SUN-EQIVSLP]
    # ESLOPE and ESLOPE_DEG (Location parameter) - slope for evapotranspiration and snowmelt, degrees.                         In code ESLOPE_DEG is in degrees, ESLOPE is in radians. ESLOPE is used for net radiation and snowmelt. See also DSLOPE. [see SUN-EQUIVSLP]
    # MELFAC (Location parameter) - degree day melt factor for open land, MJ m-2 d-1 K-1. MELFAC is the degree-day snowmelt factor for a day with daylength of 0.5 d and no plant canopy. Dividing MELFAC by the heat of fusion of water, LF (= 0.335 MJ m-2 mm-1), gives the more usual degree day factor in mm d-1 K-1. MELFAC = 1.5 MJ m-2 d-1 K-1, or 4.5 mm d-1 K-1, is a good starting value; it is in the lower end of the range given by Federer et al. (1973) and is close to the 4.2 mm d-1 K-1 used by Anderson (1976). Lower values will slow snowmelt. To turn off SMLT, set MELFAC to zero. [see SNO-SNOENRGY]
    # RSTEMP (Location parameter) - base temperature for snow-rain transition, °C. The fraction of daily precipitation as snow is (RSTEMP - TMIN) / (TMAX - TMIN) when TMIN < RSTEMP < TMAX. Lowering RSTEMP produces more snow. RSTEMP = -100 turns off SFAL and makes all precipitation rain. [see SNO-SNOFRAC]
    # RELHT(1 To 20) (Location parameter) - array of ten pairs of day of the year (DOY) and relative height between 0 and 1, dimensionless. RELHT allows canopy height to vary through the year and is a multiplier of MAXHT. The DOY values must increase; intermediate values are interpolated linearly. The first DOY used must be 1 and a later DOY value must be 366. Any remaining pairs after DOY 366 are ignored. RELHT different from 1 would generally be used for herbaceous plant covers. The phenology of RELHT must be determined externally to BROOK90, and is location as well as cover type dependent. [see PET-CANOPY]
    # RELLAI(1 To 20) (Location parameter) - array of ten pairs of DOY and relative LAI between 0 and 1, dimensionless. RELLAI allows LAI to vary through the year and is a multiplier of MAXLAI. The doy values must increase; intermediate values are interpolated linearly. The first DOY must be 1 and a later DOY value must be 366. Any remaining pairs after DOY 366 are ignored. The phenology of RELLAI must be determined externally to BROOK90, and is location as well as cover type dependent. [see PET-CANOPY]

    # L1 and L2 - latitude of equivalent slope, radians, and time shift of equivalent slope, radians
    p_L1, p_L2 = LWFBrook90.SUN.EQUIVSLP(p_LAT, p_ESLOPE, p_ASPECT)


    ## Meteo
    p_FETCH  = continuous_SPAC.params[:FETCH]  # Fetch upwind of the weather station at which wind speed was measured, m (BROOK90: FETCH is fixed at 5000 m)
    p_WNDRAT = continuous_SPAC.params[:WNDRAT] # Average ratio of nighttime to daytime wind speed, dimensionless (BROOK90: WNDRAT is fixed at 0.3)
    p_Z0G    = continuous_SPAC.params[:Z0G]    # Ground surface roughness, m ()
    p_Z0S    = continuous_SPAC.params[:Z0S]    # Roughness parameter of snow surface, m (BROOK90: Z0S is fixed at 0.001 m)
    p_Z0W    = continuous_SPAC.params[:Z0W]    # Roughness parameter at the weather station where wind speed was measured, m
    p_ZMINH  = continuous_SPAC.params[:ZMINH]  # Reference height for weather data above canopy HEIGHT, m (BROOK90: ZMINH is fixed at 2 m)
    p_ZW     = continuous_SPAC.params[:ZW]     # Weather station measurement height for wind, m (BROOK90: W is fixed at 10 m)
    p_HS     = continuous_SPAC.params[:HS]     # Lower height limit, for roughness parameter interpolation, if canopy HEIGHT is below -> CZS (BROOK90: HS is fixed at 1 m)
    p_HR     = continuous_SPAC.params[:HR]     # Upper height limit, for roughness parameter interpolation, if canopy HEIGHT is above -> CZR (BROOK90: HR is fixed at 10 m)
    p_CZS    = continuous_SPAC.params[:CZS]    # Ratio of roughness parameter to HEIGHT (canopy) for HEIGHT < HS, dimensionless (BROOK90: CZS is fixed at 0.13)
    p_CZR    = continuous_SPAC.params[:CZR]    # Ratio of roughness parameter to HEIGHT (canopy) for HEIGHT > HR, dimensionless (BROOK90: CZR is fixed at 0.05)
    p_C1     = continuous_SPAC.params[:C1]     # intercept of linear relation between the ratio of actual to potential solar radiation for the day and sunshine duration (BROOK90: C1 is fixed at 0.25 following Brutsaert (1982))
    p_C2     = continuous_SPAC.params[:C2]     # slope     of linear relation between the ratio of actual to potential solar radiation for the day and sunshine duration (BROOK90: C2 is fixed at 0.5 following Brutsaert (1982))
    p_C3     = continuous_SPAC.params[:C3]     # Ratio of net longwave radiation for overcast sky (sunshine duration = 0) to that for clear sky (sunshine duration = 1). (BROOK90: C3 is fixed at 0.2 following Brutsaert (1982))
    # Documentation from ecoshift:
    # C3 (Fixed parameter) - ratio of net longwave radiation for overcast sky (sunshine duration = 0) to that for clear sky (sunshine duration = 1). C3 is fixed at 0.2 following Brutsaert (1982). [see SUN-AVAILEN]
    # FETCH (Fixed parameter) - fetch upwind of the weather station at which wind speed was measured, m. Sensitivity to FETCH is small if it is > 1000 m. See also Z0W. FETCH is fixed at 5000 m. FETCH is ignored if Z0W = 0. [see PET-WNDADJ]
    # WNDRAT (Fixed parameter) - average ratio of nighttime to daytime wind speed, dimensionless. WNDRAT is fixed at 0.3. The value is based on Hubbard Brook data and could be changed if a better value were known. [see PET-WEATHER]
    # Z0G (Fixed parameter) - ground surface roughness, m. Z0G is the roughness parameter of the ground surface below the canopy. It controls the amount of turbulent transfer at the ground surface, and thus soil evaporation. It should generally be on the order of a few centimeters, but the value probably doesn't matter much. See also Z0S. In code, Z0GS is Z0S when snow is present and Z0G with no snow. [see PET-ROUGH] [see PET-SWGRA]
    # Z0S (Fixed parameter) - roughness parameter of snow surface, m. Z0S used in place of Z0G when there is snow on the ground. Z0S is fixed at 0.001 m. In code, Z0GS is Z0S when snow is present and Z0G with no snow. [see PET-ROUGH]
    # Z0W (Fixed parameter) - roughness parameter at the weather station where wind speed was measured, m. Used to adjust UW to wind speed at reference height above the surface being simulated. If wind speed was actually measured above the surface being simulated set Z0W = 0, then UA = UW and FETCH and ZW are ignored. Z0W of 0.005 m is appropriate for a short grass surface. [see PET-WNDADJ]
    # ZMINH (Fixed parameter) - reference height for weather data above canopy HEIGHT, m. ZA is ZMINH + RELHT * MAXHT, which can vary seasonally. ZMINH is fixed at 2 m. [see PET-ROUGH]
    # ZW (Fixed parameter) - weather station measurement height for wind, m. ZW is the height above the ground at which wind speed (UW) was measured. ZW is fixed at 10 m. Unless more is known about the wind speed weather station, ZW, FETCH, and Z0W should not be changed. ZW is ignored if Z0W = 0. [see PET-WNDADJ]
    # CZR, HR, CZS, and HS (Fixed parameters) - CZR is the ratio of roughness parameter to HEIGHT when HEIGHT is greater than HR and canopy is closed (LAI > LPC).
                                                # CZS is the ratio of roughness parameter to HEIGHT when HEIGHT is less than HS and canopy is closed.
                                                # For heights between HS and HR, the roughness parameter is interpolated.
                                                # CZR is fixed at 0.05 and HR is fixed at 10 m.
                                                # CZS is fixed at 0.13 and HS is fixed at 1 m. [see PET-ROUGH]
    # C1 and C2 (Fixed parameters) - intercept and slope of linear relation between the ratio of actual to potential solar radiation for the day and sunshine duration. C1 is fixed at 0.25 and C2 at 0.5 following Brutsaert (1982). [see SUN-AVAILEN]


    ## Snowpack
    # SNODEN is the snowpack density or the ratio of water content to depth. SNODEN is used only to correct leaf area index and stem area index to the fraction of them that is above the snow. It does not vary with time or snowpack ripeness in BROOK90. [see PET-CANOPY]
    p_SNODEN = continuous_SPAC.params[:SNODEN] # snowpack density, mm/mm (BROOK90: SNODEN is fixed at 0.3)
    p_CCFAC  = continuous_SPAC.params[:CCFAC]  # cold content factor, MJ m-2 d-1 K-1 (BROOK90: CCFAC is fixed at 0.3 MJ m-2 d-1 K-1)
    p_GRDMLT = continuous_SPAC.params[:GRDMLT] # rate of groundmelt of snowpack, mm/d (BROOK90: GRDMLT is fixed at 0.35 mm/d from Hubbard Brook (Federer 1965))
    p_LAIMLT = continuous_SPAC.params[:LAIMLT] # dependence of snowmelt on LAI, dimensionless (BROOK90: LAIMLT is fixed at 0.2)
    p_SAIMLT = continuous_SPAC.params[:SAIMLT] # dependence of snowmelt on SAI, dimensionless (BROOK90: SAIMLT is fixed at 0.5)
    p_MAXLQF = continuous_SPAC.params[:MAXLQF] # maximum liquid water fraction of SNOW, dimensionless (BROOK90: MAXLQF is fixed at 0.05)
    # Documentation from ecoshift:
    # SNODEN (Fixed parameter) - snow density, mm/mm. SNODEN is the snowpack density or the ratio of water content to depth. SNODEN is used only to correct leaf area index and stem area index to the fraction of them that is above the snow. It does not vary with time or snowpack ripeness in BROOK90. SNODEN is fixed at 0.3. [see PET-CANOPY]
    # CCFAC (Fixed parameter) - cold content factor, MJ m-2 d-1 K-1. CCFAC is a degree day factor for accumulation of cold content for a day with daylength of 0.5 d. It controls the snow energy balance when TA is less than 0°C. Larger values make the snow temperature lag farther behind the air temperature. Sensitivity of snowmelt to CCFAC is small unless CCFAC is less than 0.05 MJ m-2 d-1 K-1. CCFAC = 0 means there is no cold content and snow temperature is always 0°C. CCFAC is fixed at 0.3 MJ m-2 d-1 K-1. [see SNO-SNOENRGY]
    # GRDMLT (Fixed parameter) - rate of groundmelt of snowpack, mm/d. GRDMLT is the constant rate of melt at the bottom of the snowpack because of soil heat transfer to the snow. This parameter controls low flow during winter periods with snowpack. BROOK90 assumes that there is never any soil frost. GRDMLT is fixed at 0.35 mm/d from Hubbard Brook (Federer 1965). [see SNO-SNOWPACK]
    # LAIMLT (Fixed parameter) - dependence of snowmelt on LAI, dimensionless. Melt is linearly proportional to exp(-LAIMLT * LAI) so melt decreases exponentially as LAI increases. LAIMLT is fixed at 0.2. This value makes the exp factor 1.0, 0.67, 0.45, and 0.30 for LAI = 0, 2, 4, and 6 respectively. [see SNO-SNOENRGY]
    # SAIMLT (Fixed parameter) - snowmelt dependence on SAI, dimensionless. Melt is linearly proportional to exp (-SAIMLT * SAI); so melt decreases exponentially as SAI increases. SAIMLT is fixed at 0.5. With SAIMLT = 0.5 and LAIMLT = 0.2 , melt in deciduous forest (LAI=0, SAI=0.7) is 0.70 times that in the open, and melt in conifer forest (LAI=6, SAI=0.7) is 0.21 times that in the open. Federer et al. (1973b) generalized literature ratios to 0.5 for deciduous/open and 0.25 for conifer/open. Matching their results with the exponential algorithm would require greater difference between LAIMLT and SAIMLT, but it is difficult to see why they should differ much. [see SNO-SNOENRGY]
    # MAXLQF (Fixed parameter) - maximum liquid water fraction of SNOW, dimensionless. MAXLQF is the liquid water fraction of the snow water, SNOW, at which water drains. MAXLQF is fixed at 0.05, which is appropriate in most snow environments. [see SNO-SNOWPACK]


    ## Vegetation (canopy)
    ### Vegetation dimensions
    p_CR     = continuous_SPAC.params[:CR]     # (Canopy parameter), ratio of projected stem area index (SAI) to HEIGHT when DENSEF = 1, (SAI = CS * HEIGHT * DENSEF)
    p_RHOTP  = continuous_SPAC.params[:RHOTP]  # (Fixed  parameter), Ratio of total leaf area to projected area, dimensionless (BROOK90: RHOTP is fixed at 2)
    p_LWIDTH = continuous_SPAC.params[:LWIDTH] # (Canopy parameter), average leaf width, m
    p_LPC    = continuous_SPAC.params[:LPC]    # (Fixed  parameter), Minimum LAI defining a closed canopy, dimensionless (BROOK90: LPC is fixed at 4.0)
    ### Vegetation influence on atmosphere and meteorology
    p_KSNVP  = continuous_SPAC.params[:KSNVP]  # (Canopy parameter), reduction factor to reduce snow evaporation (SNVP), (0.05 - 1)
    p_ALBSN  = continuous_SPAC.params[:ALBSN]  # (Canopy parameter), albedo or surface reflectivity with snow on the ground, (typically 0.1-0.9)
    p_ALB    = continuous_SPAC.params[:ALB]    # (Canopy parameter), albedo or surface reflectivity without snow on the ground, (typically 0.1-0.3)
    p_NN     = continuous_SPAC.params[:NN]     # (Fixed  parameter), Wind/diffusivity extinction coefficient, dimensionless (BROOK90: NN is fixed at 2.5 following Shuttleworth and Gurney (1990))
    p_CS     = continuous_SPAC.params[:CS]     # (Canopy parameter), extinction coefficient for photosynthetically-active radiation in the canopy, (unit?)
    p_RM     = continuous_SPAC.params[:RM]     # (Fixed  parameter), Nominal maximum solar radiation possible on a leaf, W/m2 (BROOK90: RM is fixed at 1000 W/m2)
    ### Interception
    p_CINTRL = continuous_SPAC.params[:CINTRL] # (Fixed  parameter), Maximum interception storage of rain       per unit LAI, mm (BROOK90: CINTRL and CINTRS are both fixed at 0.15 mm)
    p_CINTRS = continuous_SPAC.params[:CINTRS] # (Fixed  parameter), Maximum interception storage of rain       per unit SAI, mm (BROOK90: CINTRL and CINTRS are both fixed at 0.15 mm)
    p_CINTSL = continuous_SPAC.params[:CINTSL] # (Fixed  parameter), Maximum interception storage of snow water per unit LAI, mm (BROOK90: CINTSL and CINTSS are both fixed at 0.6 mm)
    p_CINTSS = continuous_SPAC.params[:CINTSS] # (Fixed  parameter), Maximum interception storage of snow water per unit SAI, mm (BROOK90: CINTSL and CINTSS are both fixed at 0.6 mm)
    p_FRINTL = continuous_SPAC.params[:FRINTLAI] # (Fixed  parameter), Intercepted fraction of rain per unit LAI, dimensionless (BROOK90: FRINTLAI and FRINTSAI are both fixed at 0.06)
    p_FRINTS = continuous_SPAC.params[:FRINTSAI] # (Fixed  parameter), Intercepted fraction of rain per unit SAI, dimensionless (BROOK90: FRINTLAI and FRINTSAI are both fixed at 0.06)
    p_FSINTL = continuous_SPAC.params[:FSINTLAI] # (Fixed  parameter), Intercepted fraction of snow per unit LAI, dimensionless (BROOK90: FSINTLAI and FSINTSAI are both fixed at 0.04)
    p_FSINTS = continuous_SPAC.params[:FSINTSAI] # (Fixed  parameter), Intercepted fraction of snow per unit SAI, dimensionless (BROOK90: FSINTLAI and FSINTSAI are both fixed at 0.04)
    ### Vegetation conductivity
    p_MXKPL  = continuous_SPAC.params[:MXKPL]  # (Canopy parameter), maximum plant conductivity, mm d-1 MPa-1.
    p_FXYLEM = continuous_SPAC.params[:FXYLEM] # (Canopy parameter), fraction of plant resistance that is in the xylem (above ground) and not in the roots, (0-1)
    p_GLMAX  = continuous_SPAC.params[:GLMAX]  # (Canopy parameter), maximum leaf conductance, cm/s
    p_GLMIN  = continuous_SPAC.params[:GLMIN]  # (Canopy parameter), minimum leaf conductance, cm/s
    ### Stomatal regulation
    p_TH     = continuous_SPAC.params[:TH]     # (Canopy parameter), temperature controlling closing of stomates,°C
    p_T1     = continuous_SPAC.params[:T1]     # (Canopy parameter), temperature controlling closing of stomates,°C
    p_T2     = continuous_SPAC.params[:T2]     # (Canopy parameter), temperature controlling closing of stomates,°C
    p_TL     = continuous_SPAC.params[:TL]     # (Canopy parameter), temperature controlling closing of stomates,°C
    p_CVPD   = continuous_SPAC.params[:CVPD]   # (Fixed  parameter), Vapor pressure deficit at which stomatal conductance is halved, kPa (BROOK90: CVPD is fixed at 2 kPa for all cover types)
    p_R5     = continuous_SPAC.params[:R5]     # (Fixed  parameter), Solar radiation at which stomatal conductance is half of its value at RM, W/m2 (BROOK90: BROOK90 fixes R5 = 100 W/m2 as the default for all cover types)
    p_PSICR  = continuous_SPAC.params[:PSICR]  # (Canopy parameter), minimum plant leaf water potential, MPa.
    # Documentation from ecoshift:
    # MXKPL The internal resistance to water flow through the plants RPLANT = 1 / (MXKPL * RELHT * DENSEF). MXKPL is the main controller of soil-water availability and is a property of all of the plants on a unit area, not of any one plant. When the canopy is at its maximum seasonal LAI and height, and when soil is wet so that soil water potential is effectively zero, and when the leaf water potential is at its critical value, PSICR, then MXKPL is the transpiration rate divided by PSICR. Abdul-Jabbar et al. (1988) found that MXKPL ranges only from 7 to 30 mm d-1 MPa-1 over a wide range of vegetations, a surprisingly constant parameter. A transpiration rate of 0.5 mm/hr at a gradient of -1.5 MPa is typical, giving MXKPL of 8 mm d-1 MPa-1 (Hunt et al. 1991). MXKPL controls the rate of water supply to the leaves and thus the transpiration when soil water supply is limiting. Decreasing MXKPL makes soil water less available and thus reduces actual transpiration below potential transpiration at higher soil water content. [see PET-CANOPY] [see EVP-PLNTRES]
    # CR (Canopy parameter) - extinction coefficient for photosynthetically-active radiation in the canopy. Values usually range from 0.5 to 0.7. Values outside this range should be used very cautiously. The extinction coefficient can be determined from the canopy transmissivity, t, as CR = - (ln t) / (Lp + Sp). For a canopy of Lp = 6 and Sp = 0.7, PAR penetration at the ground of 1, 3, and 5% gives CR = 0.69, 0.52, and 0.45 respectively. I use CR values of 0.5 for conifer forest, 0.6 for broadleaved forest, and 0.7 for short vegetation covers. CR is also used to calculate the extinction of net radiation, though this is theoretically incorrect. [see PET-SRSC] [see SUN-AVAILEN]
    # CS (Canopy parameter) - ratio of projected stem area index (SAI) to HEIGHT when DENSEF = 1. SAI = CS * HEIGHT * DENSEF. For a total stem area index of closed canopy forests of roughly 2 for a 20 m height (Federer and Lash 1978), and cylindrical "stems" then CS is 2 / (20 π) or 0.035. A high value of CS and a small MAXHT might be used to simulate logging slash in a clearcut. If DENSEF is not reduced, CS should be reduced for a sparse canopy. [see PET-CANOPY]
    # FXYLEM (Canopy parameter) - fraction of plant resistance that is in the xylem. FXYLEM is the fraction of the internal plant resistance to water flow that is in the xylem, which is considered to be all above ground.The remaining resistance is considered to be in the root cortex and is therefore distributed among soil layers containing roots. Increasing FXYLEM reduces the dependence of layer uptake on root density in the layer and thus makes transpiration uptake more uniform with depth. FXYLEM probably should be zero for short canopies and increase linearly with height to 0.5 for forests of MAXHT = 25 m (Hunt et al. 1991). It can vary from 0 to 1. When FXYLEM = 1 it is reset to 0.99. [see EVP-PLNTRES]
    # GLMAX (Canopy parameter) - maximum leaf conductance, cm/s. In code GLMAX is in m/s and GLMAXC is in cm/s. GLMAX is the maximum leaf conductance when stomates are fully open.It is the total conductance of all sides of a leaf or needle, based on projected leaf area. This is an important parameter controlling potential transpiration. Values for all plant types should be in the range 0.2 to 2.0 cm/s. [see PET-SRSC]. This conductance is reduced by low light, low or high temperature, and high vapor pressure deficit in calculating canopy resistance (RSC) and thus potential transpiration (PTRAN). Values by cover type in the b90v44data.zip Canopy parameter files are from Körner (1994). He concluded that GLMAX does not vary among forest types, even though it differs among individual species. Guidance in selecting GLMAX for specific species can be found in Hinckley et al. (1978) and Körner (1979), or this can be a fitting parameter to give appropriate transpiration. [see PET-SRSC]
    # GLMIN (Canopy parameter) - minimum leaf conductance, cm/s. In code, GLMIN is in m/s and GLMINC is in cm/s. GLMIN is the average nighttime leaf conductance, or the conductance when stomates are closed. This value is used for the whole day when mean daily temperature is below TL or above TH. If GLMIN is set to 0 a value of 0.00001 is used to avoid possible zero divide. GLMIN is fixed at 0.0003 m/s. [see PET-SRSC]
    # KSNVP (Canopy parameter) - reduction factor between 0.05 and 1 to reduce snow evaporation (SNVP), dimensionless. It is needed for tall canopies because the Shuttleworth-Gurney aerodynamic resistances as corrected for SAI are too low and thus overestimate SNVP from forests. Set KSNVP to 1 for short canopies, and in the absence of better information, use 0.3 for tall forest and intermediate values for intermediate heights. To turn off SNVP set KSNVP to zero. [see SNO-SNOVAP]
    # LWIDTH (Canopy parameter) - average leaf width, m. LWIDTH is the average leaf width (generally its second smallest dimension) used to determine the leaf boundary resistance RAC. [see PET-SWGRA]
    # PSICR (Canopy parameter) - minimum plant leaf water potential, MPa. PSICR is the critical leaf water potential at which stomates close. BROOK90 assumes that transpiration is limited by potential transpiration (PTRAN) until water uptake at a plant water potential of PSICR is less than PTRAN. PSICR can be considered as the water potential at the turgor-loss point. PSICR varies from -1.5 to -3.0 MPa for most species and is quite species dependent (Hinckley et al. 1978). This parameter is best selected from knowledge of the water potential - diffusion resistance relation for the species involved. [see KPT-SOILPAR] [see EVP-TBYLAYER]
    # ALB and ALBSN (Canopy parameters) - albedo or surface reflectivity without and with snow on the ground, respectively. Used to calculate net radiation from solar radiation [see SUN-AVAILEN]. Except for choosing which of these two values to use, albedo does not vary with time in BROOK90. ALB generally ranges from 0.1 for bare soil and water to 0.3 for some vegetation. ALBSN ranges from 0.1 for some forests to 0.9 for bare fresh snow. In BROOK90 ALBSN only affects interception and transpiration, not snow evaporation and snow energy balance. [see SUN-AVAILEN]
    # TH, T1, T2, and TL - input parameters- temperatures controlling when stomates are closed, °C. TL, T1, T2, and TH express the reduction of leaf conductance caused by sub-optimal temperature. When the mean air temperature for the day (TA) is less than TL or greater than TH, stomates are closed and GLMIN applies. When TA is between T1 and T2, there is no stomatal closure induced by suboptimal temperature. BROOK90 uses a parabolic interpolation when TA is in the TL to T1 and T2 to TH ranges. If TL = T1 and TH = T2 the temperature response is square. When TL and T1 are very low and T2 and TH are very high there is no temperature effect on leaf conductance. There is very little data on temperature response by species as it is difficult to separate from the vapor pressure response (Körner 1994). Values are TL = 0°C, T1 = 10°C, T2 = 30°C, and TH = 40°C are reasonable choices. [see PET-SRSC]
    # CVPD (Fixed parameter) - vapor pressure deficit at which stomatal conductance is halved, kPa. Values of cD between 0.5 and 2 kPa are generally found in the literature. In BROOK90 CVPD is fixed at 2 kPa for all cover types. A very large value will remove vapor pressure dependence. Values less than 1 kPa produce too much sensitivity. After GLMAX, this parameter probably has the next largest effect on potential transpiration. Unfortunately there is very little data on how it varies among species (Körner 1994). [see PET-SRSC]
    # LPC (Fixed parameter) - minimum LAI defining a closed canopy, dimensionless. LPC is the projected leaf area index above which the canopy is considered closed in the Shuttleworth and Wallace (1985) equations. LPC is fixed at 4.0. [see PET-ROUGH]
    # NN (Fixed parameter) - wind/diffusivity extinction coefficient, dimensionless. NN is the canopy extinction coefficient for wind and eddy diffusivity. NN is fixed at 2.5 following Shuttleworth and Gurney (1990). Federer et al. (1995) show that PE is insensitive to the value of n, but they did not test sparse canopies in a wet climate. [see PET-SWGRA]
    # R5 (Fixed parameter) - solar radiation at which stomatal conductance is half of its value at RM, W/m2. Values in the range of 50 to 200 W/m2 are likely for many species. BROOK90 fixes R5 = 100 W/m2 as the default for all cover types because insufficient data are available to describe R.5 by cover type (Körner 1994). [see PET-SRSC]
    # RHOTP (Fixed parameter) - ratio of total leaf area to projected area, dimensionless. RHOTP is always 2 for broadleaves, and ranges from 2 for flat needles to π for cylindrical needles. RHOTP is fixed at 2. The difference between 2 and 3 is generally negligible. [see PET-SWGRA] [see PET-SRSC]
    # RM (Fixed parameter) - nominal maximum solar radiation possible on a leaf, W/m2, as used in calculating leaf conductance. See also R5. RM is fixed at 1000 W/m2. [see PET-SRSC]
    # CINTRL and CINTRS (Fixed parameters) - maximum interception storage of rain per unit LAI and SAI, respectively, mm. This storage is only removed by evaporation, so does not include water that drips off. CINTRL and CINTRS are both fixed at 0.15 mm. Studies of rain interception in mature forests generally indicate a capacity of 1 or 2 mm. With CINTRL and CINTRS set at 0.15, an LAI of 6 and an SAI of 0.7 gives a capacity of 1.0 mm. LAI is on a projected area basis, so the assumption is that only the projected area is wetted. [see EVP-INTER]
    # CINTSL and CINTSS (Fixed parameters) - maximum interception storage of snow water per unit LAI and SAI, respectively, mm. This storage is only removed by evaporation, so does not include water that drips or falls off. CINTSL and CINTSS are both fixed at 0.6. With LAI = 6 and SAI = 0.7, maximum snow water interception capacity is 4.0 mm corresponding to the value used by Federer and Lash (1978b). This capacity is slightly lower than the 5 and 7.5 mm that Leaf and Brink (1973) used for lodgepole pine and spruce-fir, but BROOK90 assumes that all of this will evaporate while they do not. When LAI = 0 and SAI = 0.7 the capacity is 0.4 mm. LAI is on a projected area basis, so the assumption is that only the projected area is wetted. [see EVP-INTER]
    # FRINTLAI and FRINTSAI (Fixed parameters) - intercepted fraction of rain per unit LAI and per unit SAI respectively, dimensionless. See also FSINTLAI. FRINTLAI and FRINTSAI are both fixed at 0.06. For LAI = 6 and SAI = 0.7 these values give a rain catch rate of 40% of the rainfall rate. For leafless deciduous forest with LAI = 0 and SAI = 0.7 the rain catch rate is 4%. To turn off RINT, set both FRINTLAI and FRINTSAI to zero. [see EVP-INTER]
    # FSINTLAI and FSINTSAI (Fixed parameters) - intercepted fraction of snow per unit LAI and per unit SAI respectively, dimensionless. See also FRINTLAI. FSINTLAI and FSINTSAI are both fixed at 0.04. For LAI = 6 and SAI = 0.7 these values catch snow at 27%. For leafless deciduous forest with LAI = 0 and SAI = 0.7 the snowfall catch rate is 3%. To turn off SINT, set both FSINTLAI and FSINTSAI to zero. [see EVP-INTER]

    ## Soil vegetation (roots)
    NOOUTF   = 1 == continuous_SPAC.params[:NOOUTF] # flag to prevent outflow from roots (hydraulic redistribution), (0/1)
    p_RTRAD  = continuous_SPAC.params[:RTRAD] # average root radius, mm (BROOK90: RTRAD is fixed at 0.35 mm)
    p_MXRTLN = continuous_SPAC.params[:MXRTLN] # (Canopy parameter), maximum length of fine roots per unit ground area, m/m2.
    # Documentation from ecoshift:
    # NOOUTF (Fixed parameter) - 0 to allow outflow from roots, 1 for no outflow. NOOUTF is a switch that when set to 1 prevents outflow from the plant roots to the soil when soil is dry. NOOUTF = 0 allows such outflow, so water can move from wet soil layers to dry soil layers through the roots. [see EVP-TBYLAYER]
    # RTRAD (Fixed parameter) - average root radius, mm. RTRAD is the average radius of the fine or water-absorbing roots. It is only relevant to transpiration from dry soil. RTRAD is fixed at 0.35 mm. [see EVP-PLNTRES]
    # MXRTLN Total root length per unit area (RTLEN) is MXRTLN * RELHT * DENSEF. MXRTLN is used to calculate rhizosphere resistance and is only important when soil is dry or roots are sparse. Values of MXRTLN are not frequent in the literature, especially for forests. Newman (1974) reported a range of 1700 to 11000 m/m2 for 5 woody plants. Safford (1974) found fine root masses of 1200 g/m2 for northern hardwoods, and Safford and Bell (1972) found 700 g/m2 for white spruce; with a mean diameter of 0.7 mm and density of 0.5 g/cm3, these become 6200 and 3600 m/m2. To turn off TRAN set MXRTLN to zero. [see PET-CANOPY] [see EVP-PLNTRES]

    ## Soil discretization
    NLAYER   = soil_discr["NLAYER"] # Number of soil layers used
    p_THICK  = soil_discr["THICK"]  # (Soil parameter),  layer thicknesses, mm
    # Documentation from ecoshift:
    # NLAYER (Soil parameter) - number of soil layers to be used, dimensionless. NLAYER is the number of soil layers to be used in the model run. It can vary from 1 to ML. Run time is more or less proportional to NLAYER. Soil parameter values for layers greater than NLAYER can be 0.
    # THICK(1 To ML) (Soil parameter) - layer thicknesses, mm. THICK is the vertical thickness of each soil layer. Each layer can have a different thickness, but the number of iterations goes up as the thickness of any layer goes down. THICK should probably not be less than 50 mm unless run time is not important. [see EVP-PLNTRES] [see KPT] [see WAT-VERT]

    ## Soil hydraulics
    p_RSSA   = continuous_SPAC.params[:RSSA] # Soil evaporation resistance (RSS) at field capacity, s/m (BROOK90: RSSA is fixed at 500 s/m following Shuttleworth and Gurney (1990))
    p_RSSB   = continuous_SPAC.params[:RSSB] # Exponent in relation of soil evaporation resistance (RSS) to soil water potential (PSIM) in the top layer, dimensionless, (BROOK90: RSSB is fixed at 1.0, which makes RSS directly proportional to PSIM)
    # TODO(bernharf): get rid of this and simply use the vector of AbstractSoilHydraulicParams...
    FLAG_MualVanGen = typeof(continuous_SPAC.soil_horizons[1,:shp]) == LWFBrook90.KPT.MualemVanGenuchtenSHP # 0 for Clapp-Hornberger; 1 for Mualem-van Genuchten
    if FLAG_MualVanGen == 0
        p_soil = LWFBrook90.KPT.KPT_SOILPAR_Ch1d(;
            p_THICK = p_THICK,
            p_STONEF = [shp.p_STONEF for shp in soil_discr["SHP"]], # stone volume fraction in each soil layer, dimensionless
            p_THSAT  = [shp.p_THSAT  for shp in soil_discr["SHP"]], # THETA at saturation, m3/m3
            p_THETAF = [shp.p_THETAF for shp in soil_discr["SHP"]], # volumetric water content at "field capacity" corresponding to KF and PSIF for soil layer, m3/m3
            p_KF     = [shp.p_KF     for shp in soil_discr["SHP"]], # hydraulic conductivity at field capacity corresponding to THETAF and PSIF for a soil layer, mm/d
            p_PSIF   = [shp.p_PSIF   for shp in soil_discr["SHP"]], # matric potential at "field capacity" corresponding to KF and THETAF for a soil layer, kPa
            p_BEXP   = [shp.p_BEXP   for shp in soil_discr["SHP"]], # Clapp-Hornberger exponent for ψ-θ relation
            p_WETINF = [shp.p_WETINF for shp in soil_discr["SHP"]]) # wetness at dry end of near-saturation range for a soil layer, dimensionless
    elseif FLAG_MualVanGen == 1
        # ..it's a sin......
        # Hard default of 2 [mm d-1] for saturated hydraulic conductivity at field capacity K(θ_fc)

        # Instantiate soil parameters
        hardcodedKθfc = 2. .+ 0. * similar(p_THICK) #PAR[:,"K(θ_fc)"] .= 2.
        p_soil = LWFBrook90.KPT.KPT_SOILPAR_Mvg1d(;
            p_THICK  = p_THICK,
            p_STONEF = [shp.p_STONEF for shp in soil_discr["SHP"]],
            p_THSAT  = [shp.p_THSAT  for shp in soil_discr["SHP"]],
            p_Kθfc   = 2. .+ 0. * similar(p_THICK), #hardcodedKθfc,
            p_KSAT   = [shp.p_KSAT   for shp in soil_discr["SHP"]],
            p_MvGα   = [shp.p_MvGα   for shp in soil_discr["SHP"]],
            p_MvGn   = [shp.p_MvGn   for shp in soil_discr["SHP"]],
            p_MvGl   = [shp.p_MvGl   for shp in soil_discr["SHP"]],
            p_θr     = [shp.p_θr     for shp in soil_discr["SHP"]])
    else
        error("Unsupported FLAG_MualVanGen: $FLAG_MualVanGen")
    end

    # TODO(bernhard): treat following note:
    # NOTE(bernhard) the difference between p_PSICR and p_PsiCrit:
    # p_PSICR (Brook90): PSICR (Canopy parameter) - minimum plant leaf water
    #    potential, MPa. PSICR is the critical leaf water potential at which stomates
    #    close. BROOK90 assumes that transpiration is limited by potential
    #    transpiration (PTRAN) until water uptake at a plant water potential of PSICR
    #    is less than PTRAN. PSICR can be considered as the water potential at the
    #    turgor-loss point. PSICR varies from -1.5 to -3.0 MPa for most species and is
    #    quite species dependent (Hinckley et al. 1978). This parameter is best
    #    selected from knowledge of the water potential - diffusion resistance relation
    #    for the species involved. [see KPT-SOILPAR] [see EVP-TBYLAYER]
    # p_PsiCrit (LWFBrook90):
    #    Definition: minimum soil matric potential to allow water supply for
    #    evapotranspiration , Hammel, 2001 (p_PsiCrit = f(ThCrit) = FPSIM(ThCrit))
    #    use: if (PsiM < PsiCrit) TRANI = 0 # no transpiration
    #    use: if (PsiM < PsiCrit) SLVP = 0  # no soil evaporation # TODO(bernhard): this seems incorrect

    # Documentation from ecoshift:
    # RSSA           (Fixed parameter) - soil evaporation resistance (RSS) at field capacity, s/m. RSSA is the soil surface resistance to soil evaporation at field capacity as defined by PSIF. There is essentially no information on how this parameter varies with the characteristics of the surface soil layer, particularly for forests. RSSA is fixed at 500 s/m following Shuttleworth and Gurney (1990). To eliminate SLVP set RSSA to zero. [see PET-FRSS]
    # RSSB           (Fixed parameter) - exponent in relation of soil evaporation resistance (RSS) to soil water potential (PSIM) in the top layer, dimensionless. There is essentially no information on how this parameter varies. RSSB is fixed at 1.0, which makes RSS directly proportional to PSIM. Setting RSSB = 0 makes the soil evaporation resistance always equal to RSSA. [see PET-FRSS]
    # STONEF(1 To ML) (Soil parameter) - stone volume fraction in each soil layer, dimensionless. STONEF is the fraction of the total layer that is coarse fragments, which are assumed to neither store nor conduct water. STONEF reduces water storage capacity and both downslope and vertical matric flow. STONEF is usually estimated in soil surveys. It can, of course, be set to zero, but for many forest soils, values of 0.1 to 0.4 are more appropriate. [see EVP-PLNTRES] [see KPT] [see WAT-VERT] [see WAT-DSLOP]
    # THSAT(1 To ML)  (Soil parameter) - THETA at saturation, m3/m3. THSAT is the matrix porosity (volume fraction) by layer and is also the water content by volume at saturation. THSAT is only important when the soil is nearly saturated due to impeded vertical flow. [see KPT]
    # THETAF(1 To ML) (Soil parameter) - volumetric water content at "field capacity" corresponding to KF and PSIF for soil layer, m3/m3. THETAF replaces the more usual value of saturated water content in the Clapp and Hornberger (1978) equations. The concept of field capacity is not used for matric flow of soil water, but is used to calculate SLVP, SRFL, BYFL, and ADEF. THETAF is normally obtained for the given PSIF from θ - ψ curves that are obtained from laboratory analysis of soil cores. [see KPT]
    # KF(1 To ML)     (Soil parameter) - hydraulic conductivity at field capacity corresponding to THETAF and PSIF for a soil layer, mm/d. KF replaces the more usual but conceptually questionable value of saturated hydraulic conductivity. The concept of field capacity is not used for matric flow of soil water, but is used to calculate SLVP, SRFL, BYFL, and ADEF. [see KPT for much more on KF]
    # PSIF(1 To ML)   (Soil parameter) - matric potential at "field capacity" corresponding to KF and THETAF for a soil layer, kPa. PSIF replaces the more usual but conceptually questionable value of air entry value in the Clapp and Hornberger (1978) equations. The concept of field capacity is not used for matric flow of soil water, but is used to calculate SLVP, SRFL, BYFL, and ADEF. [see KPT for much more on PSIF] [see PET-FRSS]
    # BEXP(1 To ML)   (Soil parameter) - exponent for ψ-θ relation. BEXP is the negative slope of the log ψ - log θ relationship, or the exponent in the Brooks and Corey equation as given by Clapp and Hornberger (1978). This parameter is usually b in the literature. Values of BEXP above 11.5 will not work unless WETINF is set higher than 0.92. [see KPT]
    # WETINF(1 To ML) (Soil parameter) - wetness at dry end of near-saturation range for a soil layer, dimensionless. WETINF is the wetness (saturation fraction) at the inflection point in the Clapp-Hornberger (1978) equation. This value is the lower limit of the parabolic approach to saturation. BROOK90 changes WETINF if it is outside the range from BEXP/(1+BEXP) to 0.999. [see KPT]

    ## FLOW Infiltration, groundwater, and overland flow
    ### Infiltration (incl. preferential flow)
    p_IMPERV = continuous_SPAC.params[:IMPERV] # (Flow parameter), fraction of impervious surface area generating surface or source area flow (SRFL), dimensionless
    p_INFEXP = continuous_SPAC.params[:INFEXP] # (Flow parameter), infiltration exponent that determines the distribution of infiltrated water with depth, dimensionless (from 0 to >1; 0 = all infiltration to top soil layer, 1 = uniform distribution down to ILAYER, >1 = more water in lower layers closer to ILAYER)
    p_ILAYER = soil_discr["ILAYER"] # (Flow parameter), number of layers over which infiltration is distributed
    p_QLAYER = soil_discr["QLAYER"] # (Flow parameter), number of soil layers for SRFL
    p_INFRAC = LWFBrook90.WAT.INFPAR(p_INFEXP, p_ILAYER, p_soil, NLAYER) # fraction of (preferential) infiltration to each layer
    # TODO(bernhard):switch to ILAYAER and QLAYER to IDEPTH_m and QDEPTH_m, which are independent of soil discretization.

    ### Flow generation
    p_BYPAR  = continuous_SPAC.params[:BYPAR]  # (Flow parameter), flag to activate bypass flow (BYFL), (0/1)
    p_DRAIN  = continuous_SPAC.params[:DRAIN] # (Flow parameter), continuous flag to activate drainge VRFLI(n), (between 0 and 1; 1 = gravity drainage, 0 = no drainage)
    p_DSLOPE = continuous_SPAC.params[:DSLOPE] # (Flow parameter), hillslope angle for downslope matric flow (DSFL), degrees
    p_GSC    = continuous_SPAC.params[:GSC] # (Flow parameter), fraction of groundwater storage (GWAT), that is transferred to groundwater flow (GWFL) and deep seepage (SEEP) each day, d-1
    p_GSP    = continuous_SPAC.params[:GSP] # (Flow parameter), fraction of groundwater discharge produced by GSC that goes to deep seepage (SEEP) and is not added to streamflow (FLOW), dimensionless
    p_LENGTH_SLOPE = continuous_SPAC.params[:LENGTH_SLOPE] # (Flow parameter), slope length for downslope flow (DSFL), m
    p_QFFC   = continuous_SPAC.params[:QFFC] # (Flow parameter), quick flow fraction for SRFL and BYFL at THETAF, dimensionless
    p_QFPAR  = continuous_SPAC.params[:QFPAR] # (Flow parameter), raction of the water content between field capacity (THETAF) and saturation (THSAT) at which the quick flow fraction is 1, dimensionless

    # source area parameters SRFPAR()
    p_SWATQX = sum(p_soil.p_SWATMAX[1:p_QLAYER]) # maximum water storage for layers 1 through QLAYER, mm
    p_SWATQF = sum(
        p_soil.p_THETAF[1:p_QLAYER] .*
        p_soil.p_THICK[1:p_QLAYER] .*
        (1 .- p_soil.p_STONEF[1:p_QLAYER]))     # water storage at field capacity for layers 1 through QLAYER, mm

    # Documentation from ecoshift:
    # BYPAR (Flow parameter) - either 0 to prevent bypass flow (BYFL), or 1 to allow BYFL . BYFL is zero when BYPAR = 0, unless the surface layer becomes saturated. When BYPAR = 1, a fraction of infiltration to each layer is immediately routed to bypass flow to simulate downslope macropore or pipe flow of new water. The fraction depends on QFFC and QFPAR. These are the same parameters used to determine source area flow (SRFL), so simulation with both SRFL and BYFL (QDEPTH_m > 0 and BYPAR = 1) is discouraged. The difference is that SRFL depends on the total water content down to QDEPTH_m whereas BYFL depends on water content in each layer down to IDEPTH_m. When NLAYER = 1, SRFL and BYFL are identical. [see WAT-BYFLFR]
    # DRAIN (Flow parameter) - multiplier between 0 and 1 of drainage from the lowest soil layer, VRFLI(n), for drainage to groundwater, dimensionless. DRAIN = 1 produces vertical drainage under gravity gradient. DRAIN = 0 prevents drainage from the bottom of the soil column. Values between 0 and 1 can also be used, especially to control DSFL. [see WAT-VERT] [see WAT-DSLOP]
    # DSLOPE and DSLOPED (Flow parameter) - hillslope angle for downslope matric flow (DSFL), degrees. In code DSLOPED is in degrees, DSLOPE is in radians. Because downslope flow is overparameterized, arbitrarily setting DSLOPE to 10° is satisfactory for DSFL from the bottom soil layer(s). When either DSLOPE or LENGTH_SLOPE is 0 there is no DSFL, and the other parameter is ignored. [see WAT-DSLOP]
    # GSC (Flow parameter) - fraction of groundwater storage (GWAT), that is transferred to groundwater flow (GWFL) and deep seepage (SEEP) each day, d-1. Where groundwater is being simulated, GSC should be some fraction like 0.1 d-1 or less. If GSC = 0 there is no groundwater storage and all vertical drainage from the soil profile becomes seepage or streamflow directly. See also GSP. [see WAT-GWATER]
    # GSP (Flow parameter) - fraction of groundwater discharge produced by GSC that goes to deep seepage (SEEP) and is not added to streamflow (FLOW), dimensionless. If GSC = 0, GSP applies to vertical drainage from the bottom soil layer. To eliminate SEEP set GSP to zero. [see WAT-GWATER]
    # IDEPTH_m (Flow parameter) - depth over which infiltration is distributed, mm. IDEPTH_m determines the number of soil layers over which infiltration is distributed when INFEXP is greater than 0. It should correspond to the depth of vertical macropores. IDEPTH_m does not need to correspond to the bottom of a soil layer; it is converted into the number of soil layers most closely corresponding to IDEPTH_m. [see WAT-INFPAR].
    # IMPERV (Flow parameter) - fraction of the soil surface that is impermeable and always routes water reaching it directly to streamflow as SRFL. For a watershed, IMPERV represents at least the area of the stream channel; an appropriate value is 0.01. To turn off SLFL set IMPERV = 1. To turn off SRFL set both IMPERV and QDEPTH_m to zero. [see WAT-SRFLFR]
    # INFEXP (Flow parameter) - infiltration exponent that determines the distribution of infiltrated water with depth, dimensionless. When INFEXP = 0, all infiltration goes to the top soil layer and a classic top-down wetting front is produced. Increasing INFEXP corresponds to increasing macropore-assisted infiltration, and produces an exponential depth distribution of infiltrated water down through the layer whose lower depth most closely corresponds to IDEPTH_m. With INFEXP = 1, infiltrated water is distributed uniformly down to the IDEPTH_m layer. Values above 1 put more water into lower layers than into upper layers. [see WAT-INFPAR] [see WAT-BYFLFR]
    # LENGTH_SLOPE (Flow parameter) - slope length for downslope flow (DSFL), m. LENGTH_SLOPE is conceptually the hillslope length in m as horizontal or map distance from ridge to channel. But in practice it is a fitted value to produce the desired amount of DSFL, which is roughly inversely proportional to LENGTH_SLOPE. Downslope flow from the bottom soil layer(s) is overparameterized, so LENGTH_SLOPE can be set to 10 m and DSFL can be varied by changing DRAIN. When either DSLOPE or LENGTH_SLOPE are 0 there is no downslope flow, and the other parameter is ignored. [see WAT-DSLOP]
    # QDEPTH_m (Flow parameter) - soil depth for SRFL calculation, mm. QDEPTH_m determines the number of soil layers over which wetness is calculated to determine source area fraction and SRFL. QDEPTH_m does not need to correspond to the bottom of a soil layer, but it is converted into the number of soil layers most closely corresponding to QDEPTH_m. When QDEPTH_m equals or exceeds the depth of NLAYER then all NLAYERs are used. When QDEPTH_m = 0, the source area fraction is equal to IMPERV. Smaller QDEPTH_m means a larger contrast between wet and dry conditions, and thus more responsiveness of source area fraction. See also QFPAR and QFFC. To turn off SRFL set both IMPERV and QDEPTH_m to zero. [see WAT-BYFLFR] [see WAT-SRFLFR] [see WAT-SRFPAR]
    # QFFC (Flow parameter) - quick flow fraction for SRFL and BYFL at THETAF, dimensionless. QFFC is used for both SRFL and BYFL, so normally only one or the other should be simulated. For Hubbard Brook Watershed 6, QFFC = 0.2 and QFPAR = 0.3 fit storm hydrographs well using SRFL; these values give a generally high stormflow response. Decreasing QFFC decreases SRFL and BYFL proportionally at all soil water contents. See also QFPAR, and QDEPTH_m. QFFC must be greater than 0.0001. When both BYPAR and QDEPTH_m = 0, QFFC is ignored. When QFFC = 1 there is never any infiltration into the top soil layer. [see WAT-BYFLFR] [see WAT-SRFLFR]
    # QFPAR (Flow parameter) - fraction of the water content between field capacity (THETAF) and saturation (THSAT) at which the quick flow fraction is 1, dimensionless. QFPAR is used for both SRFL and BYFL, so normally only one or the other should be simulated. When QFPAR = 0 (<= 0.01), quick flow operates like a field-capacity bucket; all water input above field capacity becomes BYFL or SRFL. For Hubbard Brook Watershed 6, QFFC = 0.2 and QFPAR = 0.3 fit storm hydrographs well using SRFL; these values give a generally high stormflow response. Increasing QFPAR increases quickflow from soil dryer than THETAF and decreases it from soil wetter than THETAF. Values > 1 are allowed. For SRFL, the average water content of layers down through QDEPTH_m controls the flow rate. For BYFL, the water content of each layer controls the flow rate from that layer. When both BYPAR and QDEPTH_m = 0, QFPAR is ignored. [see WAT-BYFLFR] [see WAT-SRFLFR]
    # INFRAC(1 To ML) - fraction of infiltration to each layer, calculated parameter. [see WAT-INFPAR] [see WAT-INFLOW]
    # SWATQF - water storage at field capacity for layers 1 to ,mm, calculated parameter. [see WAT-SRFLFR] [see WAT-SRFPAR]
    # SWATQX - maximum water storage for layers 1 to QLAYER, mm, calculated parameter. [see WAT-SRFLFR] [see WAT-SRFPAR]

    ########
    # 2) Define parameters for differential equation:
    # 2a) Constant parameters

    # p_cst_1 and p_cst_2 for both RHS and CallBack in DiffEq.jl
    p_cst_1 = p_soil
    p_cst_2 = (NLAYER, FLAG_MualVanGen, continuous_SPAC.solver_options.compute_intermediate_quantities, false, # Reset is hardcoded as false
        p_DTP, p_NPINT,

        # FOR MSBITERATE:
        p_QLAYER, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC, p_IMPERV,
        p_LENGTH_SLOPE, p_DSLOPE, LWFBrook90.CONSTANTS.p_RHOWG, p_DPSIMAX,
        p_DRAIN, p_DTIMAX, p_INFRAC, p_DSWMAX,
        p_GSC, p_GSP,

        p_BYPAR)

    # p_cst_3 only for CallBack in DiffEq.jl
    p_cst_3 = (p_LAT, p_ESLOPE, p_L1, p_L2,
        p_SNODEN, p_MXRTLN, p_MXKPL, p_CS,
        p_Z0S, p_Z0G,
        p_ZMINH, p_CZS, p_CZR, p_HS, p_HR, p_LPC,
        p_RTRAD, p_FXYLEM,
        p_WNDRAT, p_FETCH, p_Z0W, p_ZW,
        p_RSTEMP,
        LWFBrook90.CONSTANTS.p_CVICE,
        p_LWIDTH, p_RHOTP, p_NN, p_KSNVP,
        p_ALBSN, p_ALB,
        p_RSSA, p_RSSB,
        p_CCFAC, p_MELFAC, p_LAIMLT, p_SAIMLT,

        LWFBrook90.CONSTANTS.p_WTOMJ, p_C1, p_C2, p_C3, p_CR,
        p_GLMIN, p_GLMAX, p_R5, p_CVPD, p_RM, p_TL, p_T1, p_T2, p_TH,
        p_PSICR, NOOUTF,

        # for MSBPREINT:
        p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS,
        p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS,
        p_DURATN, p_MAXLQF, p_GRDMLT,

        # for isotope mixing:
        p_VXYLEM, p_DISPERSIVITY)


    p_cst_4 = (
        FLAG_MualVanGen,
        continuous_SPAC.solver_options.compute_intermediate_quantities,
        continuous_SPAC.solver_options.simulate_isotopes,
        row_idx_scalars = [],
        row_idx_SWATI   = [],
        row_idx_RWU     = [],
        row_idx_accum   = [],
        names_accum     = [],
        col_idx_d18O    = [],
        col_idx_d2H     = []) # Placeholders will be overwritten when defining u0
    p_cst = (p_cst_1, p_cst_2, p_cst_3, p_cst_4)

    # 2b) Time varying parameters (e.g. meteorological forcings)
    p_fT = (p_DOY          = (t) -> LWFBrook90.p_DOY(t,    continuous_SPAC.reference_date),
            p_MONTHN       = (t) -> LWFBrook90.p_MONTHN(t, continuous_SPAC.reference_date),
            p_GLOBRAD      = continuous_SPAC.meteo_forcing.p_GLOBRAD,
            p_TMAX         = continuous_SPAC.meteo_forcing.p_TMAX,
            p_TMIN         = continuous_SPAC.meteo_forcing.p_TMIN,
            p_VAPPRES      = continuous_SPAC.meteo_forcing.p_VAPPRES,
            p_WIND         = continuous_SPAC.meteo_forcing.p_WIND,
            p_PREC         = continuous_SPAC.meteo_forcing.p_PREC,

            p_DENSEF       = continuous_SPAC.canopy_evolution.p_DENSEF, # canopy density multiplier between 0.05 and 1, dimensionless
            p_HEIGHT       = continuous_SPAC.canopy_evolution.p_HEIGHT,
            p_LAI          = continuous_SPAC.canopy_evolution.p_LAI,
            p_SAI          = continuous_SPAC.canopy_evolution.p_SAI,
            p_AGE          = continuous_SPAC.canopy_evolution.p_AGE,
            p_RELDEN       = p_RELDEN,

            p_d18OPREC     = continuous_SPAC.meteo_iso_forcing.p_d18OPREC,
            p_d2HPREC      = continuous_SPAC.meteo_iso_forcing.p_d2HPREC,
            REFERENCE_DATE = continuous_SPAC.reference_date)
    # Documentation from ecoshift:
    # DENSEF (Fixed parameter) - canopy density multiplier between 0.05 and 1, dimensionless. DENSEF is normally 1; it should be reduced below this ONLY to simulate thinning of the existing canopy by cutting. It multiplies MAXLAI, CS, MXRTLN, and MXKPL and thus proportionally reduces LAI, SAI, and RTLEN, and increases RPLANT. However it does NOT reduce canopy HEIGHT and thus will give erroneous aerodynamic resistances if it is less than about 0.05. It should NOT be set to 0 to simulate a clearcut. [see PET-CANOPY]


    # 2c) Time varying "parameters" (depending on state variables)
    #     These need to be exchanged between CallBack and RHS in DiffEq.jl which is why they
    #     can temporarily be saved in the parameter vector to avoid computing them twice

    # Initialize placeholder for parameters that depend on solution and are computed
    p_fu = ([NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN],
            fill(NaN, NLAYER),     # see Localizing variables helps to ensure type stability. under https://nextjournal.com/sosiris-de/ode-diffeq?change-id=CkQATVFdWBPaEkpdm6vuto
            fill(NaN, (NLAYER,5)), # for du_NTFLI, aux_du_VRFLI, aux_du_DSFLI, aux_du_INFLI, u_aux_WETNES)
            # TODO(bernhard): should we combine these caches to a matrix of (NLAYER,6)?
            fill(NaN, 2))          # for du_GWFL du_SEEP
    #TODO(bernhard): what are the additional 3x NaNs needed for in isotope code???
    #Earlier it was simply: [NaN, NaN, NaN, NaN, NaN, NaN], fill(NaN, NLAYER)

    # Initialize further caches to hold computed quantities without need to allocate memory
    p_cache = ((
        # chaches for flow equation
        zeros(NLAYER), # u_aux_WETNES
        zeros(NLAYER), # u_aux_PSIM
        zeros(NLAYER), # u_aux_PSITI
        zeros(NLAYER), # u_aux_θ
        zeros(NLAYER), # u_aux_θ_tminus1
        zeros(NLAYER), # p_fu_KK
        zeros(NLAYER), # aux_du_DSFLI
        zeros(NLAYER), # aux_du_VRFLI
        zeros(NLAYER), # aux_du_VRFLI_1st_approx
        zeros(NLAYER), # aux_du_INFLI
        zeros(NLAYER), # aux_du_BYFLI
        zeros(NLAYER), # du_NTFLI
        zeros(NLAYER)), # p_fu_BYFRAC
        # chaches for advection dispersion equation
            # for quantities (all NLAYER long):
            # 7 vectors: θᵏ⁺¹, θᵏ, C_¹⁸Oᵏ⁺¹, C_¹⁸Oᵏ, C_²Hᵏ⁺¹, C_²Hᵏ, q,
            # 9 vectors: Tsoil_K, τw, Λ, D⁰_¹⁸O, D⁰_²H, D_¹⁸O_ᵏ⁺¹, D_²H_ᵏ⁺¹, C_¹⁸O_SLVP, C_²H_SLVP,
            # 8 vectors: diff¹⁸O_upp, diff²H_upp, qCᵢ¹⁸O_upp, qCᵢ²H_upp,
            #            diff¹⁸O_low, diff²H_low, qCᵢ¹⁸O_low, qCᵢ²H_low,
            # 4 vectors: du_Cᵢ¹⁸_SWATI, du_Cᵢ²H_SWATI, du_δ18O_SWATI, du_δ2H_SWATI
        Tuple(zeros(NLAYER) for i=1:28)#,
            # 4 vectors: diff¹⁸O_interfaces, diff²H_interfaces, qCᵢ¹⁸O_interfaces, qCᵢ²H_interfaces,
        # Tuple(zeros(NLAYER+1) for i=1:4)
        )

    # 3) Return different types of parameters as a single object
    return (p_cst, p_fT, p_fu, p_cache)
end




"""
    interpolate_spaceDiscretized_root_density()

Take root density distribution parameters and generates root density in time (t) and space (z, negative downward).
"""

function interpolate_spaceDiscretized_root_density(;
        timepoints,
        p_AGE,
        p_INITRDEP,
        p_INITRLEN,
        p_RGROPER,
        p_tini,
        p_frelden)

    NLAYER = length(p_frelden)

    # 2a Compute time dependent root density parameters
    # Which is a vector quantity that is dependent on time:
    p_RELDEN_2Darray = fill(NaN, length(timepoints), NLAYER)
    for i in eachindex(timepoints)
        p_RELDEN_2Darray[i,:] =
            LWFBrook90.WAT.LWFRootGrowth(;
                final_relden      = p_frelden,
                tini_yrs          = p_tini,
                INITRDEP_m        = p_INITRDEP,
                INITRLEN_m_per_m2 = p_INITRLEN,
                RGROPER_yrs       = p_RGROPER,
                age_yrs           = p_AGE(timepoints[i]),
                NLAYER            = NLAYER)
    end
    p_RELDEN =  extrapolate(interpolate((timepoints, 1:NLAYER), p_RELDEN_2Darray,
                                        (Gridded(Constant{Next}()), NoInterp()), # 1st dimension: ..., 2nd dimension NoInterp()
                                        ), Flat()) # extrapolate flat, alternative: Throw()
    return p_RELDEN
end

"""
    interpolate_meteoveg(
        input_meteoveg::DataFrame,
        input_meteoveg_reference_date::DateTime,
        input_meteoiso::DataFrame, )

Take climate and vegetation parameters in `input_meteoveg` and `input_meteoiso` and generates continuous parameters.
"""
function interpolate_meteoveg(;
    canopy_evolution::DataFrame,
    meteo_forcing::DataFrame,
    meteo_iso_forcing::Union{DataFrame,Nothing},
    p_MAXLAI,
    p_SAI_baseline_,
    p_DENSEF_baseline_,
    p_AGE_baseline_yrs,
    p_HEIGHT_baseline_m)

    @assert meteo_forcing.days == canopy_evolution.days "Discrete DataFrames: `meteo_forcing` and `canopy_evolution` are not available for the same dates."
    # @assert meteo_iso_forcing.days # NOTE: DataFrame `meteo_iso_forcing` can be on a different spacing

    # 2) Interpolate input data in time
    # ### FOR DEVELOPMENT:
    # # Test interpolation based on regular grid:
    # # using Plots
    # # time_range = range(minimum(input_meteoveg.days), maximum(input_meteoveg.days), length=length(input_meteoveg.days))
    # # ts = 0:0.01:365
    # # scatter(input_meteoveg.days[1:365], input_meteoveg.PRECIN[1:365])
    # # plot!(ts, scale(interpolate(input_meteoveg.PRECIN, (BSpline(Constant{Previous}()))), time_range)(ts), label = "PRECIN {Previous}",  xlims=(0,30))
    # # plot!(ts, scale(interpolate(input_meteoveg.PRECIN, (BSpline(Constant{Next}()))),     time_range)(ts), label = "PRECIN {Next}",      xlims=(0,30))
    # # plot!(ts .+ 1,
    # #       scale(interpolate(input_meteoveg.PRECIN, (BSpline(Constant{Next}()))),     time_range)(ts), label = "PRECIN+1 {Previous}",      xlims=(0,30))
    # # plot!(ts, scale(interpolate(input_meteoveg.PRECIN, (BSpline(Constant()))),           time_range)(ts), label = "PRECIN {Nearest}",   xlims=(0,30))
    # # plot!(ts, scale(interpolate(input_meteoveg.PRECIN, (BSpline(Linear()))),             time_range)(ts), label = "PRECIN Linear",   xlims=(0,30))

    # # Test interpolation based on irregular:
    # test_ip_mv = input_meteoveg[[1:10...,15,17:19..., 25:4326...], :]
    # time_range = range(minimum(test_ip_mv.days), maximum(test_ip_mv.days), length=length(test_ip_mv.days))
    # ts = 0:0.01:365
    # scatter(test_ip_mv.days[1:365], test_ip_mv.PRECIN[1:365])
    # # plot!(ts, scale(interpolate(test_ip_mv.PRECIN, (BSpline(Constant{Previous}()))), time_range)(ts), label = "PRECIN {Previous}",  xlims=(0,30))
    # # plot!(ts, scale(interpolate(test_ip_mv.PRECIN, (BSpline(Constant{Next}()))),     time_range)(ts), label = "PRECIN {Next}",      xlims=(0,30))
    # plot!(ts .+ 1,
    #       scale(interpolate(test_ip_mv.PRECIN, (BSpline(Constant{Next}()))),     time_range)(ts), label = "PRECIN+1 {Previous}",      xlims=(0,30))
    # # plot!(ts, scale(interpolate(test_ip_mv.PRECIN, (BSpline(Constant()))),           time_range)(ts), label = "PRECIN {Nearest}",   xlims=(0,30))
    # # plot!(ts, scale(interpolate(test_ip_mv.PRECIN, (BSpline(Linear()))),             time_range)(ts), label = "PRECIN Linear",   xlims=(0,30))
    # plot!(ts, scale(interpolate(test_ip_mv.PRECIN, (BSpline(Linear()))),             time_range)(ts), label = "PRECIN Linear",   xlims=(0,30))
    # plot!(ts, extrapolate(interpolate((test_ip_mv.days, ), test_ip_mv.PRECIN, Gridded(Linear())), Throw())(ts), label = "Gridded(Linear()",xlims=(0,30))
    # plot!(ts, extrapolate(interpolate((test_ip_mv.days, ), test_ip_mv.PRECIN, Gridded(Constant())), Throw())(ts), label = "Gridded(Constant()",xlims=(0,30))
    # plot!(ts, extrapolate(interpolate((test_ip_mv.days, ), test_ip_mv.PRECIN, Gridded(Constant{Next}())), Throw())(ts), label = "Gridded(Constant{Previous}()",xlims=(0,30))

    #     ts = 0:1:4000
    #     scatter(meteo_iso_forcing.days, meteo_iso_forcing.delta18O_permil)
    #     # # plot!(ts, scale(interpolate(meteo_iso_forcing.delta18O_permil, (BSpline(Constant{Previous}()))), time_range)(ts), label = "PRECIN {Previous}",  xlims=(0,30))
    #     # # plot!(ts, scale(interpolate(meteo_iso_forcing.delta18O_permil, (BSpline(Constant{Next}()))),     time_range)(ts), label = "PRECIN {Next}",      xlims=(0,30))
    #     # plot!(ts .+ 1,
    #     #     scale(interpolate(meteo_iso_forcing.delta18O_permil, (BSpline(Constant{Next}()))),     time_range)(ts), label = "PRECIN+1 {Previous}",      xlims=(0,30))
    #     # # plot!(ts, scale(interpolate(meteo_iso_forcing.delta18O_permil, (BSpline(Constant()))),           time_range)(ts), label = "PRECIN {Nearest}",   xlims=(0,30))
    #     # # plot!(ts, scale(interpolate(meteo_iso_forcing.delta18O_permil, (BSpline(Linear()))),             time_range)(ts), label = "PRECIN Linear",   xlims=(0,30))
    #     # plot!(ts, scale(interpolate(meteo_iso_forcing.delta18O_permil, (BSpline(Linear()))),             time_range)(ts), label = "PRECIN Linear",   xlims=(0,30))
    #     plot!(ts, extrapolate(interpolate((meteo_iso_forcing.days, ), meteo_iso_forcing.delta18O_permil, Gridded(Linear())), Throw())(ts), label = "Gridded(Linear()")
    #     plot!(ts, extrapolate(interpolate((meteo_iso_forcing.days, ), meteo_iso_forcing.delta18O_permil, Gridded(Constant())), Throw())(ts), label = "Gridded(Constant()")
    #     plot!(ts, extrapolate(interpolate((meteo_iso_forcing.days, ), meteo_iso_forcing.delta18O_permil, Gridded(Constant{Next}())), Throw())(ts), label = "Gridded(Constant{Previous}()")
    #     xlims!(0, 500)
    # ### END FOR DEVELOPMENT

    # # Interpolate regularly spaced data (`meteo_forcing` and `canopy_evolution`) with Previous.
    # # (e.g. PREC reported for 2010-01-01 is rain that has fallen on that day and rate is applied until 2010-01-02.)
    # # Note that we need to shift it slightly ahead so that p_GLOBRAD(1.0) is equal to p_GLOBRAD(1.1)
    # # Interpolate input data without modification
    # p_GLOBRAD = extrapolate(interpolate((meteo_forcing.days  .- 0.00001, ), meteo_forcing.GLOBRAD,          Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_TMAX    = extrapolate(interpolate((meteo_forcing.days  .- 0.00001, ), meteo_forcing.TMAX,             Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_TMIN    = extrapolate(interpolate((meteo_forcing.days  .- 0.00001, ), meteo_forcing.TMIN,             Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_VAPPRES = extrapolate(interpolate((meteo_forcing.days  .- 0.00001, ), meteo_forcing.VAPPRES,          Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_WIND    = extrapolate(interpolate((meteo_forcing.days  .- 0.00001, ), meteo_forcing.WIND,             Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_PREC    = extrapolate(interpolate((meteo_forcing.days  .- 0.00001, ), meteo_forcing.PRECIN,           Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # # NOTE: PRECIN is already in mm/day from the input data set. No transformation is needed for p_PREC.

    # # Interpolate input data with consideration of baselines
    # p_LAI     = extrapolate(interpolate((canopy_evolution.days  .- 0.00001, ), canopy_evolution.LAI./100 .* p_MAXLAI,              Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_SAI     = extrapolate(interpolate((canopy_evolution.days  .- 0.00001, ), canopy_evolution.SAI./100 .* p_SAI_baseline_,       Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_DENSEF  = extrapolate(interpolate((canopy_evolution.days  .- 0.00001, ), canopy_evolution.DENSEF./100 .* p_DENSEF_baseline_, Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_HEIGHT  = extrapolate(interpolate((canopy_evolution.days  .- 0.00001, ), canopy_evolution.HEIGHT./100 .* p_HEIGHT_baseline_m,Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())

    ### Quickfix:
    @assert 1==length(unique(diff(meteo_forcing.days))) """
        Error: LWFBrook90.jl expects the input data in meteoveg.csv to be regularly spaced in time.
    """
    @assert unique(diff(meteo_forcing.days)) == [1.0] """
        Error: LWFBrook90.jl expects the input data in meteoveg.csv to provided as daily values.
    """
    time_range = range(minimum(meteo_forcing.days), maximum(meteo_forcing.days), length=length(meteo_forcing.days))
    p_GLOBRAD = extrapolate(scale(interpolate(meteo_forcing.GLOBRAD, (BSpline(Constant{Previous}()))), time_range) ,0)
    p_TMAX    = extrapolate(scale(interpolate(meteo_forcing.TMAX,    (BSpline(Constant{Previous}()))), time_range) ,0)
    p_TMIN    = extrapolate(scale(interpolate(meteo_forcing.TMIN,    (BSpline(Constant{Previous}()))), time_range) ,0)
    p_VAPPRES = extrapolate(scale(interpolate(meteo_forcing.VAPPRES, (BSpline(Constant{Previous}()))), time_range) ,0)
    p_WIND    = extrapolate(scale(interpolate(meteo_forcing.WIND,    (BSpline(Constant{Previous}()))), time_range) ,0)
    p_PREC    = extrapolate(scale(interpolate(meteo_forcing.PRECIN,  (BSpline(Constant{Previous}()))), time_range) ,0)

    p_LAI     = extrapolate(scale(interpolate(canopy_evolution.LAI./100 .* p_MAXLAI,              (BSpline(Constant{Previous}()))), time_range) ,0)
    p_SAI     = extrapolate(scale(interpolate(canopy_evolution.SAI./100 .* p_SAI_baseline_,       (BSpline(Constant{Previous}()))), time_range) ,0)
    p_DENSEF  = extrapolate(scale(interpolate(canopy_evolution.DENSEF./100 .* p_DENSEF_baseline_, (BSpline(Constant{Previous}()))), time_range) ,0)
    p_HEIGHT  = extrapolate(scale(interpolate(canopy_evolution.HEIGHT./100 .* p_HEIGHT_baseline_m,(BSpline(Constant{Previous}()))), time_range) ,0)
    ###
    p_AGE     = (t) -> p_AGE_baseline_yrs + t/365
    # p_AGE     = extrapolate(interpolate((canopy_evolution.days  .- 0.00001, ), canopy_evolution.AGE, Gridded(Constant{Next}())), Flat()) #extrapolate flat, alternative: Throw())

    # Note that meteoiso does not need to be regularly spaced:
    # Interpolate irregular data with Next and shift to noon

    if (isnothing(meteo_iso_forcing))
        function p_d18OPREC(t)
            return nothing
        end
        function p_d2HPREC(t)
            return nothing
        end
    else
        ## Quickfix:
        # @assert 1==length(unique(diff(meteo_iso_forcing.days))) """
        #     Error: LWFBrook90.jl expects the input data in meeoiso.csv to be regularly spaced in time.
        # """
        # @assert unique(diff(meteo_iso_forcing.days)) == [1.0] """
        #     Error: LWFBrook90.jl expects the input data in meeoiso.csv to provided as daily values.
        # """
        # time_range_iso = range(minimum(meteo_iso_forcing.days), maximum(meteo_iso_forcing.days), length=length(meteo_iso_forcing.days))
        # p_d18OPREC= extrapolate(scale(interpolate(meteo_iso_forcing.delta18O_permil,    (BSpline(Constant{Previous}()))), time_range_iso) ,0)
        # p_d2HPREC = extrapolate(scale(interpolate(meteo_iso_forcing.delta2H_permil,     (BSpline(Constant{Previous}()))), time_range_iso) ,0)
        ## End Quickfix:
        # Here we shift the days of collection of the cumulative samples to noon. => + 12/24 days
        p_d18OPREC= extrapolate(interpolate((meteo_iso_forcing.days .+ 12/24, ), meteo_iso_forcing.delta18O_permil,    (Gridded(Constant{Next}()))), Flat()) #extrapolate flat, alternative: Throw())
        p_d2HPREC = extrapolate(interpolate((meteo_iso_forcing.days .+ 12/24, ), meteo_iso_forcing.delta2H_permil,     (Gridded(Constant{Next}()))), Flat()) #extrapolate flat, alternative: Throw())
    end

    meteo_forcing_cont = (
        p_days     = time_range,
        p_GLOBRAD  = p_GLOBRAD,
        p_TMAX     = p_TMAX,
        p_TMIN     = p_TMIN,
        p_VAPPRES  = p_VAPPRES,
        p_WIND     = p_WIND,
        p_PREC     = p_PREC)
    meteo_iso_forcing_cont = (
        p_d18OPREC = p_d18OPREC,
        p_d2HPREC  = p_d2HPREC)
    canopy_evolution_cont = (
        p_DENSEF   = p_DENSEF,
        p_HEIGHT   = p_HEIGHT,
        p_LAI      = p_LAI,
        p_SAI      = p_SAI,
        p_AGE      = p_AGE)

    return (meteo_forcing_cont, meteo_iso_forcing_cont, canopy_evolution_cont)
end