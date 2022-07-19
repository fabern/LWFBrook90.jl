using Interpolations: interpolate, extrapolate, NoInterp, Gridded, Constant, Next, Previous, Flat, Throw, scale, BSpline
"""
    define_diff_eq_parameters()

Generate vector p needed for ODE() problem in DiffEq.jl package.

# Arguments
- `NLAYER::...`: TODO argument description.
- `FLAG_MualVanGen::...`: TODO argument description.
- `constant_dt_solver::...`: TODO argument description.
- `NOOUTF::...`: TODO argument description.
- `compute_intermediate_quantities::...`: TODO argument description.
- `input_meteoveg::...`: TODO argument description.
- `input_param::...`: TODO argument description.
- `input_soil::...`: TODO argument description.
- `input_storm_durations::...`: TODO argument description.
- `soil_output_depths`: vector of depths at which state variables should be extractable (negative numeric values [in meter])
"""
function define_LWFB90_p(
    input_meteoveg,
    input_meteoiso,
    input_meteoveg_reference_date,
    input_param,
    input_storm_durations,
    input_soil_horizons,
    input_soil_discretization,
    simOption_FLAG_MualVanGen;
    Reset = false,# TODO(bernhard): remove this
    compute_intermediate_quantities = false,
    simulate_isotopes::Bool = false,
    soil_output_depths = zeros(Float64, 0), # = [],
    constant_dt_solver = 1)

    ########
    ## Discretize soil parameters and interpolate meteo and vegetation parameters
    soil_discr =
        LWFBrook90.discretize_soil_params(
            input_soil_horizons,
            input_soil_discretization,
            soil_output_depths,
            input_param[1,"IDEPTH_m"],
            input_param[1,"QDEPTH_m"],
            input_param[1,"INITRDEP"],
            input_param[1,"RGRORATE"],
            simOption_FLAG_MualVanGen)

    interpolated_meteoveg =
        LWFBrook90.interpolate_meteoveg(
            input_meteoveg,
            input_meteoveg_reference_date,
            input_meteoiso,
            soil_discr["NLAYER"],
            input_param[1,"INITRDEP"],
            input_param[1,"INITRLEN"],
            input_param[1,"RGROPER"],
            soil_discr["tini"],
            soil_discr["frelden"],
            input_param[1,"MAXLAI"],
            input_param[1,"SAI_baseline_"],
            input_param[1,"DENSEF_baseline_"],
            input_param[1,"AGE_baseline_yrs"],
            input_param[1,"HEIGHT_baseline_m"])


    # TODO(bernhard): document input parameters: INITRDEP, INITRLEN, RGROPER, tini, frelden, MAXLAI, HEIGHT_baseline_m
    ########

    ########
    ## Solver algorithm options
    p_DPSIMAX = input_param[1,"DPSIMAX"] # maximum potential difference considered "equal", kPa               (BROOK90: DPSIMAX is fixed at 0.01 kPa)
    p_DSWMAX = input_param[1,"DSWMAX"] # maximum change allowed in soil wetness SWATI, percent of SWATMAX(i) (BROOK90: DSWMAX is fixed at 2 %)
    p_DTIMAX = input_param[1,"DTIMAX"] # maximum iteration time step, d (BROOK90: DTIMAX is fixed at 0.5 d)
    # Documentation from ecoshift:
    # DPSIMAX (Fixed parameter) - maximum potential difference considered equal during soil water integration, kPa. There is no vertical flow between layers whose potentials differ by less than DPSIMAX. This reduces oscillation initiated by flows that are the product of large conductivities and large time steps, but small gradients. The number of iterations used is not at all linearly related to the three iteration parameters DPSIMAX, DSWMAX, and DTIMAX. Selection of values depends on whether the user only wants monthly or daily totals, or is concerned with behaviour at shorter time steps. Generally, faster runs are obtained by using fewer thicker soil layers rather than by using large values of DSWMAX and DPSIMAX. DPSIMAX is fixed at 0.01 kPa. [see WAT-ITER]
    # DSWMAX (Fixed parameter) - maximum change allowed in soil wetness for any layer during an iteration, percent. DSWMAX sets the maximum change in soil wetness or saturation fraction (SWATI / SWATMAX(i)) allowed for any layer in an iteration. See also DPSIMAX. DSWMAX is fixed at 2 %. [see WAT-ITER]
    # DTIMAX (Fixed parameter) - maximum iteration time step, d. DTIMAX is fixed at 0.5 d, which forces at least two iterations per day. This is the largest value that should be used. Much smaller values of DTIMAX, between 0.01 and 0.001 d, will force many iterations per day and thus smooth integration. However, a run with such a small DTIMAX takes a long time. See also DPSIMAX. [see WAT-ITER]


    ## Heat flow (unimplemented)
    # p_HEAT    = input_param[1,"HEAT"]
    # p_TopInfT = soil_discr["TopInfT"]
    # unused p_HeatCapOld = soil_discr["HeatCapOld"]

    # Isotope transport parameters
    p_VXYLEM = input_param[1,"VXYLEM_mm"] # mm, storage volume of well mixed xylem storage per ground area # TODO(bernhard): possibly link this to SAI...

    ## Location / Meteo
    p_NPINT  = 1 # Hardcoded. If p_NPINT>1, then multiple precipitation intervals would need
                 #            to be defined in an additional input data set PRECDAT.
    if p_NPINT == 1
        p_DTP = 1 / p_NPINT
    else
        error("Case with multiple precipitation intervals (using PRECDAT and precip_interval != 1) is not implemented.")
    end
    p_DURATN = input_storm_durations[1:12,"storm_durations_h"]# average storm duration for each month, hr
    p_LAT    = input_param[1,"LAT_DEG"]   /57.296  # (Location parameter), latitude, radians
    p_ASPECT = input_param[1,"ASPECT_DEG"]/57.296  # (Location parameter), aspect, radians through east from north
    p_ESLOPE = input_param[1,"ESLOPE_DEG"]/57.296  # (Location parameter), slope for evapotranspiration and snowmelt, radians
    p_MELFAC = input_param[1,"MELFAC"]    # (Location parameter), degree day melt factor for open land, MJ m-2 d-1 K-1
    p_RSTEMP = input_param[1,"RSTEMP"]    # (Location parameter), base temperature for snow-rain transition, °C
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
    p_FETCH  = input_param[1,"FETCH"]  # Fetch upwind of the weather station at which wind speed was measured, m (BROOK90: FETCH is fixed at 5000 m)
    p_WNDRAT = input_param[1,"WNDRAT"] # Average ratio of nighttime to daytime wind speed, dimensionless (BROOK90: WNDRAT is fixed at 0.3)
    p_Z0G    = input_param[1,"Z0G"]    # Ground surface roughness, m ()
    p_Z0S    = input_param[1,"Z0S"]    # Roughness parameter of snow surface, m (BROOK90: Z0S is fixed at 0.001 m)
    p_Z0W    = input_param[1,"Z0W"]    # Roughness parameter at the weather station where wind speed was measured, m
    p_ZMINH  = input_param[1,"ZMINH"]  # Reference height for weather data above canopy HEIGHT, m (BROOK90: ZMINH is fixed at 2 m)
    p_ZW     = input_param[1,"ZW"]     # Weather station measurement height for wind, m (BROOK90: W is fixed at 10 m)
    p_HS     = input_param[1,"HS"]     # Lower height limit, for roughness parameter interpolation, if canopy HEIGHT is below -> CZS (BROOK90: HS is fixed at 1 m)
    p_HR     = input_param[1,"HR"]     # Upper height limit, for roughness parameter interpolation, if canopy HEIGHT is above -> CZR (BROOK90: HR is fixed at 10 m)
    p_CZS    = input_param[1,"CZS"]    # Ratio of roughness parameter to HEIGHT (canopy) for HEIGHT < HS, dimensionless (BROOK90: CZS is fixed at 0.13)
    p_CZR    = input_param[1,"CZR"]    # Ratio of roughness parameter to HEIGHT (canopy) for HEIGHT > HR, dimensionless (BROOK90: CZR is fixed at 0.05)
    p_C1     = input_param[1,"C1"]     # intercept of linear relation between the ratio of actual to potential solar radiation for the day and sunshine duration (BROOK90: C1 is fixed at 0.25 following Brutsaert (1982))
    p_C2     = input_param[1,"C2"]     # slope     of linear relation between the ratio of actual to potential solar radiation for the day and sunshine duration (BROOK90: C2 is fixed at 0.5 following Brutsaert (1982))
    p_C3     = input_param[1,"C3"]     # Ratio of net longwave radiation for overcast sky (sunshine duration = 0) to that for clear sky (sunshine duration = 1). (BROOK90: C3 is fixed at 0.2 following Brutsaert (1982))
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
    p_SNODEN = input_param[1,"SNODEN"] # snowpack density, mm/mm (BROOK90: SNODEN is fixed at 0.3)
    p_CCFAC  = input_param[1,"CCFAC"]  # cold content factor, MJ m-2 d-1 K-1 (BROOK90: CCFAC is fixed at 0.3 MJ m-2 d-1 K-1)
    p_GRDMLT = input_param[1,"GRDMLT"] # rate of groundmelt of snowpack, mm/d (BROOK90: GRDMLT is fixed at 0.35 mm/d from Hubbard Brook (Federer 1965))
    p_LAIMLT = input_param[1,"LAIMLT"] # dependence of snowmelt on LAI, dimensionless (BROOK90: LAIMLT is fixed at 0.2)
    p_SAIMLT = input_param[1,"SAIMLT"] # dependence of snowmelt on SAI, dimensionless (BROOK90: SAIMLT is fixed at 0.5)
    p_MAXLQF = input_param[1,"MAXLQF"] # maximum liquid water fraction of SNOW, dimensionless (BROOK90: MAXLQF is fixed at 0.05)
    # Documentation from ecoshift:
    # SNODEN (Fixed parameter) - snow density, mm/mm. SNODEN is the snowpack density or the ratio of water content to depth. SNODEN is used only to correct leaf area index and stem area index to the fraction of them that is above the snow. It does not vary with time or snowpack ripeness in BROOK90. SNODEN is fixed at 0.3. [see PET-CANOPY]
    # CCFAC (Fixed parameter) - cold content factor, MJ m-2 d-1 K-1. CCFAC is a degree day factor for accumulation of cold content for a day with daylength of 0.5 d. It controls the snow energy balance when TA is less than 0°C. Larger values make the snow temperature lag farther behind the air temperature. Sensitivity of snowmelt to CCFAC is small unless CCFAC is less than 0.05 MJ m-2 d-1 K-1. CCFAC = 0 means there is no cold content and snow temperature is always 0°C. CCFAC is fixed at 0.3 MJ m-2 d-1 K-1. [see SNO-SNOENRGY]
    # GRDMLT (Fixed parameter) - rate of groundmelt of snowpack, mm/d. GRDMLT is the constant rate of melt at the bottom of the snowpack because of soil heat transfer to the snow. This parameter controls low flow during winter periods with snowpack. BROOK90 assumes that there is never any soil frost. GRDMLT is fixed at 0.35 mm/d from Hubbard Brook (Federer 1965). [see SNO-SNOWPACK]
    # LAIMLT (Fixed parameter) - dependence of snowmelt on LAI, dimensionless. Melt is linearly proportional to exp(-LAIMLT * LAI) so melt decreases exponentially as LAI increases. LAIMLT is fixed at 0.2. This value makes the exp factor 1.0, 0.67, 0.45, and 0.30 for LAI = 0, 2, 4, and 6 respectively. [see SNO-SNOENRGY]
    # SAIMLT (Fixed parameter) - snowmelt dependence on SAI, dimensionless. Melt is linearly proportional to exp (-SAIMLT * SAI); so melt decreases exponentially as SAI increases. SAIMLT is fixed at 0.5. With SAIMLT = 0.5 and LAIMLT = 0.2 , melt in deciduous forest (LAI=0, SAI=0.7) is 0.70 times that in the open, and melt in conifer forest (LAI=6, SAI=0.7) is 0.21 times that in the open. Federer et al. (1973b) generalized literature ratios to 0.5 for deciduous/open and 0.25 for conifer/open. Matching their results with the exponential algorithm would require greater difference between LAIMLT and SAIMLT, but it is difficult to see why they should differ much. [see SNO-SNOENRGY]
    # MAXLQF (Fixed parameter) - maximum liquid water fraction of SNOW, dimensionless. MAXLQF is the liquid water fraction of the snow water, SNOW, at which water drains. MAXLQF is fixed at 0.05, which is appropriate in most snow environments. [see SNO-SNOWPACK]


    ## Vegetation (canopy)
    ### Vegetation dimensions
    p_CR     = input_param[1,"CR"]     # (Canopy parameter), ratio of projected stem area index (SAI) to HEIGHT when DENSEF = 1, (SAI = CS * HEIGHT * DENSEF)
    p_RHOTP  = input_param[1,"RHOTP"]  # (Fixed  parameter), Ratio of total leaf area to projected area, dimensionless (BROOK90: RHOTP is fixed at 2)
    p_LWIDTH = input_param[1,"LWIDTH"] # (Canopy parameter), average leaf width, m
    p_LPC    = input_param[1,"LPC"]    # (Fixed  parameter), Minimum LAI defining a closed canopy, dimensionless (BROOK90: LPC is fixed at 4.0)
    ### Vegetation influence on atmosphere and meteorology
    p_KSNVP  = input_param[1,"KSNVP"]  # (Canopy parameter), reduction factor to reduce snow evaporation (SNVP), (0.05 - 1)
    p_ALBSN  = input_param[1,"ALBSN"]  # (Canopy parameter), albedo or surface reflectivity with snow on the ground, (typically 0.1-0.9)
    p_ALB    = input_param[1,"ALB"]    # (Canopy parameter), albedo or surface reflectivity without snow on the ground, (typically 0.1-0.3)
    p_NN     = input_param[1,"NN"]     # (Fixed  parameter), Wind/diffusivity extinction coefficient, dimensionless (BROOK90: NN is fixed at 2.5 following Shuttleworth and Gurney (1990))
    p_CS     = input_param[1,"CS"]     # (Canopy parameter), extinction coefficient for photosynthetically-active radiation in the canopy, (unit?)
    p_RM     = input_param[1,"RM"]     # (Fixed  parameter), Nominal maximum solar radiation possible on a leaf, W/m2 (BROOK90: RM is fixed at 1000 W/m2)
    ### Interception
    p_CINTRL = input_param[1,"CINTRL"] # (Fixed  parameter), Maximum interception storage of rain       per unit LAI, mm (BROOK90: CINTRL and CINTRS are both fixed at 0.15 mm)
    p_CINTRS = input_param[1,"CINTRS"] # (Fixed  parameter), Maximum interception storage of rain       per unit SAI, mm (BROOK90: CINTRL and CINTRS are both fixed at 0.15 mm)
    p_CINTSL = input_param[1,"CINTSL"] # (Fixed  parameter), Maximum interception storage of snow water per unit LAI, mm (BROOK90: CINTSL and CINTSS are both fixed at 0.6 mm)
    p_CINTSS = input_param[1,"CINTSS"] # (Fixed  parameter), Maximum interception storage of snow water per unit SAI, mm (BROOK90: CINTSL and CINTSS are both fixed at 0.6 mm)
    p_FRINTL = input_param[1,"FRINTLAI"] # (Fixed  parameter), Intercepted fraction of rain per unit LAI, dimensionless (BROOK90: FRINTLAI and FRINTSAI are both fixed at 0.06)
    p_FRINTS = input_param[1,"FRINTSAI"] # (Fixed  parameter), Intercepted fraction of rain per unit SAI, dimensionless (BROOK90: FRINTLAI and FRINTSAI are both fixed at 0.06)
    p_FSINTL = input_param[1,"FSINTLAI"] # (Fixed  parameter), Intercepted fraction of snow per unit LAI, dimensionless (BROOK90: FSINTLAI and FSINTSAI are both fixed at 0.04)
    p_FSINTS = input_param[1,"FSINTSAI"] # (Fixed  parameter), Intercepted fraction of snow per unit SAI, dimensionless (BROOK90: FSINTLAI and FSINTSAI are both fixed at 0.04)
    ### Vegetation conductivity
    p_MXKPL  = input_param[1,"MXKPL"]  # (Canopy parameter), maximum plant conductivity, mm d-1 MPa-1.
    p_FXYLEM = input_param[1,"FXYLEM"] # (Canopy parameter), fraction of plant resistance that is in the xylem (above ground) and not in the roots, (0-1)
    p_GLMAX  = input_param[1,"GLMAX"]  # (Canopy parameter), maximum leaf conductance, cm/s
    p_GLMIN  = input_param[1,"GLMIN"]  # (Canopy parameter), minimum leaf conductance, cm/s
    ### Stomatal regulation
    p_TH     = input_param[1,"TH"]     # (Canopy parameter), temperature controlling closing of stomates,°C
    p_T1     = input_param[1,"T1"]     # (Canopy parameter), temperature controlling closing of stomates,°C
    p_T2     = input_param[1,"T2"]     # (Canopy parameter), temperature controlling closing of stomates,°C
    p_TL     = input_param[1,"TL"]     # (Canopy parameter), temperature controlling closing of stomates,°C
    p_CVPD   = input_param[1,"CVPD"]   # (Fixed  parameter), Vapor pressure deficit at which stomatal conductance is halved, kPa (BROOK90: CVPD is fixed at 2 kPa for all cover types)
    p_R5     = input_param[1,"R5"]     # (Fixed  parameter), Solar radiation at which stomatal conductance is half of its value at RM, W/m2 (BROOK90: BROOK90 fixes R5 = 100 W/m2 as the default for all cover types)
    p_PSICR  = input_param[1,"PSICR"]  # (Canopy parameter), minimum plant leaf water potential, MPa.
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
    NOOUTF   = 1 == input_param[1,"NOOUTF"] # flag to prevent outflow from roots (hydraulic redistribution), (0/1)
    p_RTRAD  = input_param[1,"RTRAD"] # average root radius, mm (BROOK90: RTRAD is fixed at 0.35 mm)
    p_MXRTLN = input_param[1,"MXRTLN"] # (Canopy parameter), maximum length of fine roots per unit ground area, m/m2.
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
    FLAG_MualVanGen = simOption_FLAG_MualVanGen # 0 for Clapp-Hornberger; 1 for Mualem-van Genuchten
    p_RSSA   = input_param[1,"RSSA"] # Soil evaporation resistance (RSS) at field capacity, s/m (BROOK90: RSSA is fixed at 500 s/m following Shuttleworth and Gurney (1990))
    p_RSSB   = input_param[1,"RSSB"] # Exponent in relation of soil evaporation resistance (RSS) to soil water potential (PSIM) in the top layer, dimensionless, (BROOK90: RSSB is fixed at 1.0, which makes RSS directly proportional to PSIM)
    if FLAG_MualVanGen == 0
        p_soil = LWFBrook90.KPT.KPT_SOILPAR_Ch1d(;
            p_THICK = p_THICK,
            p_STONEF = soil_discr["STONEF"],           # stone volume fraction in each soil layer, dimensionless
            p_THSAT  = soil_discr["PAR"][!,"θs"],      # THETA at saturation, m3/m3
            p_THETAF = soil_discr["PAR"][!,"θf"],      # volumetric water content at "field capacity" corresponding to KF and PSIF for soil layer, m3/m3
            p_KF     = soil_discr["PAR"][!,"kf"],      # hydraulic conductivity at field capacity corresponding to THETAF and PSIF for a soil layer, mm/d
            p_PSIF   = soil_discr["PAR"][!,"ψf"],      # matric potential at "field capacity" corresponding to KF and THETAF for a soil layer, kPa
            p_BEXP   = soil_discr["PAR"][!,"bexp"],    # Clapp-Hornberger exponent for ψ-θ relation
            p_WETINF = soil_discr["PAR"][!,"wetinf"])  # wetness at dry end of near-saturation range for a soil layer, dimensionless

    elseif FLAG_MualVanGen == 1
        # Instantiate soil parameters
        p_soil = LWFBrook90.KPT.KPT_SOILPAR_Mvg1d(;
            p_THICK  = p_THICK,
            p_STONEF = soil_discr["STONEF"],
            p_THSAT  = soil_discr["PAR"][!,"θs"],
            p_Kθfc   = soil_discr["PAR"][!,"K(θ_fc)"],
            p_KSAT   = soil_discr["PAR"][!,"Ksat"],
            p_MvGα   = soil_discr["PAR"][!,"α"],
            p_MvGn   = soil_discr["PAR"][!,"n"],
            p_MvGl   = soil_discr["PAR"][!,"tort"],
            p_θr     = soil_discr["PAR"][!,"θr"])

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
    p_IMPERV = input_param[1,"IMPERV"] # (Flow parameter), fraction of impervious surface area generating surface or source area flow (SRFL), dimensionless
    p_INFEXP = input_param[1,"INFEXP"] # (Flow parameter), infiltration exponent that determines the distribution of infiltrated water with depth, dimensionless (from 0 to >1; 0 = all infiltration to top soil layer, 1 = uniform distribution down to ILAYER, >1 = more water in lower layers closer to ILAYER)
    p_ILAYER = soil_discr["ILAYER"] # (Flow parameter), number of layers over which infiltration is distributed
    p_QLAYER = soil_discr["QLAYER"] # (Flow parameter), number of soil layers for SRFL
    p_INFRAC = LWFBrook90.WAT.INFPAR(p_INFEXP, p_ILAYER, p_soil, NLAYER) # fraction of (preferential) infiltration to each layer
    # TODO(bernhard):switch to ILAYAER and QLAYER to IDEPTH_m and QDEPTH_m, which are independent of soil discretization.

    ### Flow generation
    p_BYPAR  = input_param[1,"BYPAR"]  # (Flow parameter), flag to activate bypass flow (BYFL), (0/1)
    p_DRAIN  = input_param[1,"DRAIN"] # (Flow parameter), continuous flag to activate drainge VRFLI(n), (between 0 and 1; 1 = gravity drainage, 0 = no drainage)
    p_DSLOPE = input_param[1,"DSLOPE"] # (Flow parameter), hillslope angle for downslope matric flow (DSFL), degrees
    p_GSC    = input_param[1,"GSC"] # (Flow parameter), fraction of groundwater storage (GWAT), that is transferred to groundwater flow (GWFL) and deep seepage (SEEP) each day, d-1
    p_GSP    = input_param[1,"GSP"] # (Flow parameter), fraction of groundwater discharge produced by GSC that goes to deep seepage (SEEP) and is not added to streamflow (FLOW), dimensionless
    p_LENGTH_SLOPE = input_param[1,"LENGTH_SLOPE"] # (Flow parameter), slope length for downslope flow (DSFL), m
    p_QFFC   = input_param[1,"QFFC"] # (Flow parameter), quick flow fraction for SRFL and BYFL at THETAF, dimensionless
    p_QFPAR  = input_param[1,"QFPAR"] # (Flow parameter), raction of the water content between field capacity (THETAF) and saturation (THSAT) at which the quick flow fraction is 1, dimensionless

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
    Reset = false # hardcoded
    p_cst_2 = (NLAYER, FLAG_MualVanGen, compute_intermediate_quantities, Reset,
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
        p_VXYLEM)


    p_cst_4 = (
        FLAG_MualVanGen,
        compute_intermediate_quantities,
        simulate_isotopes,
        row_idx_scalars = [],
        row_idx_SWATI   = [],
        row_idx_RWU     = [],
        row_idx_accum   = [],
        names_accum     = [],
        col_idx_d18O    = [],
        col_idx_d2H     = []) # Placeholders will be overwritten when defining u0
    p_cst = (p_cst_1, p_cst_2, p_cst_3, p_cst_4)

    # 2b) Time varying parameters (e.g. meteorological forcings)
    p_fT = (p_DOY          = (t) -> LWFBrook90.p_DOY(t,    interpolated_meteoveg.REFERENCE_DATE),
            p_MONTHN       = (t) -> LWFBrook90.p_MONTHN(t, interpolated_meteoveg.REFERENCE_DATE),
            p_GLOBRAD      = interpolated_meteoveg.p_GLOBRAD,
            p_TMAX         = interpolated_meteoveg.p_TMAX,
            p_TMIN         = interpolated_meteoveg.p_TMIN,
            p_VAPPRES      = interpolated_meteoveg.p_VAPPRES,
            p_WIND         = interpolated_meteoveg.p_WIND,
            p_PREC         = interpolated_meteoveg.p_PREC,
            p_DENSEF       = interpolated_meteoveg.p_DENSEF, # canopy density multiplier between 0.05 and 1, dimensionless
            p_HEIGHT       = interpolated_meteoveg.p_HEIGHT,
            p_LAI          = interpolated_meteoveg.p_LAI,
            p_SAI          = interpolated_meteoveg.p_SAI,
            p_AGE          = interpolated_meteoveg.p_AGE,
            p_RELDEN       = interpolated_meteoveg.p_RELDEN,
            p_d18OPREC     = interpolated_meteoveg.p_d18OPREC,
            p_d2HPREC      = interpolated_meteoveg.p_d2HPREC,
            REFERENCE_DATE = interpolated_meteoveg.REFERENCE_DATE)
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
    return ((soil_discr["PSIM_init"], soil_discr["d18O_init"],soil_discr["d2H_init"]),
            (p_cst, p_fT, p_fu, p_cache))
end






"""
    interpolate_meteoveg(
        input_meteoveg::DataFrame,
        input_meteoveg_reference_date::DateTime,
        input_meteoiso::DataFrame, )

Take climate and vegetation parameters in `input_meteoveg` and `input_meteoiso` and generates continuous parameters.
"""
function interpolate_meteoveg(
    input_meteoveg::DataFrame,
    input_meteoveg_reference_date::DateTime,
    input_meteoiso::Union{DataFrame,Nothing},
    # root density parameters
    NLAYER,
    p_INITRDEP,
    p_INITRLEN,
    p_RGROPER,
    p_tini,
    p_frelden,
    p_MAXLAI,
    p_SAI_baseline_,
    p_DENSEF_baseline_,
    p_AGE_baseline_yrs,
    p_HEIGHT_baseline_m)

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
    #     scatter(input_meteoiso.days, input_meteoiso.delta18O_permil)
    #     # # plot!(ts, scale(interpolate(input_meteoiso.delta18O_permil, (BSpline(Constant{Previous}()))), time_range)(ts), label = "PRECIN {Previous}",  xlims=(0,30))
    #     # # plot!(ts, scale(interpolate(input_meteoiso.delta18O_permil, (BSpline(Constant{Next}()))),     time_range)(ts), label = "PRECIN {Next}",      xlims=(0,30))
    #     # plot!(ts .+ 1,
    #     #     scale(interpolate(input_meteoiso.delta18O_permil, (BSpline(Constant{Next}()))),     time_range)(ts), label = "PRECIN+1 {Previous}",      xlims=(0,30))
    #     # # plot!(ts, scale(interpolate(input_meteoiso.delta18O_permil, (BSpline(Constant()))),           time_range)(ts), label = "PRECIN {Nearest}",   xlims=(0,30))
    #     # # plot!(ts, scale(interpolate(input_meteoiso.delta18O_permil, (BSpline(Linear()))),             time_range)(ts), label = "PRECIN Linear",   xlims=(0,30))
    #     # plot!(ts, scale(interpolate(input_meteoiso.delta18O_permil, (BSpline(Linear()))),             time_range)(ts), label = "PRECIN Linear",   xlims=(0,30))
    #     plot!(ts, extrapolate(interpolate((input_meteoiso.days, ), input_meteoiso.delta18O_permil, Gridded(Linear())), Throw())(ts), label = "Gridded(Linear()")
    #     plot!(ts, extrapolate(interpolate((input_meteoiso.days, ), input_meteoiso.delta18O_permil, Gridded(Constant())), Throw())(ts), label = "Gridded(Constant()")
    #     plot!(ts, extrapolate(interpolate((input_meteoiso.days, ), input_meteoiso.delta18O_permil, Gridded(Constant{Next}())), Throw())(ts), label = "Gridded(Constant{Previous}()")
    #     xlims!(0, 500)
    # ### END FOR DEVELOPMENT

    # # Interpolate regular data with Previous.
    # # (e.g. PREC reported for 2010-01-01 is rain that has fallen on that day and rate is applied until 2010-01-02.)
    # # Note that we need to shift it slightly ahead so that p_GLOBRAD(1.0) is equal to p_GLOBRAD(1.1)
    # # Interpolate input data without modification
    # p_GLOBRAD = extrapolate(interpolate((input_meteoveg.days  .- 0.00001, ), input_meteoveg.GLOBRAD,          Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_TMAX    = extrapolate(interpolate((input_meteoveg.days  .- 0.00001, ), input_meteoveg.TMAX,             Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_TMIN    = extrapolate(interpolate((input_meteoveg.days  .- 0.00001, ), input_meteoveg.TMIN,             Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_VAPPRES = extrapolate(interpolate((input_meteoveg.days  .- 0.00001, ), input_meteoveg.VAPPRES,          Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_WIND    = extrapolate(interpolate((input_meteoveg.days  .- 0.00001, ), input_meteoveg.WIND,             Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_PREC    = extrapolate(interpolate((input_meteoveg.days  .- 0.00001, ), input_meteoveg.PRECIN,           Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # # NOTE: PRECIN is already in mm/day from the input data set. No transformation is needed for p_PREC.

    # # Interpolate input data with consideration of baselines
    # p_LAI     = extrapolate(interpolate((input_meteoveg.days  .- 0.00001, ), input_meteoveg.LAI./100 .* p_MAXLAI,              Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_SAI     = extrapolate(interpolate((input_meteoveg.days  .- 0.00001, ), input_meteoveg.SAI./100 .* p_SAI_baseline_,       Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_DENSEF  = extrapolate(interpolate((input_meteoveg.days  .- 0.00001, ), input_meteoveg.DENSEF./100 .* p_DENSEF_baseline_, Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())
    # p_HEIGHT  = extrapolate(interpolate((input_meteoveg.days  .- 0.00001, ), input_meteoveg.HEIGHT./100 .* p_HEIGHT_baseline_m,Gridded(Constant{Previous}())), Throw()) #extrapolate flat, alternative: Throw())

    ### Quickfix:
    @assert 1==length(unique(diff(input_meteoveg.days))) """
        Error: LWFBrook90.jl expects the input data in meteoveg.csv to be regularly spaced in time.
    """
    @assert unique(diff(input_meteoveg.days)) == [1.0] """
        Error: LWFBrook90.jl expects the input data in meteoveg.csv to provided as daily values.
    """
    time_range = range(minimum(input_meteoveg.days), maximum(input_meteoveg.days), length=length(input_meteoveg.days))
    p_GLOBRAD = extrapolate(scale(interpolate(input_meteoveg.GLOBRAD, (BSpline(Constant{Previous}()))), time_range) ,0)
    p_TMAX    = extrapolate(scale(interpolate(input_meteoveg.TMAX,    (BSpline(Constant{Previous}()))), time_range) ,0)
    p_TMIN    = extrapolate(scale(interpolate(input_meteoveg.TMIN,    (BSpline(Constant{Previous}()))), time_range) ,0)
    p_VAPPRES = extrapolate(scale(interpolate(input_meteoveg.VAPPRES, (BSpline(Constant{Previous}()))), time_range) ,0)
    p_WIND    = extrapolate(scale(interpolate(input_meteoveg.WIND,    (BSpline(Constant{Previous}()))), time_range) ,0)
    p_PREC    = extrapolate(scale(interpolate(input_meteoveg.PRECIN,  (BSpline(Constant{Previous}()))), time_range) ,0)

    p_LAI     = extrapolate(scale(interpolate(input_meteoveg.LAI./100 .* p_MAXLAI,              (BSpline(Constant{Previous}()))), time_range) ,0)
    p_SAI     = extrapolate(scale(interpolate(input_meteoveg.SAI./100 .* p_SAI_baseline_,       (BSpline(Constant{Previous}()))), time_range) ,0)
    p_DENSEF  = extrapolate(scale(interpolate(input_meteoveg.DENSEF./100 .* p_DENSEF_baseline_, (BSpline(Constant{Previous}()))), time_range) ,0)
    p_HEIGHT  = extrapolate(scale(interpolate(input_meteoveg.HEIGHT./100 .* p_HEIGHT_baseline_m,(BSpline(Constant{Previous}()))), time_range) ,0)
    ###
    # p_AGE     = extrapolate(interpolate((input_meteoveg.days  .- 0.00001, ), input_meteoveg.AGE, Gridded(Constant{Next}())), Flat()) #extrapolate flat, alternative: Throw())
    p_AGE     = (t) -> p_AGE_baseline_yrs + t/365

    # Note that meteoiso does not need to be regularly spaced:
    # Interpolate irregular data with Next and shift to noon

    if (isnothing(input_meteoiso))
        function p_d18OPREC(t)
            return nothing
        end
        function p_d2HPREC(t)
            return nothing
        end
    else
        ## Quickfix:
        # @assert 1==length(unique(diff(input_meteoiso.days))) """
        #     Error: LWFBrook90.jl expects the input data in meeoiso.csv to be regularly spaced in time.
        # """
        # @assert unique(diff(input_meteoiso.days)) == [1.0] """
        #     Error: LWFBrook90.jl expects the input data in meeoiso.csv to provided as daily values.
        # """
        # time_range_iso = range(minimum(input_meteoiso.days), maximum(input_meteoiso.days), length=length(input_meteoiso.days))
        # p_d18OPREC= extrapolate(scale(interpolate(input_meteoiso.delta18O_permil,    (BSpline(Constant{Previous}()))), time_range_iso) ,0)
        # p_d2HPREC = extrapolate(scale(interpolate(input_meteoiso.delta2H_permil,     (BSpline(Constant{Previous}()))), time_range_iso) ,0)
        ## End Quickfix:
        # Here we shift the days of collection of the cumulative samples to noon. => + 12/24 days
        p_d18OPREC= extrapolate(interpolate((input_meteoiso.days .+ 12/24, ), input_meteoiso.delta18O_permil,    (Gridded(Constant{Next}()))), Flat()) #extrapolate flat, alternative: Throw())
        p_d2HPREC = extrapolate(interpolate((input_meteoiso.days .+ 12/24, ), input_meteoiso.delta2H_permil,     (Gridded(Constant{Next}()))), Flat()) #extrapolate flat, alternative: Throw())
    end

    # 2a Compute time dependent root density parameters
    # Which is a vector quantity that is dependent on time:
    p_RELDEN_2Darray = fill(NaN, nrow(input_meteoveg), NLAYER)
    for i in 1:nrow(input_meteoveg)
        p_RELDEN_2Darray[i,:] =
            LWFBrook90.WAT.LWFRootGrowth(
                p_frelden,
                p_tini,
                p_AGE(input_meteoveg.days[i]),
                p_RGROPER,
                p_INITRDEP,
                p_INITRLEN,
                NLAYER)
    end
    p_RELDEN =  extrapolate(interpolate((input_meteoveg.days, 1:NLAYER), p_RELDEN_2Darray,
                                        (Gridded(Constant{Next}()), NoInterp()), # 1st dimension: ..., 2nd dimension NoInterp()
                                        ), Flat()) # extrapolate flat, alternative: Throw()
    return (p_GLOBRAD      = p_GLOBRAD,
            p_TMAX         = p_TMAX,
            p_TMIN         = p_TMIN,
            p_VAPPRES      = p_VAPPRES,
            p_WIND         = p_WIND,
            p_PREC         = p_PREC,
            p_d18OPREC     = p_d18OPREC,
            p_d2HPREC      = p_d2HPREC,
            p_DENSEF       = p_DENSEF,
            p_HEIGHT       = p_HEIGHT,
            p_LAI          = p_LAI,
            p_SAI          = p_SAI,
            p_AGE          = p_AGE,
            p_RELDEN       = p_RELDEN,
            REFERENCE_DATE = input_meteoveg_reference_date)
end

"""
    discretize_soil_params(
        input_soil_horizons,
        input_soil_discretization,
        soil_output_depths,
        IDEPTH_m,
        QDEPTH_m,
        INITRDEP,
        RGRORATE,
        FLAG_MualVanGen)

Discretize soil domain into computational layers and attribute soil parameters based on the defined horizons.
Densify discretization whenever an interface or an additional layer is needed. This is the case for interfaces
when a new soil horizons begins or an additional computational layer (consisting of upper and lower interfaces)
when a state variable needs to be extracted at a specified output depth.
"""
function discretize_soil_params(
    input_soil_horizons,
    input_soil_discretization,
    soil_output_depths,
    IDEPTH_m,
    QDEPTH_m,
    INITRDEP,
    RGRORATE,
    FLAG_MualVanGen)

    # This function maps the input parameters: (input_soil_horizons, ILAYER, QLAYER, INITRDEP, RGRORATE)
    # onto the soil discretization (input_soil_discretization)

    # input_soil_horizons: A matrix of the 8 soil materials parameters.
    # input_soil_discretization:
    # IDEPTH_m depth over which infiltration is distributed, [m] "IDEPTH_m determines the
    #        number of soil layers over which infiltration is distributed when INFEXP is
    #        greater than 0."
    # QDEPTH_m soil depth for SRFL calculation, 0 to prevent SRFL, [m] "QDEPTH_m determines the
    #        layers over which wetness is calculated to determine source area (SRFL)
    #        parameters."
    # INITRDEP
    # RGRORATE

    ############
    # 1) Check if discretization needs to be refined
    # 1a) Find out which layers need to be added
    layers_to_insert = soil_output_depths

    # ε = 0.001 # thickness of layer to be inserted, [m]
    # ε = 0.025 # thickness of layer to be inserted, [m]
    ε = 0.050 # thickness of layer to be inserted, [m]

    needed_interfaces_for_additional_layers = zeros(Float64, 0)
    for layer in layers_to_insert
        println(layer)
        if any(layer .∈ input_soil_discretization.Lower_m)
            # specific depth is already an interface in the discretization
            # -> only add a lower interface ε below the upper
            append!(needed_interfaces_for_additional_layers, round(layer - ε; digits=3))
        else
            # specific depth is not yet an interface
            # -> add upper and lower interface
            append!(needed_interfaces_for_additional_layers, round(layer;     digits=3))
            append!(needed_interfaces_for_additional_layers, round(layer - ε; digits=3))
            # TODO: would technically need to check that no interface between layer and layer + ε
        end
    end

    # 1b) Find out which interfaces need to be added
    all_needed_interfaces = unique(sort(
        [
         # interfaces for change in horizons or inflow depths or QDEPTH_m
         input_soil_horizons[!,"Upper_m"];
         input_soil_horizons[!,"Lower_m"];
         -IDEPTH_m; # m (positive) to m (negative)
         -QDEPTH_m;  # m (positive) to m (negative)

         # interfaces for layers determined in 1a)
         needed_interfaces_for_additional_layers
         ];
        rev = true))

    existing_interfaces = [0; input_soil_discretization[!,"Lower_m"]]
    to_add = all_needed_interfaces[(!).(all_needed_interfaces .∈ (existing_interfaces,))]

    # Add them to the DataFrame
    soil_discretization = copy(input_soil_discretization) # otherwise input argument is modified in-place
    # soil_discretization = input_soil_discretization
    for interface_to_add in to_add
        for i in 1:nrow(soil_discretization)
            condition = interface_to_add > soil_discretization[i, "Lower_m"]
            #println("$i, interface_to_add:$interface_to_add, $(soil_discretization[i, "Lower_m"]) $condition")
            if (condition)
                interface_to_add
                pos = i+1
                insert!.(eachcol(soil_discretization),
                        pos,
                        Matrix(soil_discretization)[pos-1,:])
                soil_discretization[pos-1, "Lower_m"] = interface_to_add[1]
                soil_discretization[pos,   "Upper_m"] = interface_to_add[1]
                break
            elseif (i == nrow(soil_discretization))
                # If interface to add is below the lowest layer: add it only in case it is about ε lower
                if (interface_to_add >= soil_discretization[i, "Lower_m"] - ε)
                    pos = i+1
                    insert!.(eachcol(soil_discretization),
                        pos,
                        Matrix(soil_discretization)[pos-1,:])
                    soil_discretization[pos, "Upper_m"] = soil_discretization[pos-1, "Lower_m"]
                    soil_discretization[pos, "Lower_m"] = interface_to_add[1]
                    break
                end
            end
        end
    end

    # Define ILAYER and QLAYER (to be used internally instead of IDEPTH_m and QDEPTH_m)
    is_infiltration_layer_BOOLEAN = -IDEPTH_m .<= soil_discretization[!,"Lower_m"]
    is_SRFL_layer_BOOLEAN         = -QDEPTH_m .<= soil_discretization[!,"Lower_m"]
    ILAYER = sum(is_infiltration_layer_BOOLEAN) # lowest node where Bottom_m is below IDEPTH_m
    QLAYER = sum(is_SRFL_layer_BOOLEAN)         # lowest node where Bottom_m is below QDEPTH_m

    if (-IDEPTH_m < soil_discretization[end,"Lower_m"]) ||
        (-QDEPTH_m < soil_discretization[end,"Lower_m"])
        @error "QDEPTH_m or IDEPTH_m were defined deeper than the lowest simulation element."
    end
    ############

    ############
    # 2) Append soil horizon data to discretized soil domain
    # Assert expectations
    @assert input_soil_horizons[1,"Upper_m"] ≈ 0
    @assert soil_discretization[1,"Upper_m"] ≈ 0
    @assert input_soil_horizons[1,"Lower_m"] < 0
    @assert soil_discretization[1,"Lower_m"] < 0

    # Find for each soil_discretization node in which horizon it lies
   which_horizon = fill(0, nrow(soil_discretization))
    for i = 1:nrow(soil_discretization)
        # For each soil discretizations (i) check to which which horizon it belongs
        idx = findfirst(soil_discretization[i,"Lower_m"] .>= input_soil_horizons[:,"Lower_m"])

        if (isnothing(idx) && i == nrow(soil_discretization))
            idx = nrow(input_soil_horizons)
        end

        which_horizon[i] = input_soil_horizons[idx,"HorizonNr"]
    end
    soil_discretization[:,"HorizonNr"] = which_horizon

    # If lowest discretized layer is only ε lower than defined soil horizons extrapolate
    # the properties of the lowest soil horizon to that infinitesimal thin layer
    if ( (soil_discretization[end,"HorizonNr"] == 0) &
        (abs(soil_discretization[end,"Lower_m"] - input_soil_horizons[end,"Lower_m"]) < ε) )
        soil_discretization[end,"HorizonNr"] = nrow(input_soil_horizons)
    end

    nlayer_before_join = nrow(soil_discretization)
    soil_discretization = innerjoin(soil_discretization,
                                    input_soil_horizons, makeunique=true,
                                    #select(input_soil_horizons, Not([:Upper_m,:Lower_m])),
                                    on = :HorizonNr)
                                    # NOTE: not using leftjoin becaus it transforms types to Union{Missing, Float64} instead of only Float64
                                    # Therefore, use innerjoin and check manually that no rows in soil_discretization are lost
    @assert nrow(soil_discretization) == nlayer_before_join """
        When merging soil properties defined in input CSV-file containing soil horizons, some layers from
        the soil discretization were lost. This is most likely due to the soil horizons not covering the entire discretization domain.
        (Note that the discretized domain is defined in a separate input CSV-file.)

        Please check and correct the input files.
        """
        # (Note that the discretized domain is defined in a separate input CSV-file and is further extended to include internal variables IDEPTH_m and QDEPTH_m.)
    ############

    ############

    # hardcoded:
    HEAT = 0 # flag for heat balance; not implemented; input_param[1,"HEAT"],

    THICK_m   = soil_discretization[!,"Upper_m"] - soil_discretization[!,"Lower_m"] # thickness of soil layer [m]
    THICK     = 1000*(THICK_m)                                    # thickness of soil layer [mm]
    PSIM_init = soil_discretization[!,"uAux_PSIM_init_kPa"]       # initial condition PSIM [kPa]
    d18O_soil_init = soil_discretization[!,"u_delta18O_init_permil"] # initial condition soil water δ18O [‰]
    d2H_soil_init  = soil_discretization[!,"u_delta2H_init_permil"]  # initial condition soil water δ2H [‰]

    @assert all(PSIM_init .<= 0) "Initial matrix psi must be negative or zero"

    NLAYER = nrow(soil_discretization)

    for i = 1:NLAYER
        if (HEAT != 0)
            @error "HEAT must be set to zero, as heat balance is not implemented."
        end
        # if (HEAT == 1)
        #     # TemperatureNew(i)       = soil_discretization[i,7] we don't have it in the input file!!!
        #     if i > NLAYER
        #         MUE[i] = THICK[i] / ( THICK[i] + THICK(I+1) )
        #         ZL[i]  = 0.5 * ( THICK[i] + THICK(I+1) )
        #     else
        #         MUE[i] = 0.5
        #         ZL[i]  = THICK[i]
        #     end
        #     TMean[i]   = 0.
        # end
    end
    ############

    ############
    # 3) Bring to simple vecotr and matrix form

    # Parse the stonefraction and hydraulic parameters for each layer
    # depending on wheter we use FLAG_MualVanGen ==1 or ==2
    STONEF = soil_discretization[:,"gravel_volFrac"]

    if FLAG_MualVanGen == 1
        # Mualem-van Genuchten
        PAR = select(soil_discretization,
                :ths_volFrac    => :θs,
                :ksat_mmDay     => :Ksat,
                :alpha_perMeter => :α,
                :npar_          => :n,
                :tort_          => :tort,
                :thr_volFrac    => :θr
                )
        # ..it's a sin......
        # Hard default of 2 [mm d-1] for saturated hydraulic conductivity at field capacity K(θ_fc)
        PAR[:,"K(θ_fc)"] .= 2. #insertcols!(copy(soil_discretization), ("K(θ_fc)" => 2.)),
    end
    if FLAG_MualVanGen == 0
        # Clapp-Hornberger
        PAR = select(soil_discretization,
                :thsat    => :θs,
                :thetaf   => :θf,
                :kf_mmD   => :kf,
                :psif_kPa => :ψf,
                :bexp     => :bexp,
                :wetinf   => :wetinf
                )
    end

    # find thickness of maximum root zone
    # frelden: relative values of final root density per unit volume
    frelden = soil_discretization[!,"Rootden_"]            # root density
    dep     = soil_discretization[!,"Upper_m"] - THICK/1000/2 # soil depth [m] (midpoint, i.e. average depth of layer)
    depmax  = dep[1] - THICK[1] / 1000.

    i1 = findfirst(frelden .> 1.e-6) # first layer where frelden[i] is >=1.e-6
    i2 = findlast(frelden .> 1.e-6)  # last layer (in 1:NLAYER) where frelden[i] is >=1.e-6
    if !isnothing(i1) && isnothing(i2)
        i2 = NLAYER
    end

    tini = fill(NaN, NLAYER)
    for i = 1:NLAYER
        tini[i] = 1.e+20 # initial time for root growth in layer
        if i >= i1 && i <= i2
            frelden[i] = max( frelden[i], 1.01e-6)
        end
        if frelden[i] >= 1.e-6 && (depmax-dep[i]) <= INITRDEP
            tini[i] = 0.
        end
        if frelden[i] >= 1.e-6 && (depmax-dep[i]) > INITRDEP
            if RGRORATE > 0
                tini[i] = (depmax-dep[i]-INITRDEP)/RGRORATE
            end
        end
        # write(*,*)'dep= ',dep[i],' tini= ',tini[i]
    end


    # heat flow -------
    # nmat   = nrow(input_soil_horizons)
    if (HEAT != 0) @error "HEAT must be set to zero, as heat balance is not implemented." end
    #       if (HEAT == 1)
    #        READ (12,*) Comment
    #        READ (12,*) tTop, Comment
    #        tTop=tTop
    #        READ (12,*) tBot, Comment
    #        tBot=tBot
    #        READ (12,*) TopInfT, Comment
    #        READ (12,*) BotInfT, Comment
    #        READ (12,*) kTopT, Comment
    #        READ (12,*) kBotT, Comment
    #        DO 207 I = 1, 7
    #         READ (12,*) Comment
    # 207    CONTINUE
    #        DO 208 I = 1, nmat
    #         READ (12,*) ilay, SV[i], OV[i], HB1[i], HB2[i], HB3[i]
    #         TPar[1,I) = SV[i]
    #         TPar[2,I) = OV[i]
    #         TPar[3,I) = THDis
    # C        thermal conductivities -- transfer from [J m-1 s-1 K-1] to  [J mm-1 d-1 K-1]
    #         TPar[4,I) = HB1[i] * 86.400
    #         TPar[5,I) = HB2[i] * 86.400
    #         TPar[6,I) = HB3[i] * 86.400
    # C         volumetric heat capacities for solid, water and organic -- transfer from [MJ m-2 mm-1 K-1] to [J mm-3 K-1]
    #         TPar[7,I) = p_CVSOL   # To define in module_CONSTANTS.jl: p_CVSOL = ?  # CVSOL  - volumetric heat capacity of solid (MJ m-2 mm-1 K-1)
    #         TPar[8,I) = p_CVORG   # To define in module_CONSTANTS.jl: p_CVORG = ?  # volumetric heat capacity of organic material (MJ m-2 mm-1 K-1) (hillel98)
    #         TPar[9,I) = LWFBrook90.CONSTANTS.p_CVLQ
    # 208    CONTINUE
    #        READ (12,*) C
    #       end

    # mat       = input_soil_horizons[!,"mat"]
    # ### # from LWFBrook90R:md_brook90.f95 (lines 105)
    # TPar = fill(NaN, (nmat, 10))
    # TPar .= -1
    # HeatCapOld = fill(NaN, NLAYER)
    # for i = 1:NLAYER
    #     HeatCapOld[i] = TPar[mat[i],7] * TPar[mat[i],1] + TPar[mat[i],8]
    #                 # TPar[2,mat[i])+TPar[9,mat[i])*SWATI[i]/THICK[i]
    # end
    ###

    return Dict([
                ("NLAYER",NLAYER),
                ("ILAYER",ILAYER),
                ("QLAYER",QLAYER),
                ("THICK",THICK),
                ("PSIM_init",PSIM_init),
                ("d18O_init",d18O_soil_init),
                ("d2H_init",d2H_soil_init),
                ("frelden",frelden),
                ("PAR",PAR),
                ("STONEF",STONEF),
                ("tini",tini)])#,
                #
                #("HeatCapOld",HeatCapOld),
                #("TopInfT", TopInfT)])
end