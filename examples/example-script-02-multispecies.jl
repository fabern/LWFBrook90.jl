# Proposed input shape; this script constructs data only.
# The species, site_params, canopy_mixing and ground_snowpack keywords are proposed,
# and are not implemented API.
# See ../docs/multispecies-design.md for ownership, area units and solver changes.

multispecies_input = (
    canopy_mixing = :tiles,
    ground_snowpack = :shared, # separate canopy/interception; common ground snow
    site_params = (
        LAT_DEG = 47.0, ESLOPE_DEG = 0.0, ASPECT_DEG = 0.0,
        IDEPTH_m = 0.4, QDEPTH_m = 0.0, SLVPDEPTH_m = 0.0,
        DISPERSIVITY_mm = 40.0,
    ),
    Δz_thickness_m = fill(0.1, 11),
    storm_durations_h = fill(4.0, 12),
    IC_soil = (
        PSIM_init_kPa = -6.3,
        delta18O_init_permil = -13.0,
        delta2H_init_permil = -95.0,
    ),
    IC_scalar = (
        # Ground stores are per total stand area, with one shared snow history.
        amount = (
            u_GWAT_init_mm = 1.0, u_SNOW_init_mm = 0.0,
            u_CC_init_MJ_per_m2 = 0.0, u_SNOWLQ_init_mm = 0.0,
        ),
        d18O = (u_GWAT_init_permil = -13.0, u_SNOW_init_permil = -13.0),
        d2H = (u_GWAT_init_permil = -95.0, u_SNOW_init_permil = -95.0),
    ),
    species = (
        Beech = (
            cover_fraction = 0.60,
            # All extensive plant inputs are per square metre of Beech tile.
            params = (
                MAXLAI = 3.0, SAI_baseline_ = 1.0, DENSEF_baseline_ = 1.0,
                HEIGHT_baseline_m = 25.0, AGE_baseline_yrs = 100.0,
                LWIDTH = 0.08, RHOTP = 2.0, CR = 0.60,
                MXRTLN = 3000.0, RTRAD = 0.35,
                INITRLEN = 12.0, INITRDEP = 0.25,
                RGRORATE = 0.03, RGROPER = 30.0, NOOUTF = 1,
                MXKPL = 15.64, FXYLEM = 0.50, PSICR = -1.50,
                GLMAX = 0.00868, GLMIN = 0.0003,
                R5 = 287.0, CVPD = 2.0, RM = 1000.0,
                TL = 0.0, T1 = 10.0, T2 = 30.0, TH = 40.0,
                CINTRL = 0.15, CINTRS = 0.15, CINTSL = 0.60, CINTSS = 0.60,
                FRINTLAI = 0.06, FRINTSAI = 0.06,
                FSINTLAI = 0.04, FSINTSAI = 0.04,
                VXYLEM_mm = 20.0,
            ),
            root_distribution = (beta = 0.98, z_rootMax_m = -1.1),
            canopy_evolution = (
                DENSEF_rel = 100.0, HEIGHT_rel = 100.0, SAI_rel = 100.0,
                LAI_rel = (
                    DOY_Bstart = 120, Bduration = 21,
                    DOY_Cstart = 270, Cduration = 60,
                    LAI_perc_BtoC = 100.0, LAI_perc_CtoB = 0.0,
                ),
            ),
            initial_conditions = (
                INTR = (mm = 0.0, d18O = -13.0, d2H = -95.0),
                INTS = (mm = 0.0, d18O = -13.0, d2H = -95.0),
                # Xylem volume is prescribed by VXYLEM_mm, not a second amount IC.
                XYLEM = (d18O = -13.0, d2H = -95.0),
            ),
        ),
        Spruce = (
            cover_fraction = 0.40,
            # Spruce owns a complete configuration with independent defaults.
            params = (
                MAXLAI = 4.5, SAI_baseline_ = 1.5, DENSEF_baseline_ = 1.0,
                HEIGHT_baseline_m = 28.0, AGE_baseline_yrs = 80.0,
                LWIDTH = 0.04, RHOTP = 3.0, CR = 0.50,
                MXRTLN = 4500.0, RTRAD = 0.25,
                INITRLEN = 12.0, INITRDEP = 0.25,
                RGRORATE = 0.02, RGROPER = 30.0, NOOUTF = 1,
                MXKPL = 10.50, FXYLEM = 0.40, PSICR = -2.20,
                GLMAX = 0.00550, GLMIN = 0.0002,
                R5 = 200.0, CVPD = 1.5, RM = 1000.0,
                TL = 0.0, T1 = 8.0, T2 = 25.0, TH = 35.0,
                CINTRL = 0.18, CINTRS = 0.18, CINTSL = 0.60, CINTSS = 0.60,
                FRINTLAI = 0.06, FRINTSAI = 0.06,
                FSINTLAI = 0.04, FSINTSAI = 0.04,
                VXYLEM_mm = 25.0,
            ),
            root_distribution = (beta = 0.92, z_rootMax_m = -0.4),
            canopy_evolution = (
                DENSEF_rel = 100.0, HEIGHT_rel = 100.0, SAI_rel = 100.0,
                LAI_rel = (
                    DOY_Bstart = 120, Bduration = 21,
                    DOY_Cstart = 270, Cduration = 60,
                    LAI_perc_BtoC = 100.0, LAI_perc_CtoB = 100.0,
                ),
            ),
            initial_conditions = (
                INTR = (mm = 0.0, d18O = -13.0, d2H = -95.0),
                INTS = (mm = 0.0, d18O = -13.0, d2H = -95.0),
                XYLEM = (d18O = -13.0, d2H = -95.0),
            ),
        ),
    ),
)

# Once implemented, the intended loading call would be:
# model = loadSPAC("shared-input-folder", "site-prefix"; multispecies_input...)
# sim = setup(model)
# simulate!(sim)
