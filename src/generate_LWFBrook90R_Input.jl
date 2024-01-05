#################
# Task: Generate input files for LWFBrook90R from LWFBrook90.jl
# Description: This code defines a julia function that takes a DiscretizedSPAC
#              and generates the corrsponding input files for LWFBrook90R::run_LWFB90()
#################

# using RCall: # required for running LWFBrook90R directly see ?prepare_for_LWFBrook90R

# Note that LWFBrook90.jl does not cover all the use cases of LWFBrook90R and vice versa.
# Please use this translation function with caution and manually check the generated output.
"""
    prepare_for_LWFBrook90R(s::DiscretizedSPAC; return_value = "inputs", out_dir = ".", file_suffix = "")

Transforms a SPAC model to input for the R package LWFBrook90R. Its return value depends on the argument return_value ∈ ["results", "inputs", "csv_paths"]
The options are:
- "results":    run LWFBrook90R and return the simulation results in Julia without writing to disk
- "inputs":     return the input as DataFrames in Julia without writing to disk
- "csv_paths":  write the input as .csv files to disk and return the paths ()
the current working directory for running with R.

Note that if `return_value = "results"` a working installation of R and
LWFBrook90R v0.5.3 is required.
"""
function prepare_for_LWFBrook90R(spac_; return_value = "inputs", out_dir = ".", file_suffix = "")
    @assert spac_ isa DiscretizedSPAC && LWFBrook90.is_setup(spac_.parametrizedSPAC) "Please provide a `discretized` SPAC, generated with `setup(SPAC)`."
    @assert return_value ∈ ["results", "inputs", "csv_paths"]
    parSPAC = spac_.parametrizedSPAC

    # Variant 1): prepare arguments to run_LWFB90()
        # run_LWFB90(options_b90 = opts_target,
        #                          param_b90 = parms_target,
        #                          climate = meteo_target,
        #                          soil = soil_target)
    # NOT IMPLEMENTED Variant 2): prepare arguments to r_lwfbrook90()
        # r_lwfbrook90(siteparam, climveg, param, pdur, soil_materials, soil_nodes,
        #              precdat = NULL, output_log = TRUE, chk_input = TRUE, timelimit = Inf)
        # r_lwfBrook90Args = list(
        #   siteparam = data.frame(simyears[1],
        #                          as.integer(format(options_b90$startdate, "%j")),
        #                          param_b90$coords_y, param_b90$snowini, param_b90$gwatini,
        #                          options_b90$prec_interval),
        #   climveg = cbind(climate[, c("yr", "mo", "da", "globrad", "tmax", "tmin", "vappres",
        #                               "windspeed", "prec", "mesfl")],
        #                   standprop_daily[, c("densef", "height", "lai", "sai", "age")]),
        #   precdat = precip[, c("yr", "mo", "da", "ii", "prec", "mesfl")],
        #   param = param_to_rlwfbrook90(param_b90, options_b90$imodel),
        #   pdur = param_b90$pdur,
        #   soil_materials = param_b90$soil_materials,
        #   soil_nodes = param_b90$soil_nodes[, c("layer", "midpoint", "thick", "mat", "psiini", "rootden")],
        #   output_log = verbose,
        #   chk_input = chk_input,
        #   timelimit = timelimit)

        # NOT IMPLEMENTED Variant 2): prepare arguments to r_lwfbrook90()
                    # # r_lwfBrook90Args$siteparam
                    # simyears.1. as.integer.format.options_b90.startdate....j... param_b90.coords_y param_b90.snowini param_b90.gwatini options_b90.prec_interval
                    # 1        2015                                               1             46.583                 0                 0                         1
                    # # r_lwfBrook90Args$siteparam
                    # simyears.1. as.integer.format.options_b90.startdate....j... param_b90.coords_y param_b90.snowini param_b90.gwatini options_b90.prec_interval
                    # 1        2015                                               1             46.583                 0                 0                         1
                    # # r_lwfBrook90Args$climveg
                    # yr mo da globrad tmax tmin vappres windspeed prec mesfl densef height lai sai      age
                    # 1  2015  1  1    6.31 -2.1 -7.4    0.38       2.0  0.0     0      1     40   0   1 200.0000
                    # 2  2015  1  2    3.69  4.1 -4.7    0.56       2.5  0.7     0      1     40   0   1 200.0027
                    # 3  2015  1  3    1.51  7.6 -0.1    0.76       5.0  5.0     0      1     40   0   1 200.0055
                    # 61 2015  3  2    8.01  7.0  3.0    0.73       6.3  7.5     0      1     40   0   1 200.1648
                    # 62 2015  3  3   14.34  7.4  0.5    0.58       2.9  3.0     0      1     40   0   1 200.1676
                    # 63 2015  3  4    7.14  4.3 -0.3    0.54       3.2  1.5     0      1     40   0   1 200.1703
                    # 64 2015  3  5   15.91  2.9 -1.9    0.42       7.0  0.0     0      1     40   0   1 200.1731
                    # 65 2015  3  6   16.48  5.4 -1.9    0.42       5.5  0.0     0      1     40   0   1 200.1758
                    # 66 2015  3  7   15.46  8.5 -1.8    0.50       3.1  0.0     0      1     40   0   1 200.1786
                    # [ reached 'max' / getOption("max.print") -- omitted 2856 rows ]
                    # # r_lwfBrook90Args$precdat
                    # NULL
                    # # r_lwfBrook90Args$param
                    # [1] 2922.00000    0.00000    4.01000   45.00000    0.20000    0.50000    0.25000    0.50000    0.20000    0.30000 5000.00000    0.00500    2.00000    0.10000    0.00325    0.00100    4.00000    0.35000    0.13000    0.05000
                    # [21]    1.00000   10.00000    2.00000    2.00000    2.50000   -0.50000    0.00000    0.00000    0.06000    0.04000    0.06000    0.04000    0.15000    0.15000    0.60000    0.60000    1.50000    0.30000    0.20000    0.50000
                    # [41]    0.35000    0.05000    0.30000    0.30000    0.00530    0.50000    0.00030 1000.00000  100.00000    2.00000    0.00000   10.00000   30.00000   40.00000    8.00000 3000.00000   12.00000    0.25000    0.03000    0.00000
                    # [61]    0.50000   -2.00000    0.35000    1.00000   19.00000    3.00000    2.00000    0.00000    1.00000  100.00000    1.00000    0.00000    0.00000    1.00000    0.00000    0.00000    0.00000  200.00000    1.00000    0.00000
                    # [81]    0.00000    0.50000    0.05000    0.00050
                    # # r_lwfBrook90Args$soil_materials
                    # mat   ths thr  alpha  npar    ksat   tort gravel
                    # 1:   1 0.410   0 10.486 1.184 1412.95 -3.236   0.01
                    # 2:   2 0.410   0 10.486 1.184 1412.95 -3.236   0.05
                    # 3:   3 0.379   0 20.387 1.235 2854.91 -3.339   0.01
                    # # r_lwfBrook90Args$soil_nodes
                    # layer midpoint thick mat psiini      rootden
                    # 1:     1   -0.025    50   1     -6 0.0282531949
                    # 2:     2   -0.075    50   1     -6 0.0242619798
                    # 3:     3   -0.125    50   1     -6 0.0208345876
                    # 4:     4   -0.175    50   1     -6 0.0178913692
                    # 5:     5   -0.250   100   1     -6 0.0142787274
                    # 6:     6   -0.350   100   1     -6 0.0105294781
                    # 7:     7   -0.450   100   1     -6 0.0077646912
                    # 8:     8   -0.550   100   2     -6 0.0057258706
                    # 9:     9   -0.650   100   2     -6 0.0042223952
                    # 10:    10   -0.750   100   2     -6 0.0031136961
                    # 11:    11   -0.850   100   2     -6 0.0022961146
                    # 12:    12   -0.950   100   2     -6 0.0016932103
                    # 13:    13   -1.050   100   2     -6 0.0012486141
                    # 14:    14   -1.150   100   2     -6 0.0009207582
                    # 15:    15   -1.250   100   1     -6 0.0006789893
                    # 16:    16   -1.325    50   1     -6 0.0005387571
                    # 17:    17   -1.375    50   3     -6 0.0004626491
                    # 18:    18   -1.450   100   3     -6 0.0003692305
                    # 19:    19   -1.550   100   3     -6 0.0000000000
                    # # r_lwfBrook90Args$output_log
                    # [1] FALSE
                    # # r_lwfBrook90Args$chk_input
                    # [1] TRUE
                    # # r_lwfBrook90Args$timelimit
                    # [1] Inf

        # Variant 1): prepare arguments to run_LWFB90()
                    # # source: https://pschmidtwalter.github.io/LWFBrook90R/
                    # library(tibble)
                    # # load package and sample data
                    # library(LWFBrook90R)
                    # utils::packageVersion("LWFBrook90R")
                    # data(slb1_meteo, slb1_soil)

                    # # set up default model control options and parameters
                    # opts <- set_optionsLWFB90()
                    # parms <- set_paramLWFB90()

                    # # Derive soil hydraulic properties from soil physical properties
                    # # using a pedotransfer function:
                    # soil <- cbind(slb1_soil, hydpar_wessolek_tab(texture = slb1_soil$texture))

                    # # run the model and capture results
                    # lwfb90_res <- run_LWFB90(options_b90 = opts,
                    #                          param_b90 = parms,
                    #                          climate = slb1_meteo,
                    #                          soil = soil)
                    ###############
                    # dput(opts)
                    # # list(startdate        = structure(11688, class = "Date"),
                    # #      enddate          = structure(12417, class = "Date"),
                    # #      fornetrad        = "globrad",
                    # #      prec_interval    = 1,
                    # #      correct_prec     = FALSE,
                    # #      budburst_method  = "fixed",
                    # #      leaffall_method  = "fixed",
                    # #      standprop_input  = "parameters",
                    # #      standprop_interp = "constant",
                    # #      use_growthperiod = FALSE,
                    # #      lai_method       = "b90",
                    # #      imodel           = "MvG",
                    # #      root_method      = "betamodel")
                    # dput(parms)
                    # list(maxlai = 5, sai = 1, sai_ini = 1,
                    #      height = 25, height_ini = 25, densef = 1, densef_ini = 1, age_ini = 100,
                    #      standprop_table = NULL,
                    #      winlaifrac = 0,
                    #      budburst_species = "Fagus sylvatica", budburstdoy = 121, leaffalldoy = 279,
                    #      shp_budburst = 0.3, shp_leaffall = 3, shp_optdoy = 210,

                    #      emergedur = 28, leaffalldur = 58, lai_doy = NULL,
                    #      lai_frac = NULL, alb = 0.2, albsn = 0.5, ksnvp = 0.3, fxylem = 0.5,
                    #      mxkpl = 8, lwidth = 0.1, psicr = -2, nooutf = 1, lpc = 4,
                    #      cs = 0.035, czs = 0.13, czr = 0.05, hs = 1, hr = 10, rhotp = 2,
                    #      nn = 2.5, maxrlen = 3000, initrlen = 12, initrdep = 0.25,
                    #      rrad = 0.35, rgrorate = 0.03, rgroper = 0,

                    #      maxrootdepth = -1.5, betaroot = 0.97, rootden_table = NULL,

                    #      radex = 0.5, glmax = 0.0053, glmin = 3e-04, rm = 1000, r5 = 100, cvpd = 2, tl = 0,
                    #      t1 = 10, t2 = 30, th = 40, frintlai = 0.06, frintsai = 0.06, fsintlai = 0.04,
                    #      fsintsai = 0.04, cintrl = 0.15, cintrs = 0.15, cintsl = 0.6, cintss = 0.6,
                    #      infexp = 0, bypar = 0, qfpar = 1, qffc = 0, imperv = 0, drain = 1, gsc = 0, gsp = 0,
                    #      ilayer = 1, qlayer = 0,
                    #      z0s = 0.001, rstemp = -0.5, melfac = 1.5, ccfac = 0.3, laimlt = 0.2, saimlt = 0.5,
                    #      grdmlt = 0.35, maxlqf = 0.05, snoden = 0.3, obsheight = 0.025,
                    #      correct_prec_statexp = "mg", rssa = 100, rssb = 1,
                    #      soil_nodes = NULL, soil_materials = NULL,
                    #      dtimax = 0.5, dswmax = 0.05, dpsimax = 5e-04,
                    #      wndrat = 0.3, fetch = 5000, z0w = 0.005, zw = 2, zminh = 2,
                    #      coords_x = 9.9095, coords_y = 51.544,
                    #      c1 = 0.25, c2 = 0.5, c3 = 0.2, pdur = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4),
                    #      eslope = 0, aspect = 0, dslope = 0, slopelen = 200,
                    #      intrainini = 0, intsnowini = 0, gwatini = 0, snowini = 0, psiini = -6.3)
                    ###############

                    # tibble(slb1_meteo[1:10,])
                    # # > tibble(slb1_meteo[1:10,])
                    # # # A tibble: 10 × 9
                    # # dates       tmin  tmax tmean  prec relhum globrad windspeed vappres
                    # # <date>     <dbl> <dbl> <dbl> <dbl>  <int>   <dbl>     <dbl>   <dbl>
                    # #   1 1960-01-01   3.8   7.9   6.2 10.4      87   0.976       1.2   0.812
                    # # 2 1960-01-02   4.6   6.6   5.2  3.08     95   0.104       0.7   0.866
                    # # 3 1960-01-03   2.9   4.9   3.5 11.6     100   0.121       0.7   0.809
                    # # 4 1960-01-04   1.2   3.9   2.8  0        98   0.130       0.7   0.722
                    # # 5 1960-01-05   2.1   5.3   3.8 14.8      96   0.138       1.3   0.769
                    # # 6 1960-01-06  -0.1   3.2   1.3  5.00     95   0.397       1.1   0.653
                    # # 7 1960-01-07   1     3.1   2.7 10.6     100   0.164       0.8   0.710
                    # # 8 1960-01-08  -2.9   2.6  -0.6  1.81     93   0.657       0.9   0.572
                    # # 9 1960-01-09  -8.9  -2.8  -7.5  1.70     83   1.13        0.9   0.336
                    # # 10 1960-01-10 -11.1  -6    -8.5  4.04     85   2.36        0.9   0.277

                    ###############
                    # tibble(soil)
                    # # > tibble(soil)
                    # # # A tibble: 21 × 17
                    # # horizon     upper lower texture    bd gravel  sand  silt  clay c_org   ths    thr alpha  npar  mpar  ksat  tort
                    # # <chr>       <dbl> <dbl> <chr>   <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
                    # #   1 I Ah         0    -0.01 Ut3       1     0.04  11.2  74.6  14.2  9.55 0.403 0.0053  1.68  1.21 0.171  277. -1.20
                    # # 2 I Aeh       -0.01 -0.03 Ut3       1.1   0.04  11.2  74.6  14.2  6.62 0.403 0.0053  1.68  1.21 0.171  277. -1.20
                    # # 3 I Aeh       -0.03 -0.05 Ut3       1.1   0.04  11.2  74.6  14.2  6.62 0.403 0.0053  1.68  1.21 0.171  277. -1.20
                    # # 4 I Aeh       -0.05 -0.08 Ut3       1.1   0.04  11.2  74.6  14.2  6.62 0.403 0.0053  1.68  1.21 0.171  277. -1.20
                    # # 5 I (Bv)Sw-Bv -0.08 -0.12 Ut3       1.4   0.04  11.2  74.6  14.2  1.25 0.403 0.0053  1.68  1.21 0.171  277. -1.20
                    # # 6 I (Bv)Sw-Bv -0.12 -0.16 Ut3       1.4   0.04  11.2  74.6  14.2  1.25 0.403 0.0053  1.68  1.21 0.171  277. -1.20
                    # # 7 I (Bv)Sw-Bv -0.16 -0.22 Ut3       1.4   0.04  11.2  74.6  14.2  1.25 0.403 0.0053  1.68  1.21 0.171  277. -1.20
                    # # 8 I (Bv)Sw-Bv -0.22 -0.28 Ut3       1.4   0.04  11.2  74.6  14.2  1.25 0.403 0.0053  1.68  1.21 0.171  277. -1.20
                    # # 9 I (Bv)Sw-Bv -0.28 -0.36 Ut3       1.4   0.04  11.2  74.6  14.2  1.25 0.403 0.0053  1.68  1.21 0.171  277. -1.20
                    # # 10 I (Bv)Sw-Bv -0.36 -0.44 Ut3       1.4   0.04  11.2  74.6  14.2  1.25 0.403 0.0053  1.68  1.21 0.171  277. -1.20
                    # # # ℹ 11 more rows
                    # # # ℹ Use `print(n = ...)` to see more rows
                    ###############

    # A) opts:
    startdate, enddate = parSPAC.reference_date .+ Day.(parSPAC.tspan)
    enddate = enddate - Day(1)

    fornetrad        = "globrad"
    @assert unique(diff(parSPAC.forcing.meteo["p_days"])) == [1.0]
    prec_interval    = 1
    correct_prec     = false # hardcoded
    @assert !(parSPAC.pars.canopy_evolution.LAI_rel isa Vector) "LAI as vector is not supported."
    @assert parSPAC.pars.canopy_evolution.LAI_rel isa NamedTuple "LAI needs to be a NamedTuple of parameters."
    @assert keys(parSPAC.pars.canopy_evolution.LAI_rel) == (:DOY_Bstart, :Bduration, :DOY_Cstart, :Cduration, :LAI_perc_BtoC, :LAI_perc_CtoB)
    lai_method       = "b90"
    budburst_method  = "fixed"
    leaffall_method  = "fixed"
    standprop_input  = "parameters"
    standprop_interp = "constant"
    @assert eltype(parSPAC.pars.soil_horizons.shp) == LWFBrook90.KPT.MualemVanGenuchtenSHP
    imodel           = "MvG"
    use_growthperiod = false
    @assert parSPAC.pars.root_distribution isa NamedTuple "root distribution needs to be a NamedTuple of parameters."
    @assert keys(parSPAC.pars.root_distribution) == (:beta, :z_rootMax_m)
    root_method      = "soilvar"
    opts = (;
        startdate        = startdate,
        enddate          = enddate,
        fornetrad        = fornetrad,
        prec_interval    = prec_interval,
        correct_prec     = correct_prec,
        budburst_method  = budburst_method,
        leaffall_method  = leaffall_method,
        standprop_input  = standprop_input,
        standprop_interp = standprop_interp,
        use_growthperiod = use_growthperiod,
        lai_method       = lai_method,
        imodel           = imodel,
        root_method      = root_method)

    # B) parms:
    paramsjl = parSPAC.pars.params
    # warning(sprintf("Note that LWFBrook90R has a parameter 'cs' (ratio of height to SAI) of %s. LWFBrook90.jl uses a hardcoded default value of 0.035."))
    # warning(sprintf("Note that LWFBrook90R has a parameter 'correct_prec_statexp' (), that does not appear in LWFBrook90.jl."))
    @assert parSPAC.pars.IC_scalar[1,:u_CC_init_MJ_per_m2] <= 0.0001
    @assert parSPAC.pars.IC_scalar[1,:u_SNOWLQ_init_mm] <= 0.0001

    parms = (;
        maxlai = paramsjl.MAXLAI,
        sai = paramsjl.SAI_baseline_,
        sai_ini = paramsjl.SAI_baseline_,

        height     = paramsjl.HEIGHT_baseline_m,
        height_ini = paramsjl.HEIGHT_baseline_m,
        densef     = paramsjl.DENSEF_baseline_,
        densef_ini = paramsjl.DENSEF_baseline_,
        age_ini    = paramsjl.AGE_baseline_yrs,
        winlaifrac = parSPAC.pars.canopy_evolution.LAI_rel.LAI_perc_CtoB / parSPAC.pars.canopy_evolution.LAI_rel.LAI_perc_BtoC,
        budburst_species = missing,#"Fagus sylvatica",
        budburstdoy = parSPAC.pars.canopy_evolution.LAI_rel.DOY_Bstart,
        leaffalldoy = parSPAC.pars.canopy_evolution.LAI_rel.DOY_Cstart,
        emergedur = parSPAC.pars.canopy_evolution.LAI_rel.Bduration,
        leaffalldur = parSPAC.pars.canopy_evolution.LAI_rel.Cduration,
        shp_optdoy = missing,   #parSPAC.pars.canopy_evolution.LAI_rel.DOY_Bstart + parSPAC.pars.canopy_evolution.LAI_rel.Bduration,
        shp_budburst = missing,
        shp_leaffall = missing,

        alb = paramsjl.ALB, albsn = paramsjl.ALBSN, ksnvp = paramsjl.KSNVP, fxylem = paramsjl.FXYLEM,
        mxkpl = paramsjl.MXKPL, lwidth = paramsjl.LWIDTH, psicr = paramsjl.PSICR, nooutf = paramsjl.NOOUTF,
        lpc = paramsjl.LPC, cs = 0.35, czs = paramsjl.CZS, czr = paramsjl.CZR, hs = paramsjl.HS,
        hr = paramsjl.HR, rhotp = paramsjl.RHOTP, nn = paramsjl.NN, maxrlen = paramsjl.MXRTLN,
        initrlen = paramsjl.INITRLEN, initrdep = paramsjl.INITRDEP,
        rrad = paramsjl.RTRAD, rgrorate = paramsjl.RGRORATE, rgroper = paramsjl.RGROPER,
        # maxrootdepth = -1.5, betaroot = 0.97,
        maxrootdepth = missing, betaroot = missing,
        radex = paramsjl.CR, glmax = paramsjl.GLMAX, glmin = paramsjl.GLMIN, rm = paramsjl.RM, r5 = paramsjl.R5, cvpd = paramsjl.CVPD,  tl = paramsjl.TL,
        t1 = paramsjl.T1, t2 = paramsjl.T2, th = paramsjl.TH, frintlai = paramsjl.FRINTLAI, frintsai = paramsjl.FRINTSAI, fsintlai = paramsjl.FSINTLAI, fsintsai = paramsjl.FSINTSAI,
        cintrl = paramsjl.CINTRL, cintrs = paramsjl.CINTRS, cintsl = paramsjl.CINTSL, cintss = paramsjl.CINTSS,
        infexp = paramsjl.INFEXP, bypar = paramsjl.BYPAR, qfpar = paramsjl.QFPAR, qffc = paramsjl.QFFC, imperv = paramsjl.IMPERV, drain = paramsjl.DRAIN, gsc = paramsjl.GSC, gsp = paramsjl.GSP,
        ilayer = sum(paramsjl.IDEPTH_m .>= cumsum(parSPAC.soil_discretization.Δz)),
        qlayer = sum(paramsjl.QDEPTH_m .>= cumsum(parSPAC.soil_discretization.Δz)),
        z0s = paramsjl.Z0S, rstemp = paramsjl.RSTEMP, melfac = paramsjl.MELFAC, ccfac = paramsjl.CCFAC, laimlt = paramsjl.LAIMLT, saimlt = paramsjl.SAIMLT,
        grdmlt = paramsjl.GRDMLT, maxlqf = paramsjl.MAXLQF, snoden = paramsjl.SNODEN,
        obsheight = paramsjl.Z0G / paramsjl.CZS,

        correct_prec_statexp = "mg",
        rssa = paramsjl.RSSA, rssb = paramsjl.RSSB,
        dtimax  = parSPAC.solver_options.DTIMAX,
        dswmax  = parSPAC.solver_options.DSWMAX,
        dpsimax = parSPAC.solver_options.DPSIMAX,
        wndrat =paramsjl.WNDRAT, fetch =paramsjl.FETCH, z0w =paramsjl.Z0W, zw =paramsjl.ZW, zminh =paramsjl.ZMINH,
        coords_x = 0, # has no effect on sim results
        coords_y = paramsjl.LAT_DEG,
        c1 = paramsjl.C1, c2 = paramsjl.C2, c3 = paramsjl.C3, pdur = parSPAC.forcing.storm_durations.storm_durations_h,
        eslope = paramsjl.ESLOPE_DEG, aspect = paramsjl.ASPECT_DEG, dslope = paramsjl.DSLOPE, slopelen = paramsjl.LENGTH_SLOPE,
        intrainini = parSPAC.pars.IC_scalar[1,:u_INTR_init_mm], intsnowini = parSPAC.pars.IC_scalar[1,:u_INTS_init_mm], gwatini = parSPAC.pars.IC_scalar[1,:u_GWAT_init_mm],
        snowini = parSPAC.pars.IC_scalar[1,:u_SNOW_init_mm],
        psiini = parSPAC.pars.IC_soil.PSIM_init_kPa,
        # in R these are NULL
        standprop_table = missing,
        lai_doy = missing, lai_frac = missing,
        rootden_table = missing,
        soil_nodes = missing, soil_materials = missing,
    )

    # C) meteo
    datetimes = LWFBrook90.RelativeDaysFloat2DateTime.(collect(parSPAC.forcing.meteo["p_days"]), parSPAC.reference_date)
    tmin  = parSPAC.forcing.meteo["p_TMIN"].itp.itp.coefs
    tmax  = parSPAC.forcing.meteo["p_TMAX"].itp.itp.coefs
    meteo = DataFrame(
        dates = Date.(datetimes),
        tmin  = tmin,
        tmax  = tmax,
        tmean = NaN, #(tmin .+ tmax)/2,
        prec = parSPAC.forcing.meteo["p_PREC"].itp.itp.coefs,
        relhum_dummy = 100, #parSPAC.forcing.meteo["p_VAPPRES"],
        globrad = parSPAC.forcing.meteo["p_GLOBRAD"].itp.itp.coefs,
        windspeed = parSPAC.forcing.meteo["p_WIND"].itp.itp.coefs,
        vappres = parSPAC.forcing.meteo["p_VAPPRES"].itp.itp.coefs,
    )

    # D) soil
    shps = parSPAC.soil_discretization.df.shp
    sd_df = parSPAC.soil_discretization.df
    soil = DataFrame(
        horizon = "NA",
        upper   = sd_df.Upper_m,
        lower   = sd_df.Lower_m,
        texture = "NA",
        bd      = missing,
        gravel  = [shp.p_STONEF for shp in shps],
        sand    = missing,
        silt    = missing,
        clay    = missing,
        c_org   = missing,
        ths     = [shp.p_THSAT for shp in shps],
        thr     = [shp.p_θr for shp in shps],
        alpha   = [shp.p_MvGα for shp in shps],
        npar    = [shp.p_MvGn for shp in shps],
        mpar    = [shp.p_MvGm for shp in shps],
        ksat    = [shp.p_KSAT for shp in shps],
        tort    = [shp.p_MvGl for shp in shps],

        rootden = parSPAC.soil_discretization.df.Rootden_)

    @assert opts.root_method == "soilvar"
    @assert "rootden" ∈ names(soil)

    if return_value ∈ ["inputs", "csv_paths"] # return_value ∈ ["results", "inputs", "csv_paths"]
        # prepare DataFrames out of NamedTuples `parms` and `opts`
        df_parms_a = DataFrame([k => v for (k,v) in pairs(parms) if k != :pdur]) # single row
        df_parms_b = DataFrame([k => v for (k,v) in pairs(parms) if k == :pdur]) # 12 rows
        # df_parms = DataFrame([k => v for (k,v) in pairs(parms)])
        df_opts  = DataFrame([k => v for (k,v) in pairs(opts)])

        if return_value == "inputs"
            return (args    = (meteo = meteo, soil = soil, parms = parms, opts = opts),
                    df_args = (meteo = meteo, soil = soil, df_parms_a = df_parms_a, df_parms_b = df_parms_b, df_opts = df_opts))
        else # "csv_paths"
            mkpath(out_dir)
            # R"print('Hello from the other side!')"
            # R"saveRDS($parms,paste0($out_dir, '/parms',$file_suffix,'.RDS'));" # keeping NamedTuple to directly import in R
            # R"saveRDS($opts, paste0($out_dir, '/opts', $file_suffix,'.RDS'));" # keeping NamedTuple to directly import in R
            # CSV.write("meteo.csv", meteo)
            # CSV.write("soil.csv", soil)
            # println("Successfully wrote: 'parms.RDS', 'opts.RDS', 'meteo.csv', 'soil.csv'.")

            meteo_path      = joinpath(out_dir, "meteo$(file_suffix).csv")      ; CSV.write(meteo_path,      meteo)
            soil_path       = joinpath(out_dir, "soil$(file_suffix).csv")       ; CSV.write(soil_path,       soil)
            df_opts_path    = joinpath(out_dir, "df_opts$(file_suffix).csv")    ; CSV.write(df_opts_path,    df_opts)
            df_parms_a_path = joinpath(out_dir, "df_parms_a$(file_suffix).csv") ; CSV.write(df_parms_a_path, df_parms_a)
            df_parms_b_path = joinpath(out_dir, "df_parms_b$(file_suffix).csv") ; CSV.write(df_parms_b_path, df_parms_b)

            println("Sucessfully wrote: '$(meteo_path)'.",)
            println("Sucessfully wrote: '$(soil_path)'.",)
            println("Sucessfully wrote: '$(df_opts_path)'.",)
            println("Sucessfully wrote: '$(df_parms_a_path)'.",)
            println("Sucessfully wrote: '$(df_parms_b_path)'.")
            return (; df_paths = (meteo = meteo_path, soil = soil_path, df_parms_a = df_parms_a_path, df_parms_b = df_parms_b_path, df_opts = df_opts_path))
        end
    else # "results"
        # TODO: this is working, but currently deactivated, as this requires RCall in the dependencies

        # # run LWFBrook90R and return results
        # # run directly without files
        # R"""
        # library(LWFBrook90R)
        # library(readr)
        # library(lubridate)
        # library(dplyr)

        # # check for version
        # stopifnot(packageVersion("LWFBrook90R") == "0.5.3")

        # # a) load from Julia RCall
        # meteo <- $meteo
        # soil <- $soil
        # parms <- $parms
        # opts <- $opts
        # # b) load from files
        #     # meteo.csv.path  = "meteo.csv"
        #     # soil.csv.path  = "soil.csv"
        #     # opts.RDS.path  = "opts.RDS"
        #     # parms.RDS.path = "parms.RDS"
        #     # parms <- readRDS(parms.RDS.path)
        #     # opts  <- readRDS(opts.RDS.path)
        #     # meteo <- readr::read_csv(meteo.csv.path)
        #     # soil  <- readr::read_csv(soil.csv.path)

        # # Fix some errors in output
        # opts$startdate <- lubridate::ymd(as.character(opts$startdate))
        # opts$enddate   <- lubridate::ymd(as.character(opts$enddate))

        # sim_out <- run_LWFB90(options_b90 = opts,
        #                     param_b90   = parms,
        #                     climate     = meteo,
        #                     soil        = soil,
        #                     rtrn_input = TRUE)

        # output <- sim_out$output
        # layer_output <- sim_out$layer_output

        # depths_to_read_out <- -c(150,500,800)/1000
        # idx_to_read_out <- c()
        # for (it in seq_along(depths_to_read_out)) {
        # idx_to_read_out[it] <- which.min(abs(depths_to_read_out[it] - sim_out$model_input$param_b90$soil_nodes$midpoint))
        # }
        # depths_actually_read_out <-
        # sim_out$model_input$param_b90$soil_nodes$midpoint[idx_to_read_out]

        # theta_specific_depths <- {select(sim_out$layer_output, yr, mo, da, nl, theta) %>%
        #         tidyr::unite(col = "date", yr, mo, da) %>%
        #         mutate(date = lubridate::ymd(date)) %>%
        #         left_join({sim_out$model_input$param_b90$soil_nodes %>%
        #             mutate(lower = -0 + -cumsum(thick)) %>%
        #             mutate(upper = c(0, lower[-n()])) %>%
        #             select(nl = layer, midpoint, upper, lower)},
        #             by = "nl") }
        # psi_specific_depths <- {select(sim_out$layer_output, yr, mo, da, nl, psimi) %>%
        #         tidyr::unite(col = "date", yr, mo, da) %>%
        #         mutate(date = lubridate::ymd(date)) %>%
        #         left_join({sim_out$model_input$param_b90$soil_nodes %>%
        #             mutate(lower = -0 + -cumsum(thick)) %>%
        #             mutate(upper = c(0, lower[-n()])) %>%
        #             select(nl = layer, midpoint, upper, lower)},
        #             by = "nl") }
        # """

        # @rget output
        # @rget layer_output
        # @rget psi_specific_depths
        # @rget theta_specific_depths

        # return (; output = output, layer_output = layer_output,
        #         psi_specific_depths = psi_specific_depths, theta_specific_depths = theta_specific_depths)
    end
end

##########################
# Below commented out lines show use cases of the code above:

        # # from run_example()
        # input_prefix = "DAV2020-full";
        # input_path   = "examples/"*input_prefix;
        # model        = loadSPAC(input_path, input_prefix; simulate_isotopes = true);
        # base_simulation   = LWFBrook90.setup(model; requested_tspan = (0,300));

        # Δz = [fill(0.05, 4);  fill(0.10, 14)]
        # mod_model = loadSPAC(input_path, input_prefix;
        #         simulate_isotopes = false,
        #         Δz_thickness_m    = Δz,
        #         root_distribution = (beta = 0.98, z_rootMax_m = -sum(Δz)), # use whole domain as root zone (beta parameter regulates distribution)
        #         IC_soil           = (PSIM_init_kPa = -6.0,
        #                             delta18O_init_permil = -9.0,
        #                             delta2H_init_permil = -11.0),
        #         canopy_evolution  = (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel    = 100,
        #                             LAI_rel = (DOY_Bstart = 110,    Bduration  = 20,
        #                                         DOY_Cstart = 270,    Cduration  = 60,
        #                                         LAI_perc_BtoC = 100, LAI_perc_CtoB = 0)),
        #         storm_durations_h = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
        #         IC_scalar         = (amount = (u_GWAT_init_mm = 0.0,         u_INTS_init_mm = 0.0,
        #                                     u_INTR_init_mm = 0.0,         u_SNOW_init_mm = 0.0,
        #                                     u_CC_init_MJ_per_m2 = 0.0001, u_SNOWLQ_init_mm =  0.),
        #                             d18O    = (u_GWAT_init_permil = -11.111, u_INTS_init_permil = -12.222,
        #                                     u_INTR_init_permil = -13.333, u_SNOW_init_permil = -14.444),
        #                             d2H     = (u_GWAT_init_permil = -95.111, u_INTS_init_permil = -95.222,
        #                                     u_INTR_init_permil = -95.333, u_SNOW_init_permil = -95.444)));
        # mod_simulation = LWFBrook90.setup(mod_model)
        # curr_simulation = mod_simulation




        # # # from LWFBrook90.jl-Calibration:
        # # using DataFrames
        # # simulate_isotopes = true
        # # Δz = [fill(0.05, 4);  fill(0.10, 14)]
        # # input_path = "../../../../LWF-Brook90.jl-calibration/pgm-calibrate-HYPERION/input-2023-10/";
        # # input_path = "../../../LWF-Brook90.jl-calibration/pgm-calibrate-HYPERION/input-2023-10/";

        # # input_prefix = "LAU-2023-10-input"
        # # # generate base model (by specifiying options here, less input files are needed):
        # # loadSPAC_args = (input_path = input_path, input_prefix = input_prefix,
        # #         simulate_isotopes = simulate_isotopes,
        # #         Δz_thickness_m    = Δz,
        # #         root_distribution = (beta = 0.98, z_rootMax_m = -sum(Δz)), # use whole domain as root zone (beta parameter regulates distribution)
        # #         IC_soil           = (PSIM_init_kPa = -6.0,
        # #                             delta18O_init_permil = -9.0,
        # #                             delta2H_init_permil = -11.0),
        # #         canopy_evolution  = (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel    = 100,
        # #                             LAI_rel = (DOY_Bstart = 110,    Bduration  = 20,
        # #                                         DOY_Cstart = 270,    Cduration  = 60,
        # #                                         LAI_perc_BtoC = 100, LAI_perc_CtoB = 0)),
        # #         storm_durations_h = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
        # #         IC_scalar         = (amount = (u_GWAT_init_mm = 0.0,         u_INTS_init_mm = 0.0,
        # #                                     u_INTR_init_mm = 0.0,         u_SNOW_init_mm = 0.0,
        # #                                     u_CC_init_MJ_per_m2 = 0.0001, u_SNOWLQ_init_mm =  0.),
        # #                             d18O    = (u_GWAT_init_permil = -11.111, u_INTS_init_permil = -12.222,
        # #                                     u_INTR_init_permil = -13.333, u_SNOW_init_permil = -14.444),
        # #                             d2H     = (u_GWAT_init_permil = -95.111, u_INTS_init_permil = -95.222,
        # #                                     u_INTR_init_permil = -95.333, u_SNOW_init_permil = -95.444)));
        # # base_model = loadSPAC(
        # #     loadSPAC_args[:input_path],
        # #     loadSPAC_args[:input_prefix];
        # #     loadSPAC_args[Not([:input_path, :input_prefix])]...);
        # # base_model_tspan_dates = LWFBrook90.RelativeDaysFloat2DateTime.(base_model.tspan, base_model.reference_date)

        # # base_simulation = setup(base_model)#, requested_tspan = simulation_tspan_dates);
        # # # simulation_tspan       = LWFBrook90.DateTime2RelativeDaysFloat.(simulation_tspan_dates, base_model.reference_date)
        # # # remakeSPAC(base_simulation; NamedTuple(k => v for (k, v) in pairs(opt_pars[sim_id]))...)
        # # base_simulation

        # curr_simulation = base_simulation

        # # variant a) get paths to input data as csv:
        # res = prepare_for_LWFBrook90R(curr_simulation,
        #     return_value = "csv_paths",# ∈ ["results", "inputs", "csv_paths"],
        #     out_dir = "generated_LWFBrook90R_input_files",
        #     file_suffix = "_specific_filenname_");

        # # a)
        # res.df_paths.meteo
        # res.df_paths.soil
        # res.df_paths.df_opts
        # res.df_paths.df_parms_a
        # res.df_paths.df_parms_b

        # # variant b) get input data as DataFrames in Julia:
        # args, derived_args = prepare_for_LWFBrook90R(curr_simulation,
        #     return_value = "inputs");

        # # b1)
        # args.meteo
        # args.soil
        # args.opts
        # args.parms
        # # b2)
        # derived_args.meteo
        # derived_args.soil
        # derived_args.df_opts
        # derived_args.df_parms_a
        # derived_args.df_parms_b


        # # variant c) get results from running R-package
        # R_res = prepare_for_LWFBrook90R(curr_simulation, return_value = "results");

        # # c)
        # R_res.output
        # R_res.layer_output
        # R_res.psi_specific_depths
        # R_res.theta_specific_depths




        #     # Make plots directly in R #####################
        #     # R"""
        #     # library(ggplot2)
        #     # select(sim_out$output, yr, mo, da, doy, rfal, sfal, #prec, # rfal+sfal instead of prec # rthr,
        #     #          flow, evap, seep, snow, swat, gwat, intr, ints, sint, isvp) %>%
        #     #   select(-doy) %>%
        #     #   tidyr::unite(col = "date", yr, mo, da) %>%
        #     #   mutate(date = lubridate::ymd(date)) %>%
        #     #   tidyr::pivot_longer(cols = -date,
        #     #                       # keep facets in ggplot in the same order as in data.frame
        #     #                       names_transform = list(name = forcats::fct_inorder)) %>%
        #     #   ggplot(aes(x=date, y=value, color = name)) +
        #     #   geom_line() + theme_bw() + facet_grid(name~., scales = "free_y") +
        #     #   theme(axis.text.y = element_text(size=5), legend.key.height = unit(1,"line")) +
        #     #   labs(x=NULL)
        #     # """

        #     # R"""
        #     # library(ggplot2)
        #     # library(dplyr)
        #     # log_modulus_trans <- function(base = exp(1)) {
        #     #   scales::trans_new(name = "log_modulus",
        #     #                     transform = function(x) sign(x) * log(abs(x) + 1, base = base),
        #     #                     inverse   = function(x) sign(x) * ( base^(abs(x)) - 1 ))}
        #     # log_modulus_rev_trans <- function(base = exp(1)) {
        #     #   scales::trans_new(name = "log_modulus_rev",
        #     #                     transform = function(x) -sign(x) * log(abs(x) + 1, base = base),
        #     #                     inverse   = function(x) -sign(x) * ( base^(abs(x)) - 1 ))}
        #     # scale_y_logModulus <- function(...) {
        #     #   scale_y_continuous(trans = "log_modulus",
        #     #                      labels=function(n) format(n, scientific=FALSE),
        #     #                      breaks = c((+1)*(10^(0:10) %*% t(c(1, 3)) ),
        #     #                                 0,
        #     #                                 (-1)*(10^(0:10) %*% t(c(1, 3)) )),
        #     #                      minor_breaks = c((+1)*(10^(-0:10) %*% t(c(1:9)) ),
        #     #                                       (-1)*(10^(-0:10) %*% t(c(1:9)) )),
        #     #                      ...)}
        #     # {select(sim_out$layer_output, yr, mo, da, nl, psimi) %>%
        #     #     tidyr::unite(col = "date", yr, mo, da) %>%
        #     #     mutate(date = lubridate::ymd(date)) %>%
        #     #     left_join({sim_out$model_input$param_b90$soil_nodes %>%
        #     #         mutate(lower = -0 + -cumsum(thick)) %>%
        #     #         mutate(upper = c(0, lower[-n()])) %>%
        #     #         select(nl = layer, midpoint, upper, lower)},
        #     #         by = "nl") } %>%
        #     # ggplot(aes(x=date, y=midpoint, fill = psimi)) +
        #     # geom_rect(aes(xmin=date, xmax=date+1, ymin=lower, ymax=upper)) +
        #     # theme_bw() + labs(fill = "ψ_m (kPa)", y = "Depth (m)", x=NULL) +
        #     # # define color palette separating saturated from unsaturated and using log_modulus transform
        #     # scale_fill_gradientn(
        #     # # sharp break between blue (unsaturated) and red (saturated) at 0:
        #     # # https://github.com/tidyverse/ggplot2/issues/3738#issuecomment-1336166757
        #     # # colours = c("white", scales::muted("blue"), scales::muted("red"), "black"), values = c(0, 0.5, 0.5001, 1.0),
        #     # colours = c("white", scales::muted("blue"), "red", "black"), values = c(0, 0.4999, 0.5001, 1.0), # This highlights saturation = 1 in read already
        #     # rescaler = ~ scales::rescale_mid(.x, mid = 0),
        #     # trans ="log_modulus",
        #     # labels=function(n) format(n, scientific=FALSE),
        #     # breaks = c((+1)*(10^(0:10) %*% t(c(1, 3)) ),
        #     #             0,
        #     #             (-1)*(10^(0:10) %*% t(c(1, 3)) )),
        #     # minor_breaks = c((+1)*(10^(-0:10) %*% t(c(1:9)) ),
        #     #                     (-1)*(10^(-0:10) %*% t(c(1:9)) )))
        #     # """

        #     # R"""
        #     # depths_to_read_out <- -c(150,500,800)/1000
        #     # idx_to_read_out <- c()
        #     # for (it in seq_along(depths_to_read_out)) {
        #     # idx_to_read_out[it] <- which.min(abs(depths_to_read_out[it] - sim_out$model_input$param_b90$soil_nodes$midpoint))
        #     # }
        #     # depths_actually_read_out <-
        #     # sim_out$model_input$param_b90$soil_nodes$midpoint[idx_to_read_out]

        #     # pl2_belowground_theta <- {select(sim_out$layer_output, yr, mo, da, nl, theta) %>%
        #     #     tidyr::unite(col = "date", yr, mo, da) %>%
        #     #     mutate(date = lubridate::ymd(date)) %>%
        #     #     left_join({sim_out$model_input$param_b90$soil_nodes %>%
        #     #         mutate(lower = -0 + -cumsum(thick)) %>%
        #     #         mutate(upper = c(0, lower[-n()])) %>%
        #     #         select(nl = layer, midpoint, upper, lower)},
        #     #         by = "nl") } %>%
        #     # filter(nl %in% c(idx_to_read_out)) %>%
        #     # ggplot(aes(x=date, y=theta, color = as.factor(midpoint))) +
        #     # geom_line() + theme_bw() +
        #     # theme_bw() + labs(y = "θ (-)", color = "Depth (m)", x=NULL) + theme(legend.position = c(0.87, 0.25))
        #     # pl2_belowground_psi<- {select(sim_out$layer_output, yr, mo, da, nl, psimi) %>%
        #     #     tidyr::unite(col = "date", yr, mo, da) %>%
        #     #     mutate(date = lubridate::ymd(date)) %>%
        #     #     left_join({sim_out$model_input$param_b90$soil_nodes %>%
        #     #         mutate(lower = -0 + -cumsum(thick)) %>%
        #     #         mutate(upper = c(0, lower[-n()])) %>%
        #     #         select(nl = layer, midpoint, upper, lower)},
        #     #         by = "nl") } %>%
        #     # filter(nl %in% c(idx_to_read_out)) %>%
        #     # ggplot(aes(x=date, y=psimi, color = as.factor(midpoint))) +
        #     # geom_line() + theme_bw() +
        #     # theme_bw() + labs(y = "ψ_m (kPa)", color = "Depth (m)", x=NULL) + theme(legend.position = c(0.87, 0.25)) +
        #     # # use log modulus scale
        #     # scale_y_logModulus()
        #     # # gridExtra::grid.arrange(pl2_belowground_theta, pl2_belowground_psi)
        #     # pl2_belowground_theta
        #     # pl2_belowground_psi
        #     # """
#################