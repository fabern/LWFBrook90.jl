#################
# Author: fabian.bernhard
# Date: 2021-02-11
# Task: Generate input files for LWFBrook90.jl
# Description: This code defines an R function that takes the same arguements as
#              LWFBrook90R::run_LWFB90() and generates the corresponding Julia
#              input functions.
#
#              If desired LWFBrook90R can also be run and output stored in csv
#              files for comparison with LWFBrook90.
#################

# Note that currently only basic functions are supported by LWFBrook90.


####### Define function
generate_LWFBrook90Julia_Input <- function(Julia_target_dir = NA,
                                           Julia_prefix = NA,
                                           run_LWFBrook90R = FALSE,
                                           ...){

  if(is.na(Julia_target_dir)) {stop("Please provide a target directory `Julia_target_dir`.")}
  if(is.na(Julia_prefix))     {stop("Please provide a prefix for input files `Julia_prefix`.")}

  require(dplyr)

  original_arguments <- list(...)
  # generate input_data by running run_LWFB90() with argument run = FALSE
  input_Data <- LWFBrook90R::run_LWFB90(..., run = FALSE)

  # A) time-dependent parameters: (meteo and stand properties)
  veg_daily   <- input_Data$standprop_daily %>% select(dates, densef, height, lai, sai, age)
  clim_daily  <- original_arguments$climate %>%
    select(dates, globrad, tmax, tmin, vappres, windspeed, prec) %>%
    filter(dates >= input_Data$options_b90$startdate) %>%
    filter(dates <= input_Data$options_b90$enddate)

  out_csv_meteoveg <- clim_daily %>%
    dplyr::left_join(veg_daily, by="dates") %>%
    select("dates",
           "globrad_MJDayM2" = globrad,
           "tmax_degC" = tmax,
           "tmin_degC" = tmin,
           "vappres_kPa" = vappres,
           "windspeed_ms" = windspeed,
           "prec_mmDay" = prec,
           "densef_" = densef,
           "height_m" = height,
           "lai_" = lai,
           "sai_" = sai,
           "age_yrs" = age)

  # B) time-independent parameters:
  # B1) site climate data
  out_pdur <- input_Data$param_b90$pdur
  out_csv_pdur <- data.frame(month = c("### Average durations of a single storm event for each month -------",
                                       "January", "Februrary", "March", "April",
                                       "May", "June", "July", "August",
                                       "September", "October", "November", "December"),
                             average_storm_duration_h = c(NA,out_pdur))
  # out_csv_precdat <- NULL # TODO(bernhard): currently not implemented

  # B2) i)  soil discretization and parameters and initial conditions of continuous quantities
  out_csv_soil_discretization <- input_Data$param_b90$soil_nodes %>%
    select(#SimulationNode = layer,
           Upper_m = upper, Lower_m = lower,
           Rootden = rootden, Psiini_kPa = psiini) %>%
    mutate(delta18O_mUr = NA, delta2H_mUr = NA)

  # B2) ii) soil horizon parameters in discrete horizons (different from numerical node discretization)
  out_csv_soil_horizons <-left_join(
      select(input_Data$param_b90$soil_nodes, HorizonNr = layer, Upper_m = upper, Lower_m = lower, mat),
      select(input_Data$param_b90$soil_materials,
             mat,
             ths_volFrac = ths,
             thr_volFrac = thr,
             alpha_perMeter = alpha,
             npar_ = npar,
             #mpar_ = mpar,
             ksat_mmDay = ksat,
             tort_ = tort,
             gravel_volFrac = gravel),
      by = "mat"
    ) %>% select(-mat)



  # B3) initial conditions of scalar state variables
  out_initial_conditions <- with(input_Data$param_b90,c(
    "### Initial conditions (except for depth-dependent u_aux_PSIM) -------" = NA,
    "u_GWAT_init"=gwatini,
    "u_INTS_init"=intsnowini,
    "u_INTR_init"=intrainini,
    "u_SNOW_init"=snowini,
    "u_CC_init"    =0,
    "u_SNOWLQ_init"=0
  ))
  out_csv_initial_conditions <- data.frame(param_id = names(out_initial_conditions),
                                           x        = unname(out_initial_conditions))


  # B4) other model parameters
  # Define idepth and qdepth
  idepth = -1000*out_csv_soil_discretization[input_Data$param_b90$ilayer,"Lower_m"] # in mm (positive depth)
  if (input_Data$param_b90$qlayer==0) {
    qdepth = 0
  } else {
    qdepth = -1000*out_csv_soil_discretization[input_Data$param_b90$qlayer,"Lower_m"] # in mm (positive depth)
  }
  # Save all parameters
  out_param <- with(input_Data$param_b90,
                        c("### Meteorologic site parameters -------" = NA,
                          "LAT_DEG"=coords_y, "ESLOPE_DEG"=eslope,"ASPECT_DEG"=aspect,
                          "ALB"=alb,          "ALBSN"=albsn,
                          "C1"=c1,            "C2"=c2,          "C3"=c3,
                          "WNDRAT"=wndrat,    "FETCH"=fetch,    "Z0W"=z0w,   "ZW"=zw,
                          "### Canopy parameters -------" = NA,
                          "LWIDTH"=lwidth,    "Z0G" = obsheight * czs,  "Z0S"=z0s,
                          "LPC"=lpc,          "CS"=cs,                  "CZS"=czs,
                          "CZR"=czr,          "HS"=hs,                  "HR"=hr,
                          "ZMINH"=zminh,      "RHOTP"=rhotp,            "NN"=nn,
                          # Note: u_aux_PSIM_init is not defined here as it is a vector quantity
                          "### Interception parameters -------" = NA,
                          "FRINTLAI"=frintlai, "FSINTLAI"=fsintlai,
                          "FRINTSAI"=frintsai, "FSINTSAI"=fsintsai,
                          "CINTRL"=cintrl,     "CINTRS"=cintrs,
                          "CINTSL"=cintsl,     "CINTSS"=cintss,
                          "RSTEMP"=rstemp,
                          "### Snowpack parameters -------" = NA,
                          "MELFAC"=melfac,     "CCFAC"=ccfac,
                          "LAIMLT"=laimlt,     "SAIMLT"=saimlt,
                          "GRDMLT"=grdmlt,     "MAXLQF"=maxlqf,
                          "KSNVP"=ksnvp,       "SNODEN"=snoden,
                          "### Leaf evaporation parameters (affecting PE) -------" = NA,
                          "GLMAX"=glmax,  "GLMIN"=glmin,  "CR"=radex,     "RM"=rm,
                          "R5"=r5,        "CVPD"=cvpd,    "TL"=tl,        "T1"=t1,
                          "T2"=t2,        "TH"=th,
                          "### Plant parameters (affecting soil-water supply) -------" = NA,
                          "MXKPL"=mxkpl,       "MXRTLN"=maxrlen,
                          "INITRLEN"=initrlen, "INITRDEP"=initrdep,
                          "RGRORATE"=rgrorate, "RGROPER"=rgroper,   "FXYLEM"=fxylem,
                          "PSICR"=psicr,       "RTRAD"=rrad,        "NOOUTF"=nooutf,
                          "### Soil parameters -------" = NA,
                          "IDEPTH"=idepth, "QDEPTH"=qdepth,
                          "RSSA"=rssa,     "RSSB"=rssb,     "INFEXP"=infexp,
                          "BYPAR"=bypar,   "QFPAR"=qfpar,   "QFFC"=qffc,
                          "IMPERV"=imperv, "DSLOPE"=dslope, "LENGTH_SLOPE"=slopelen,
                          "DRAIN"=drain,   "GSC"=gsc,       "GSP"=gsp,
                          "### Numerical solver parameters -------" = NA,
                          "DTIMAX"=dtimax, "DSWMAX"=dswmax, "DPSIMAX"=dpsimax))
  out_csv_param <- data.frame(param_id = names(out_param),
                              x        = unname(out_param))

  # Write out CSVs:
  input_folder   <- file.path(Julia_target_dir,paste0(Julia_prefix,"-input"))
  results_folder <- file.path(Julia_target_dir,paste0(Julia_prefix,"-LWFBrook90R_result"))
  dir.create(input_folder, showWarnings = TRUE, recursive = T)

  withr::with_options(c(scipen=100), { # temporarily switches off scientific notation
    require(dplyr)
    out_csv_soil_discretization

    out_csv_soil_discretization%>% mutate(across(where(is.numeric), round, 5)) %>% write.csv(file.path(input_folder, paste0(Julia_prefix, "_soil_discretization.csv")), row.names=FALSE, quote=FALSE)
    out_csv_soil_horizons      %>% mutate(across(where(is.numeric), round, 5)) %>% write.csv(file.path(input_folder, paste0(Julia_prefix, "_soil_horizons.csv")),       row.names=FALSE, quote=FALSE)
    out_csv_param              %>% mutate(across(where(is.numeric), round, 5)) %>% write.csv(file.path(input_folder, paste0(Julia_prefix, "_param.csv")),               row.names=FALSE, quote=FALSE)
    out_csv_initial_conditions %>% mutate(across(where(is.numeric), round, 5)) %>% write.csv(file.path(input_folder, paste0(Julia_prefix, "_initial_conditions.csv")),  row.names=FALSE, quote=FALSE)
    out_csv_meteoveg           %>% mutate(across(where(is.numeric), round, 5)) %>% write.csv(file.path(input_folder, paste0(Julia_prefix, "_meteoveg.csv")),            row.names=FALSE, quote=FALSE)
    out_csv_pdur               %>%                                                 write.csv(file.path(input_folder, paste0(Julia_prefix, "_pdur.csv")),                row.names=FALSE, quote=FALSE)
    # out_csv_precdat          %>%                                                 write.csv(file.path(input_folder, paste0(Julia_prefix, "_precdat.csv")),       row.names=FALSE, quote=FALSE)
  })

  # If run is requested also run LWFBrook90R and save results
  if(run_LWFBrook90R) {
    dir.create(results_folder, showWarnings = TRUE, recursive = T)
    start <- Sys.time()
    res <- LWFBrook90R::run_LWFB90(..., run = TRUE)
    stop <- Sys.time()

    # Write out result files and benchmark
    for(result_name in grep(".ASC", names(res),value = TRUE)) {
      # write.csv(res[[result_name]],
      #           file = file.path(results_folder, paste0(Julia_prefix, "_OUTPUT-", result_name, ".csv")))
      round_digits <- 3
      res[[result_name]] %>%
        mutate(across(where(is.numeric), round, round_digits)) %>%
        # mutate(across(!yr & !mo & !da & !doy, round, round_digits)) %>%
        write.csv(file = file.path(results_folder,
                                   paste0(Julia_prefix, "_OUTPUT-", result_name, ".csv")))
    }
    writeLines(text = c(capture.output(Sys.time()),capture.output(stop-start)),
               file.path(results_folder, paste0(Julia_prefix, "_benchmark.txt")))
  }
}
#######

# ####### Test function
# # BEA
# # 1) run first "LWFBrook90R-LWFSite-BEA-testdata.R" to prepare arguments
# # 2) then run following:
# generate_LWFBrook90Julia_Input(Julia_target_dir = ".", Julia_prefix = paste0("BEA2016-reset-", Reset),
#                                run_LWFBrook90R = TRUE,
#
#                                options_b90 = options,
#                                param_b90 = param,
#                                climate = meteo,
#                                soil = soilPuh,
#                                output= output,
#                                rtrn.output = TRUE, verbose = FALSE)
#
# # ALV8101
# # 1) run first "LWFBrook90R-TestScript-for-LWFSites_based_on_ALV8101_sen2.R" to prepare arguments
# # 2) then run following:
# generate_LWFBrook90Julia_Input(Julia_target_dir = ".", Julia_prefix = paste0("ALV8101_sen2-reset-", Reset),
#                                run_LWFBrook90R = TRUE,
#
#                     options_b90 = options,
#                     param_b90 = param,
#                     climate = meteo,
#                     soil = soilPuh,
#                     output= output,
#                     rtrn.output = TRUE, verbose = FALSE)
#
# # KAUFENRING
# # 1) run first "Schmidt-Walter-2020-Agric_For_Meteorol_results_fig_1.R" to prepare arguments
# # 2) then run following:
#
# generate_LWFBrook90Julia_Input(Julia_target_dir = ".", Julia_prefix = paste0("SW-2020_fig_1_evergreen-reset-", Reset),
#                                run_LWFBrook90R = TRUE,
#
#                                project.dir = "evergreen/",
#                                options_b90 = options_forest,
#                                #options_b90 = options_forest2, # FB modified
#                                param_b90 = parms_evergreen,
#                                climate = dplyr::rename(kau_meteo_d, windspeed=wind),
#                                soil = soil,
#                                chk_input = TRUE,
#                                output = output)
# #######


