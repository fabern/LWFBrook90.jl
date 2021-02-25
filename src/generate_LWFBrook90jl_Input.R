#################
# Author: fabian.bernhard@wsl.ch
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
  mesfl_daily <- data.frame(dates = veg_daily$dates, mesfl = 0) # TODO(bernhard): impelement mesfl
  clim_daily  <- original_arguments$climate %>%
    select(dates, globrad, tmax, tmin, vappres, windspeed, prec) %>%
    filter(dates >= input_Data$options_b90$startdate) %>%
    filter(dates <= input_Data$options_b90$enddate)

  out_csv_climveg <- clim_daily %>%
    dplyr::left_join(mesfl_daily, by="dates") %>%
    dplyr::left_join(veg_daily, by="dates")

  # B) time-independent parameters:
  # B1) climate data
  out_csv_pdur <- input_Data$param_b90$pdur
  # out_csv_precdat <- NULL # TODO(bernhard): currently not implemented

  # B2) soil parameters
  out_csv_soil_nodes     <- input_Data$param_b90$soil_nodes %>%
    select(layer, midpoint, thick, mat, psiini, rootden)
  out_csv_soil_materials <- input_Data$param_b90$soil_materials

  # B3) site parameters
  simyears <- seq(from = as.integer(format(input_Data$options_b90$startdate,"%Y")),
                  to = as.integer(format(input_Data$options_b90$enddate,"%Y")),
                  by = 1)
  out_csv_siteparam <- data.frame(simyears[1],
                                  as.integer(format(input_Data$options_b90$startdate, "%j")),
                                  input_Data$param_b90$coords_y,
                                  input_Data$param_b90$snowini,
                                  input_Data$param_b90$gwatini,
                                  input_Data$options_b90$prec_interval)

  # B4) other model parameters
  # Variant 1: (unfinished)
  # NOTE: from input_Data$param_b90 remove:
  #       1. soil_nodes
  #       2. soil_materials
  # out_csv_param <- input_Data$param_b90
  # names(out_csv_param) %in% c("soil_nodes", "soil_materials")

  # Variant 2: (ugly but finished and compatible with current Julia code)
  out_param <- with(input_Data$param_b90,
                        c("ndays"=ndays, "0_heat"=0, "eslope"=eslope, "aspect"=aspect,
                          "alb"=alb, "albsn"=albsn, "c1"=c1, "c2"=c2, "c3"=c3,
                          "wndrat"=wndrat, "fetch"=fetch, "z0w"=z0w,
                          "zw"=zw, "lwidth"=lwidth, "obsheight_x_czs" = obsheight * czs, "z0s"=z0s, "lpc"=lpc, "cs"=cs,
                          "czs"=czs, "czr"=czr, "hs"=hs, "hr"=hr,
                          "zminh"=zminh, "rhotp"=rhotp, "nn"=nn, "rstemp"=rstemp,
                          "intrainini"=intrainini, "intsnowini"=intsnowini, "frintlai"=frintlai,
                          "fsintlai"=fsintlai, "frintsai"=frintsai, "fsintsai"=fsintsai,
                          "cintrl"=cintrl, "cintrs"=cintrs, "cintsl"=cintsl,
                          "cintss"=cintss, "melfac"=melfac, "ccfac"=ccfac,
                          "laimlt"=laimlt, "saimlt"=saimlt, "grdmlt"=grdmlt,
                          "maxlqf"=maxlqf, "ksnvp"=ksnvp, "snoden"=snoden,
                          "glmax"=glmax, "radex"=radex, "glmin"=glmin, "rm"=rm,
                          "r5"=r5, "cvpd"=cvpd, "tl"=tl, "t1"=t1,
                          "t2"=t2, "th"=th, "mxkpl"=mxkpl, "maxrlen"=maxrlen,
                          "initrlen"=initrlen, "initrdep"=initrdep, "rgrorate"=rgrorate,
                          "rgroper"=rgroper, "fxylem"=fxylem, "psicr"=psicr,
                          "rrad"=rrad, "nooutf"=nooutf, "N_soil_nodes"=nrow(soil_nodes),
                          "N_soil_materials"=nrow(soil_materials), "ilayer"=ilayer, "qlayer"=qlayer,
                          "is_MvG_aka_iModel"=ifelse(input_Data$options_b90$imodel == "MvG", 1, 0),
                          "rssa"=rssa, "rssb"=rssb,
                          "infexp"=infexp, "bypar"=bypar, "qfpar"=qfpar, "qffc"=qffc,
                          "imperv"=imperv, "dslope"=dslope, "slopelen"=slopelen,
                          "drain"=drain, "gsc"=gsc, "gsp"=gsp, "dtimax"=dtimax,
                          "dswmax"=dswmax, "dpsimax"=dpsimax))
  out_csv_param <- data.frame(x=unname(out_param),
                              param_id=names(out_param))

  # Write out CSVs:
  input_folder   <- file.path(Julia_target_dir,paste0(Julia_prefix,"-input"))
  results_folder <- file.path(Julia_target_dir,paste0(Julia_prefix,"-LWFBrook90R_result"))
  dir.create(input_folder, showWarnings = TRUE, recursive = T)

  write.csv(out_csv_soil_nodes,     file.path(input_folder, paste0(Julia_prefix, "_soil_nodes.csv")),     row.names = FALSE)
  write.csv(out_csv_soil_materials, file.path(input_folder, paste0(Julia_prefix, "_soil_materials.csv")), row.names = FALSE)
  write.csv(out_csv_pdur,           file.path(input_folder, paste0(Julia_prefix, "_pdur.csv")),           row.names = FALSE)
  write.csv(out_csv_param,          file.path(input_folder, paste0(Julia_prefix, "_param.csv")),          row.names = FALSE)
  # write.csv(out_csv_precdat,        file.path(input_folder, paste0(Julia_prefix, "_precdat.csv")),        row.names = FALSE)
  write.csv(out_csv_climveg,        file.path(input_folder, paste0(Julia_prefix, "_climveg.csv")),        row.names = FALSE)
  write.csv(out_csv_siteparam,      file.path(input_folder, paste0(Julia_prefix, "_siteparam.csv")),      row.names = FALSE)

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

####### Test function
# BEA
# 1) run first "LWFBrook90R-LWFSite-BEA-testdata.R" to prepare arguments
# 2) then run following:
generate_LWFBrook90Julia_Input(Julia_target_dir = ".", Julia_prefix = paste0("BEA2016-reset-", Reset),
                               run_LWFBrook90R = TRUE,

                               options_b90 = options,
                               param_b90 = param,
                               climate = meteo,
                               soil = soilPuh,
                               output= output,
                               rtrn.output = TRUE, verbose = FALSE)

# ALV8101
# 1) run first "LWFBrook90R-TestScript-for-LWFSites_based_on_ALV8101_sen2.R" to prepare arguments
# 2) then run following:
generate_LWFBrook90Julia_Input(Julia_target_dir = ".", Julia_prefix = paste0("ALV8101_sen2-reset-", Reset),
                               run_LWFBrook90R = TRUE,

                    options_b90 = options,
                    param_b90 = param,
                    climate = meteo,
                    soil = soilPuh,
                    output= output,
                    rtrn.output = TRUE, verbose = FALSE)

# KAUFENRING
# 1) run first "Schmidt-Walter-2020-Agric_For_Meteorol_results_fig_1.R" to prepare arguments
# 2) then run following:

generate_LWFBrook90Julia_Input(Julia_target_dir = ".", Julia_prefix = paste0("SW-2020_fig_1_evergreen-reset-", Reset),
                               run_LWFBrook90R = TRUE,

                               project.dir = "evergreen/",
                               options_b90 = options_forest,
                               #options_b90 = options_forest2, # FB modified
                               param_b90 = parms_evergreen,
                               climate = dplyr::rename(kau_meteo_d, windspeed=wind),
                               soil = soil,
                               chk_input = TRUE,
                               output = output)
#######


