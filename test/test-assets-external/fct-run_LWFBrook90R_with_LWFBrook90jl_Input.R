################# #
# Author: fabian.bernhard@wsl.ch
# Date: 2022-03-15
# Task: Run LWFBrook90R using input files for LWFBrook90.jl
# Description: This code defines an R function that takes as arguments a folder-
#              path and a prefix, which define an input data set for LWFBrook90.jl.
#              
#              It reads this input data set and simulates with LWFBrook90R instead
#              of LWFBrook90.jl. This is useful for validation purposes.
#
#              Note that this is provided as-is for convenience. It is likely 
#              that the approach is not suited for edge cases.
#
################# #

simulate_LWFBrook90.jl_input_folder_withLWFBrook90R <-
  function(folder="input_LWFBrook90R",prefix="BEA2016-reset-FALSE",
           refine_discr_factor = 1,
           verbose = TRUE){
    # refine_discr_factor: a factor that divides each discretized cell into that many subcells
    # folder="input_LWFBrook90R"
    # prefix="BEA2016-reset-FALSE"
    # folder = "../test-assets/DAV-2020/input-files"
    # prefix = "DAV_LW1_def"
    # refine_discr_factor = 1
    # verbose = TRUE

    require(LWFBrook90R)
    require(data.table)
    require(dplyr)
    
    # read input files ###########################################################
    meteo_file              <- file.path(folder, paste0(prefix, "_meteoveg.csv"))
    soil_material_file      <- file.path(folder, paste0(prefix, "_soil_horizons.csv"))
    soil_layers_file        <- file.path(folder, paste0(prefix, "_soil_discretization.csv"))
    initial_conditions_file <- file.path(folder, paste0(prefix, "_initial_conditions.csv"))
    param_file              <- file.path(folder, paste0(prefix, "_param.csv"))
    storm_durations_file    <- file.path(folder, paste0(prefix, "_meteo_storm_durations.csv"))
    
    # Define function to read in csv files with a secondary subheader (ignored)
    read_table_skip_subheader <- function(file, ...) {
      header <- read.table(file, nrows = 1, header = TRUE, stringsAsFactors = FALSE, ...)
      dat    <- read.table(file, skip = 2,  header = FALSE, ...)
      colnames( dat ) <- colnames(header)
      return(dat)
    }
    
    # Read in 1) meteo
    meteo <- read_table_skip_subheader(
      file = meteo_file, 
      sep=",") %>%
      rename(dates=dates, 
             globrad=globrad_MJDayM2, 
             tmax=tmax_degC, 
             tmin=tmin_degC, 
             vappres=vappres_kPa, 
             windspeed=windspeed_ms, 
             prec=prec_mmDay, 
             densef=densef_percent, 
             height=height_percent, 
             lai=lai_percent, 
             sai=sai_percent) %>%
      mutate(dates = as.Date(dates)) %>% tibble()
    
    # Read in 2) soil
    soil_materials <- read_table_skip_subheader(file = soil_material_file, sep=",") %>%
      rename(ths=ths_volFrac,
             thr=thr_volFrac,
             alpha=alpha_perMeter,
             npar = npar_,
             ksat = ksat_mmDay,
             tort = tort_,
             gravel = gravel_volFrac,
             upper = Upper_m, 
             lower = Lower_m)
    
    soil_layers = read_table_skip_subheader(file = soil_layers_file, sep=",")  %>%
      rename(upper = Upper_m, 
             lower = Lower_m) %>%
      mutate(thick = -(lower - upper)*1000) %>% # thickness in mm
      mutate(midpoint = upper - thick/2/1000) %>%
      mutate(#mat = 1:n(),  TODO: mat is not correct
             mat = NA,
             layer = 1:n()) %>%
      select(layer, upper, lower, midpoint, thick, 
             mat, 
             psiini=uAux_PSIM_init_kPa, 
             rootden = Rootden_)
    
    # define which material is in which layer
    for (mat in seq_len(nrow(soil_materials))) {
      layers_in_current_material <- soil_layers$lower >= soil_materials[mat,]$lower & soil_layers$upper <= soil_materials[mat,]$upper  
      if (any(layers_in_current_material)){
        soil_layers[layers_in_current_material,]$mat <- mat
      }
    }
    
    # 2a) Increase resolution of soil discretization
    soil_layers_increased <- soil_layers %>%
      arrange(desc(upper)) %>%
      # multiply the layers:
      slice(rep(1:n(), each = refine_discr_factor)) %>%
      # Update upper, lower, thickness
      group_by(layer) %>%
      mutate(new_thick = thick/refine_discr_factor) %>%
      mutate(new_upper    = upper - new_thick/1000 * ((row_number())-1),
             new_midpoint = upper - new_thick/1000 * ((row_number())-1) - new_thick/1000,
             new_lower    = upper - new_thick/1000 * ((row_number())) ) %>%
      # cleanup
      ungroup() %>%
      mutate(layer = 1:n()) %>%
      select(layer,
             upper = new_upper, lower = new_lower, midpoint = new_midpoint, thick = new_thick,
             mat, psiini, rootden)
    
    
    soil_layers_increased <- soil_layers_increased %>%
      # select(-upper, -lower)
      select(layer, midpoint, thick, 
             mat, psiini, rootden)
    
    # Read in 3) site properties affecting model parameters
    IC <- read.table(initial_conditions_file, header = TRUE, sep = ",") %>%
      # get rid of Isotope information:
      select(param_id, amount) %>%
      # make input structure as named vector
      {structure(.$amount, .Names = .$param_id)} %>%
      as.list()
    
    param <- read.table(param_file, header = TRUE, comment.char = "#", sep = ",") %>%
      {structure(.$x, .Names = .$param_id)} %>%
      as.list()

    storm_durat <- read.table(storm_durations_file, header = TRUE, sep = ",")
    
    
    # 4) Parse the parameters back into LWFBrook90R format
    # test_params <- set_paramLWFB90()
    # Parse the read in back into LWFBrook90R format
    test_params <- list(
      eslope = param$ESLOPE_DEG,
      aspect = param$ASPECT_DEG,
      alb = param$ALB,
      albsn = param$ALBSN,
      c1 = param$C1,
      c2 = param$C2,
      c3 = param$C3,
      wndrat = param$WNDRAT,
      fetch = param$FETCH,
      z0w = param$Z0W,
      zw = param$ZW,
      lwidth = param$LWIDTH,
      z0s = param$Z0S,
      lpc = param$LPC,
      cs =  param$CS,
      czs = param$CZS,
      czr = param$CZR,
      hs =  param$HS,
      hr =  param$HR,
      zminh = param$ZMINH,
      rhotp = param$RHOTP,
      nn =  param$NN,
      rstemp = param$RSTEMP,
      frintlai = param$FRINTLAI,
      fsintlai = param$FSINTLAI,
      frintsai = param$FRINTSAI,
      fsintsai = param$FSINTSAI,
      cintrl = param$CINTRL,
      cintrs = param$CINTRS,
      cintsl = param$CINTSL,
      cintss = param$CINTSS,
      melfac = param$MELFAC,
      ccfac = param$CCFAC,
      laimlt = param$LAIMLT,
      saimlt = param$SAIMLT,
      grdmlt = param$GRDMLT,
      maxlqf = param$MAXLQF,
      ksnvp = param$KSNVP,
      snoden = param$SNODEN,
      glmax = param$GLMAX,
      glmin = param$GLMIN,
      rm = param$RM,
      r5 = param$R5,
      cvpd = param$CVPD,
      tl = param$TL,
      t1 = param$T1,
      t2 = param$T2,
      th = param$TH,
      mxkpl = param$MXKPL,
      maxrlen = param$MXRTLN,
      initrlen = param$INITRLEN,
      initrdep = param$INITRDEP,
      rgrorate = param$RGRORATE,
      rgroper = param$RGROPER,
      fxylem = param$FXYLEM,
      psicr = param$PSICR,
      rrad = param$RTRAD,
      nooutf = param$NOOUTF,
      rssa = param$RSSA,
      rssb = param$RSSB,
      infexp = param$INFEXP,
      bypar = param$BYPAR,
      qfpar = param$QFPAR,
      qffc = param$QFFC,
      imperv = param$IMPERV,
      dslope = param$DSLOPE,
      slopelen = param$LENGTH_SLOPE,
      drain = param$DRAIN,
      gsc = param$GSC,
      gsp = param$GSP,
      dtimax = param$DTIMAX,
      dswmax = param$DSWMAX,
      dpsimax = param$DPSIMAX)
    # Parse further parameters
    test_params$obsheight = param$Z0G / test_params$czs
    test_params$radex <- param$CR
    
    test_params$ilayer <- first(which(soil_materials$lower <= -param$IDEPTH_m))
    test_params$qlayer <- ifelse(param$QDEPTH_m == 0, 
                                 0,
                                 first(which(soil_materials$lower <= -param$QDEPTH_m)))
    
    test_params$intrainini <- IC$u_INTR_init_mm
    test_params$intsnowini <- IC$u_INTS_init_mm
    
    test_params$ndays <- as.numeric(diff(range(meteo$dates))) + 1
    
    # No need to have the actual materials and layers, just the nrow(.) should work
    test_params$soil_materials <- data.frame(rep(1,nrow(soil_materials)))
    test_params$soil_nodes <- data.frame(rep(1,nrow(soil_layers_increased)))
    test_params$intrainini
    # Fix meteo by using absolute values of densef, height, lai, and add age
    meteo2 <- meteo %>%
      mutate(age = param$AGE_baseline_yrs + as.numeric(dates - dates[1],"days")/365,
             densef=densef/100*param$DENSEF_baseline_,
             height=height/100*param$HEIGHT_baseline_m,
             lai=lai/100*param$MAXLAI,
             sai=sai/100*param$SAI_baseline_)

    # prepare simulation input
    sim_in <- list(siteparam = data.frame(as.integer(format(min(meteo2$dates),"%Y")),
                                         as.integer(format(min(meteo2$dates),"%j")),
                                         param$LAT_DEG, 
                                         IC$u_SNOW_init_mm,
                                         IC$u_GWAT_init_mm,
                                         1.), # hardcoded: options.b90$prec.interval
                  climveg = 
                    meteo2 %>%
                    mutate(yr = as.integer(format(dates,"%Y")),
                           mo = as.integer(format(dates,"%m")),
                           da = as.integer(format(dates,"%d")),
                           mesfl = 0.) %>% # hardcoded: MESFL = 0
                    select(yr, mo, da, globrad, tmax, tmin, vappres, windspeed, prec, mesfl,
                           densef, height, lai, sai, age),
                  precdat = NULL,
                  param = param_to_rlwfbrook90(test_params, "MvG"), #hardcoded: options_b90$imodel as "MvG"
                  pdur = as.numeric(storm_durat$average_storm_duration_h),
                  soil_materials = rename(soil_materials, layer = HorizonNr) %>% select(-upper, -lower),
                  soil_nodes = soil_layers_increased)

    # run simulation #############################################################
    sim_out <- r_lwfbrook90(siteparam = sim_in$siteparam,
                           climveg = sim_in$climveg,
                           precdat = sim_in$precdat,
                           param = sim_in$param,
                           pdur = sim_in$pdur,
                           soil_materials = sim_in$soil_materials,
                           soil_nodes = sim_in$soil_nodes,
                           output_log = verbose#,
                           # timelimit = Inf
    )
    
    # postprocess simulation #####################################################
    # daily outputs
    sim_out$daily_output <- data.table::data.table(sim_out$daily_output)
    data.table::setnames(sim_out$daily_output, 
                         names(sim_out$daily_output),
                         c('yr','mo','da','doy','rfal','rint','sfal','sint','rthr','rsno',
                           'rnet','smlt','snow','swat','gwat','intr', 'ints','evap','tran','irvp',
                           'isvp','slvp','snvp','pint','ptran','pslvp','flow','seep',
                           'srfl','slfl','byfl','dsfl','gwfl','vrfln','safrac',
                           'stres','adef','awat','relawat','nits','balerr', 'slrad',
                           'solnet', 'lngnet', 'aa', 'asubs'))
    # layer outputs
    sim_out$layer_output <- data.table::rbindlist(lapply(seq(dim(sim_out$layer_output)[3]),
                                                        function(x) data.frame(sim_out$layer_output[ , , x])),
                                                 idcol = "nl")
    data.table::setnames(sim_out$layer_output, paste0("X", 1:16),
                         c('yr','mo','da','doy','swati','theta','wetnes','psimi','psiti','infl',
                           'byfl','tran','slvp','vrfl','dsfl','ntfl'))
    
    sim_out$layer_output <- sim_out$layer_output[order(sim_out$layer_output$yr,
                                                     sim_out$layer_output$doy,
                                                     sim_out$layer_output$nl),]
    
    return(list(simout=sim_out, simin=sim_in))
  }