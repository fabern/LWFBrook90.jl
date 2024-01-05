################# #
# Author: fabian.bernhard@wsl.ch
# Date: 2023-11-28
# Task: Run LWFBrook90R using input files generated with LWFBrook90.jl from a DiscretizedSPAC
# Description: This code defines an R function that takes as arguments five paths
#              to five files: meteo.csv, soil.csv, df_opts.csv, df_parms_a.csv, df_parms_b.csv
#              which were exported with a function from a DiscretizedSPAC.
#
#              It reads this input data set and simulates with LWFBrook90R instead
#              of LWFBrook90.jl. This is useful for validation purposes.
#
#              Note that this is provided as-is for convenience. It is likely
#              that the approach is not suited for edge cases.
#
################# #

# # source: https://pschmidtwalter.github.io/LWFBrook90R/
# library(LWFBrook90R)
# data(slb1_meteo, slb1_soil)
#
# opts_target  <- set_optionsLWFB90()
# parms_target <- set_paramLWFB90()
# soil_target  <- cbind(slb1_soil, hydpar_wessolek_tab(texture = slb1_soil$texture))
# meteo_target <- slb1_meteo
#
# # run the model and capture results
# lwfb90_res <- run_LWFB90(options_b90 = opts_target,
#                          param_b90 = parms_target,
#                          climate = meteo_target,
#                          soil = soil_target)


simulate_LWFBrook90.jl_DiscretizedSPAC_withLWFBrook90R <-
    function(meteo.csv.path, soil.csv.path, opts.csv.path, parms_a.csv.path, parms_b.csv.path){

    library(LWFBrook90R)
    library(readr)
    library(lubridate)
    library(dplyr)

    # check for version
    stopifnot(packageVersion("LWFBrook90R") == "0.5.3")

    # read input files ###########################################################
    meteo <- readr::read_csv(meteo.csv.path)
    soil  <- readr::read_csv(soil.csv.path)
    df_opts   <- readr::read_csv(opts.csv.path)
    df_parms_a  <- readr::read_csv(parms_a.csv.path)
    df_parms_b  <- readr::read_csv(parms_b.csv.path)

    opts = as.list(as.data.frame(df_opts))
    parms = append(as.list(as.data.frame(df_parms_a)),
                       # insert below at specific position
                       as.list(as.data.frame(df_parms_b)), 95)
        # which(names(parms_old) == "pdur")
        # which(names(parms) == "pdur")

        # compare_two_lists <- function(l1, l2){
        #   N = length(names(l1))
        #   comparisons = structure(.Data = rep(FALSE, N), .Names = names(l1))
        #   for (name in names(l1)) {
        #     comparisons[name] = identical(l1[name], l2[name])
        #   }
        #   if (!all(comparisons)) {print("Some are unequal:"); print(comparisons[!comparisons])}
        #   return(all(comparisons))
        # }
        #   # identical(parms_old, parms_old)
        #   # identical(parms_old, parms)
        # compare_two_lists(opts_old, opts)
        #   # typeof(opts_old$prec_interval) # integer
        #   # typeof(opts$prec_interval) # double
        # compare_two_lists(parms_old, parms)

    # Fix some errors in output
    opts$startdate <- lubridate::ymd(as.character(opts$startdate))
    opts$enddate   <- lubridate::ymd(as.character(opts$enddate))

    parms$shp_budburst <- NaN
    parms$shp_leaffall <- NaN
    parms$shp_optdoy   <- NaN

    # undebug(run_LWFB90)
    sim_out <- run_LWFB90(options_b90 = opts,
                          param_b90   = parms,
                          climate     = meteo,
                          soil        = soil,
                          rtrn_input = TRUE)
    stopifnot(sim_out$error_code == 0)
        # parms$rootden_table # TODO: alternatively use inputs for r_lwfbrook90 (see below, and above sim_out$model_input)
        # if (FALSE) {
        #     # soil_layers = read_table_skip_subheader(file = soil_layers_file, sep=",")  %>%
        #     #   rename(upper = Upper_m,
        #     #          lower = Lower_m) %>%
        #     #   mutate(thick = -(lower - upper)*1000) %>% # thickness in mm
        #     #   mutate(midpoint = upper - thick/2/1000) %>%
        #     #   mutate(#mat = 1:n(),  TODO: mat is not correct
        #     #     mat = NA,
        #     #     layer = 1:n()) %>%
        #     #   select(layer, upper, lower, midpoint, thick,
        #     #          mat,
        #     #          psiini=uAux_PSIM_init_kPa,
        #     #          rootden = Rootden_)
        #     # # Fix meteo by using absolute values of densef, height, lai, and add age
        #     # meteo2 <- meteo %>%
        #     #   mutate(age = param$AGE_baseline_yrs + as.numeric(dates - dates[1],"days")/365,
        #     #          densef=densef/100*param$DENSEF_baseline_,
        #     #          height=height/100*param$HEIGHT_baseline_m,
        #     #          lai=lai/100*param$MAXLAI,
        #     #          sai=sai/100*param$SAI_baseline_)
        #     #
        #     # # prepare simulation input
        #     # # run simulation #############################################################
        #     # sim_out <- r_lwfbrook90(siteparam = sim_in$siteparam,
        #     #                        climveg = sim_in$climveg,
        #     #                        precdat = sim_in$precdat,
        #     #                        param = sim_in$param,
        #     #                        pdur = sim_in$pdur,
        #     #                        soil_materials = sim_in$soil_materials,
        #     #                        soil_nodes = sim_in$soil_nodes,
        #     #                        output_log = verbose#,
        #     #                        # timelimit = Inf
        #     # )
        #
        #     a_siteparam <- data.frame(lubridate::year(opts$startdate), #as.integer(format(min(meteo$dates),"%Y")),
        #                              lubridate::year(opts$enddate),
        #                              param$coords_y, # param$LAT_DEG,
        #                              parms$snowini,  # IC$u_SNOW_init_mm,
        #                              parms$gwatini,  # IC$u_GWAT_init_mm,
        #                              1.) # hardcoded: options.b90$prec.interval
        #     a_climveg <- meteo %>% rename(vappres = vappress) %>%
        #       mutate(densef = NA, height = NA, lai = NA, sai = NA, age = NA) %>%
        #       mutate(yr = lubridate::year(dates),#as.integer(format(dates,"%Y")),
        #              mo = lubridate::month(dates),#as.integer(format(dates,"%m")),
        #              da = lubridate::day(dates),#as.integer(format(dates,"%d")),
        #              mesfl = 0.) %>% # hardcoded: MESFL = 0
        #       select(yr, mo, da, globrad, tmax, tmin, vappres, windspeed, prec, mesfl,
        #              densef, height, lai, sai, age)
        #     a_param <- param_to_rlwfbrook90(parms, opts$imodel)
        #     a_pdur = parms$pdur
        #     parms[c("soil_nodes", "soil_materials")] <- soil_to_param(soil, opts$imodel)
        #     # parms$soil_materials <- parms$soil_materials[,  c("mat", "ths", "thr", "alpha", "npar", "ksat", "tort", "gravel")]
        #     a_soil_materials <- parms$soil_materials
        #     a_soil_nodes <- parms$soil_nodes
        #     sim_in <- list(siteparam     = a_siteparam,
        #                     climveg        = a_climveg,
        #                     precdat        = NULL,
        #                     param          = a_param,
        #                     pdur           = a_pdur,
        #                     soil_materials = a_soil_materials,
        #                     soil_nodes     = a_soil_nodes)
        #
        #       # postprocess simulation #####################################################
        #       # # note: copied from LWFBrook90R (github.com/pschmidtwalter/LWFBrook90R)
        #       # TODO... see fct-run_LWFBrook90R_with_LWFBrook90jl_Input.R
        #
        #     # In and out of r_lwfbrook90:
        #     list(simout=sim_out, simin=sim_in)
        # }

    return(sim_out)
}

postprocess_simulation <- function(sim_out, depths_to_read_out_m = c()){

    library(ggplot2)
    library(tidyr)
    library(dplyr)

    ##### Postprocess simulation ##############
    # sim_out$output
    # sim_out$layer_output

    # treat aboveground output
    dat_aboveground <- select(sim_out$output, yr, mo, da, doy, rfal, sfal, #prec, # rfal+sfal instead of prec # rthr,
                              flow, evap, seep, snow, swat, gwat, intr, ints, sint, isvp) %>%
      select(-doy) %>%
      tidyr::unite(col = "date", yr, mo, da) %>%
      mutate(date = lubridate::ymd(date)) %>%
      tidyr::pivot_longer(cols = -date,
                          # keep facets in ggplot in the same order as in data.frame
                          names_transform = list(name = forcats::fct_inorder))
    # treat layer output
    dat_theta <- select(sim_out$layer_output, yr, mo, da, nl, theta) %>%
      tidyr::unite(col = "date", yr, mo, da) %>%
      mutate(date = lubridate::ymd(date)) %>%
      left_join({sim_out$model_input$param_b90$soil_nodes %>%
          mutate(lower = -0 + -cumsum(thick)) %>%
          mutate(upper = c(0, lower[-n()])) %>%
          select(nl = layer, midpoint, upper, lower)},
          by = "nl")

    dat_psi <- select(sim_out$layer_output, yr, mo, da, nl, psimi) %>%
      tidyr::unite(col = "date", yr, mo, da) %>%
      mutate(date = lubridate::ymd(date)) %>%
      left_join({sim_out$model_input$param_b90$soil_nodes %>%
          mutate(lower = -0 + -cumsum(thick)) %>%
          mutate(upper = c(0, lower[-n()])) %>%
          select(nl = layer, midpoint, upper, lower)},
          by = "nl")

    if (length(depths_to_read_out_m) == 0) {
      return(list(theta = dat_theta,
                  psi   = dat_psi))
    }

    # Define log modulus transformation for plotting (for negative (and positive) log axis)
    log_modulus_trans <- function(base = exp(1)) {
      scales::trans_new(name = "log_modulus",
                        transform = function(x) sign(x) * log(abs(x) + 1, base = base),
                        inverse   = function(x) sign(x) * ( base^(abs(x)) - 1 ))}
    log_modulus_rev_trans <- function(base = exp(1)) {
      scales::trans_new(name = "log_modulus_rev",
                        transform = function(x) -sign(x) * log(abs(x) + 1, base = base),
                        inverse   = function(x) -sign(x) * ( base^(abs(x)) - 1 ))}
    scale_y_logModulus <- function(...) {
      scale_y_continuous(trans = "log_modulus",
                         labels=function(n) format(n, scientific=FALSE),
                         breaks = c((+1)*(10^(0:10) %*% t(c(1, 3)) ),
                                    0,
                                    (-1)*(10^(0:10) %*% t(c(1, 3)) )),
                         minor_breaks = c((+1)*(10^(-0:10) %*% t(c(1:9)) ),
                                          (-1)*(10^(-0:10) %*% t(c(1:9)) )),
                         ...)}

    # # Plotting
    # # a) plotting aboveground
    # pl0_aboveground <- ggplot(dat_aboveground,
    #                           aes(x=date, y=value, color = name)) +
    #   geom_line() + theme_bw() + facet_grid(name~., scales = "free_y") +
    #   theme(axis.text.y = element_text(size=5), legend.key.height = unit(1,"line")) +
    #   labs(x=NULL)
    #

    # # b) plotting belowground
    # pl1_belowground_theta <-ggplot(dat_theta, aes(x=date, y=midpoint, fill = theta)) +
    #   geom_rect(aes(xmin=date, xmax=date+1, ymin=lower, ymax=upper)) +
    #   theme_bw() + labs(fill = "θ (-)", y = "Depth (m)", x=NULL)
    # pl1_belowground_psimi <- ggplot(dat_psi, aes(x=date, y=midpoint, fill = psimi)) +
    #   geom_rect(aes(xmin=date, xmax=date+1, ymin=lower, ymax=upper)) +
    #   theme_bw() + labs(fill = "ψ_m (kPa)", y = "Depth (m)", x=NULL) +
    #   # define color palette separating saturated from unsaturated and using log_modulus transform
    #   scale_fill_gradientn(
    #     # sharp break between blue (unsaturated) and red (saturated) at 0:
    #     # https://github.com/tidyverse/ggplot2/issues/3738#issuecomment-1336166757
    #     # colours = c("white", scales::muted("blue"), scales::muted("red"), "black"), values = c(0, 0.5, 0.5001, 1.0),
    #     colours = c("white", scales::muted("blue"), "red", "black"), values = c(0, 0.4999, 0.5001, 1.0), # This highlights saturation = 1 in read already
    #     rescaler = ~ scales::rescale_mid(.x, mid = 0),
    #     trans ="log_modulus",
    #     labels=function(n) format(n, scientific=FALSE),
    #     breaks = c((+1)*(10^(0:10) %*% t(c(1, 3)) ),
    #                0,
    #                (-1)*(10^(0:10) %*% t(c(1, 3)) )),
    #     minor_breaks = c((+1)*(10^(-0:10) %*% t(c(1:9)) ),
    #                      (-1)*(10^(-0:10) %*% t(c(1:9)) )))
    #     # Testing color scale:
    #     # ggplot(data.frame(x = rep(1:5, each = 5),
    #     #                   y = rep(-6:3, 10)) %>%
    #     #          mutate(z = y + sign(y)*x),
    #     #        aes(x=x, y=y, fill=z)) +
    #     #   geom_rect(aes(xmin=x, xmax=x+1, ymin=y, ymax=y+1)) +
    #     #   scale_fill_gradientn(
    #     #     # sharp break between blue (unsaturated) and red (saturated) at 0:
    #     #     # https://github.com/tidyverse/ggplot2/issues/3738#issuecomment-1336166757
    #     #     colours = c("white", scales::muted("blue"), scales::muted("red"), "black"), values = c(0, 0.5, 0.5001, 1.0),
    #     #     rescaler = ~ scales::rescale_mid(.x, mid = 0))

    # specific depths
    # depths_to_read_out_m <- -c(150,500,800)/1000 # units: (mm) => (m)
    idx_to_read_out <- c()
    for (it in seq_along(depths_to_read_out_m)) {
      idx_to_read_out[it] <- which.min(abs(depths_to_read_out_m[it] - sim_out$model_input$param_b90$soil_nodes$midpoint))
    }
    depths_actually_read_out <-
      sim_out$model_input$param_b90$soil_nodes$midpoint[idx_to_read_out]

    # pl2_belowground_theta <- ggplot(filter(dat_theta, nl %in% c(idx_to_read_out)),
    #                                 aes(x=date, y=theta, color = as.factor(midpoint))) +
    #   geom_line() + theme_bw() +
    #   theme_bw() + labs(y = "θ (-)", color = "Depth (m)", x=NULL) + theme(legend.position = c(0.87, 0.25))
    # pl2_belowground_psi <-ggplot(filter(dat_psi, nl %in% c(idx_to_read_out)),
    #                              aes(x=date, y=psimi, color = as.factor(midpoint))) +
    #   geom_line() + theme_bw() +
    #   theme_bw() + labs(y = "ψ_m (kPa)", color = "Depth (m)", x=NULL) + theme(legend.position = c(0.87, 0.25)) +
    #   # use log modulus scale
    #   scale_y_logModulus()
    # gridExtra::grid.arrange(pl2_belowground_theta, pl2_belowground_psi)

    if (length(depths_to_read_out_m) > 0) {
      return(list(theta = filter(dat_theta, nl %in% c(idx_to_read_out)),
                  psi   = filter(dat_psi,   nl %in% c(idx_to_read_out))))
    }


    # # Plot Mass Balance Error
    # BALERR_computations <- sim_out$output %>%
    #   tibble() %>%
    #   # consider soil domain as control volume
    #   mutate(BALERR_SWATI_cumInflow  = cumsum(slfl - byfl),
    #          BALERR_SWATI_cumOutflow = cumsum(vrfln + dsfl + tran + slvp)) %>%
    #   mutate(BALERR_SWATI_cumInflow  = BALERR_SWATI_cumInflow  - BALERR_SWATI_cumInflow[1],
    #          BALERR_SWATI_cumOutflow = BALERR_SWATI_cumOutflow - BALERR_SWATI_cumOutflow[1]) %>%
    #   # consider entire model domain as control volume
    #   # TODO
    #   select(yr, mo, da, doy, swat, BALERR_SWATI_cumInflow, BALERR_SWATI_cumOutflow) %>%
    #   # compute error
    #   mutate(error_mm =
    #            (BALERR_SWATI_cumInflow - BALERR_SWATI_cumOutflow) -
    #            # Note: at position [1] cumInflow/cumOutflow brings us from swat[1] to swat[2] -> hence the lag:
    #            #       However the errors seems smaller without the lag.
    #            # lag((BALERR_SWATI_cumInflow - BALERR_SWATI_cumOutflow), default=0) -
    #            (swat - swat[1]))
    # pl_error <- BALERR_computations %>%
    #   tidyr::pivot_longer(cols = c(BALERR_SWATI_cumInflow, BALERR_SWATI_cumOutflow, error_mm, swat)) %>%
    #   ggplot(data = ., aes(x=doy, color = name)) +
    #   geom_line(aes(y=value)) +
    #   facet_grid((name=="error_mm") ~ ., scales = "free_y") +
    #   theme_bw() + labs(y="Cumulative Fluxes or Storage [mm]", x = "Day of Year", color = NULL) +
    #   theme(strip.background = element_blank(),
    #         strip.text = element_blank(), legend.background = element_blank(),
    #         legend.position = c(0.98,0.51), legend.justification = c(1,0))
    # ggsave(filename = gsub(".png","error.png",fnamePNG),
    #        width = 1200, height = 1000, units = "px", dpi = 220, pl_error)
  }


meteo.csv.path   = "../../generated_LWFBrook90R_input_files/meteo.csv"
soil.csv.path    = "../../generated_LWFBrook90R_input_files/soil.csv"
opts.csv.path    = "../../generated_LWFBrook90R_input_files/df_opts.csv"
parms_a.csv.path = "../../generated_LWFBrook90R_input_files/df_parms_a.csv"
parms_b.csv.path = "../../generated_LWFBrook90R_input_files/df_parms_b.csv"
sim_out <- simulate_LWFBrook90.jl_DiscretizedSPAC_withLWFBrook90R(meteo.csv.path, soil.csv.path, opts.csv.path, parms_a.csv.path, parms_b.csv.path)

dat_postproc <- postprocess_simulation(sim_out)
dat_postproc <- postprocess_simulation(sim_out, depths_to_read_out_m = -c(150,500,800)/1000) # units: (mm) => (m)
dat_postproc$theta
dat_postproc$psi
