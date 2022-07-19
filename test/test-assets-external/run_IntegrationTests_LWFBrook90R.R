# F.Bernhard 2022-02-02

rm(list=ls())
# install.packages("LWFBrook90R")
library(LWFBrook90R)
# packageVersion("LWFBrook90R")
# version
# sessionInfo()

library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)


# Load function that reads .jl input files to run it with LWFBrook90R
source("fct-run_LWFBrook90R_with_LWFBrook90jl_Input.R")
simulate_LWFBrook90.jl_input_folder_withLWFBrook90R

debug(simulate_LWFBrook90.jl_input_folder_withLWFBrook90R)
res1 <- simulate_LWFBrook90.jl_input_folder_withLWFBrook90R(folder="../test-assets/BEA-2016/input-files/", prefix="BEA2016-reset-FALSE", 
                                                            refine_discr_factor = 1, verbose = TRUE)
res2 <- simulate_LWFBrook90.jl_input_folder_withLWFBrook90R(folder="../test-assets/BEA-2016/input-files/", prefix="BEA2016-reset-FALSE", 
                                                            refine_discr_factor = 2, verbose = TRUE)
res3 <- simulate_LWFBrook90.jl_input_folder_withLWFBrook90R(folder="../test-assets/BEA-2016/input-files/", prefix="BEA2016-reset-FALSE", 
                                                            refine_discr_factor = 3, verbose = TRUE)
res10 <- simulate_LWFBrook90.jl_input_folder_withLWFBrook90R(folder="../test-assets/BEA-2016/input-files/", prefix="BEA2016-reset-FALSE", 
                                                            refine_discr_factor = 10, verbose = TRUE)

resB1 <- simulate_LWFBrook90.jl_input_folder_withLWFBrook90R(folder="../test-assets/DAV-2020/input-files/", prefix="DAV_LW1_def", 
                                                            refine_discr_factor = 1, verbose = TRUE)
resB2 <- simulate_LWFBrook90.jl_input_folder_withLWFBrook90R(folder="../test-assets/DAV-2020/input-files/", prefix="DAV_LW1_def", 
                                                            refine_discr_factor = 2, verbose = TRUE)
resB3 <- simulate_LWFBrook90.jl_input_folder_withLWFBrook90R(folder="../test-assets/DAV-2020/input-files/", prefix="DAV_LW1_def", 
                                                            refine_discr_factor = 3, verbose = TRUE)
resB10 <- simulate_LWFBrook90.jl_input_folder_withLWFBrook90R(folder="../test-assets/DAV-2020/input-files/", prefix="DAV_LW1_def", 
                                                             refine_discr_factor = 10, verbose = TRUE)

### Outputs
for (iter in list(list(res = res1, out_dir = "BEA-2016-IntegrationTests-LWFBrook90/output_LWFBrook90R", fname_base = "BEA2016-reset-FALSE"),
                  list(res = res2, out_dir = "BEA-2016-IntegrationTests-LWFBrook90/output_LWFBrook90R", fname_base = "BEA2016-reset-FALSE"),
                  list(res = res3, out_dir = "BEA-2016-IntegrationTests-LWFBrook90/output_LWFBrook90R", fname_base = "BEA2016-reset-FALSE"),
                  list(res = res10, out_dir = "BEA-2016-IntegrationTests-LWFBrook90/output_LWFBrook90R", fname_base = "BEA2016-reset-FALSE"),

                  list(res = resB1, out_dir = "DAV-2020-IntegrationTests-LWFBrook90/output_LWFBrook90R", fname_base = "DAV_LW1_def"),
                  list(res = resB2, out_dir = "DAV-2020-IntegrationTests-LWFBrook90/output_LWFBrook90R", fname_base = "DAV_LW1_def"),
                  list(res = resB3, out_dir = "DAV-2020-IntegrationTests-LWFBrook90/output_LWFBrook90R", fname_base = "DAV_LW1_def"),
                  list(res = resB10, out_dir = "DAV-2020-IntegrationTests-LWFBrook90/output_LWFBrook90R", fname_base = "DAV_LW1_def"))) {
  res <- iter$res
  out_dir    <- iter$out_dir    #"output_LWFBrook90R"
  fname_base <- iter$fname_base #"BEA2016-reset-FALSE"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  fnamePNG <- file.path(out_dir,
                        sprintf("%s_NLAYER%d_LWFBrook90R-%s.png", 
                                fname_base,
                                nrow(res$simin$soil_nodes),
                                packageVersion("LWFBrook90R")))
  fnameCSV <- file.path(out_dir,
                        sprintf("%s_NLAYER%d_LWFBrook90R-%s.csv", 
                                fname_base,
                                nrow(res$simin$soil_nodes),
                                packageVersion("LWFBrook90R")))
  dir.create(out_dir, showWarnings = F, recursive = T)
  
  # Output *.csvs
  round_digits <- 3
  for (result_name in list("daily_output","layer_output")) {
    if (result_name == "daily_output") {
      res_to_output <- res$simout[[result_name]] %>% tibble()
    } else {
      res_to_output <- res$simout[[result_name]] %>% tibble() %>%
        # append midpoint, upper, lower
        left_join(res$simin$soil_nodes %>%
                    mutate(lower = -0 + -cumsum(thick)) %>%
                    mutate(upper = c(0, lower[-n()])) %>%
                    select(nl = layer, midpoint, upper, lower), 
                  by = "nl")  
    }
    
    # Save RDS for internal use with R
    saveRDS(res_to_output,
            file = gsub(".csv", paste0(result_name, ".RDS") ,fnameCSV))
    
    # Save CSV for external use with Julia (unit tests in test suite)
    res_to_output %>%
      mutate(across(where(is.numeric), round, round_digits)) %>%
      mutate(across(!yr & !mo & !da & !doy, round, round_digits)) %>%
      write.csv(file = gsub(".csv", paste0(result_name, ".csv") ,fnameCSV))
  }

  # Plotting
  pl0_aboveground <-
    select(res$simout$daily_output, yr, mo, da, doy, rfal, sfal, #prec, # rfal+sfal instead of prec
           # rthr,
           flow, evap, seep, snow, swat, gwat, intr, ints, sint, isvp) %>%
    select(-doy) %>%
    tidyr::unite(col = "date", yr, mo, da) %>%
    mutate(date = lubridate::ymd(date)) %>%
    tidyr::pivot_longer(cols = -date,
                        # keep facets in ggplot in the same order as in data.frame
                        names_transform = list(name = forcats::fct_inorder)) %>%
    ggplot(aes(x=date, y=value, color = name)) +
    geom_line() + theme_bw() + facet_grid(name~., scales = "free_y") +
    theme(#strip.text = element_text(size=5),
      axis.text.y = element_text(size=5),
      legend.key.height = unit(1,"line")) +
    labs(x=NULL)
  
  pl1_belowground_theta <-
    {select(res$simout$layer_output, yr, mo, da, nl, theta) %>%
        tidyr::unite(col = "date", yr, mo, da) %>%
        mutate(date = lubridate::ymd(date)) %>%
        left_join({res$simin$soil_nodes %>%
            mutate(lower = -0 + -cumsum(thick)) %>%
            mutate(upper = c(0, lower[-n()])) %>%
            select(nl = layer, midpoint, upper, lower)}, 
            by = "nl") } %>%
    ggplot(aes(x=date, y=midpoint, fill = theta)) +
    geom_rect(aes(xmin=date, xmax=date+1, ymin=lower, ymax=upper)) +
    theme_bw() + labs(fill = "θ (-)", y = "Depth (m)", x=NULL)
  pl1_belowground_psimi <-
    {select(res$simout$layer_output, yr, mo, da, nl, psimi) %>%
        tidyr::unite(col = "date", yr, mo, da) %>%
        mutate(date = lubridate::ymd(date)) %>%
        left_join({res$simin$soil_nodes %>%
            mutate(lower = -0 + -cumsum(thick)) %>%
            mutate(upper = c(0, lower[-n()])) %>%
            select(nl = layer, midpoint, upper, lower)}, 
            by = "nl") } %>%
    ggplot(aes(x=date, y=midpoint, fill = psimi)) +
    geom_rect(aes(xmin=date, xmax=date+1, ymin=lower, ymax=upper)) +
    theme_bw() + labs(fill = "ψ_m (kPa)", y = "Depth (m)", x=NULL)
  
  depths_to_read_out <- -c(100,500,1000,1500,1900)/1000
  idx_to_read_out <- c()
  for (it in seq_along(depths_to_read_out)) {
    idx_to_read_out[it] <- which.min(abs(depths_to_read_out[it] - res$simin$soil_nodes$midpoint))
  }
  depths_actually_read_out <- 
    res$simin$soil_nodes$midpoint[idx_to_read_out]
  
  pl2_belowground_theta <- {select(res$simout$layer_output, yr, mo, da, nl, theta) %>%
      tidyr::unite(col = "date", yr, mo, da) %>%
      mutate(date = lubridate::ymd(date)) %>%
      left_join({res$simin$soil_nodes %>%
          mutate(lower = -0 + -cumsum(thick)) %>%
          mutate(upper = c(0, lower[-n()])) %>%
          select(nl = layer, midpoint, upper, lower)}, 
          by = "nl") } %>%
    filter(nl %in% c(idx_to_read_out)) %>%
    ggplot(aes(x=date, y=theta, color = as.factor(midpoint))) +
    geom_line() + theme_bw() +
    theme_bw() + labs(y = "θ (-)", color = "Depth (m)", x=NULL) + theme(legend.position = c(0.87, 0.25))
  pl2_belowground_psi<- {select(res$simout$layer_output, yr, mo, da, nl, psimi) %>%
      tidyr::unite(col = "date", yr, mo, da) %>%
      mutate(date = lubridate::ymd(date)) %>%
      left_join({res$simin$soil_nodes %>%
          mutate(lower = -0 + -cumsum(thick)) %>%
          mutate(upper = c(0, lower[-n()])) %>%
          select(nl = layer, midpoint, upper, lower)}, 
          by = "nl") } %>%
    filter(nl %in% c(idx_to_read_out)) %>%
    ggplot(aes(x=date, y=psimi, color = as.factor(midpoint))) +
    geom_line() + theme_bw() +
    theme_bw() + labs(y = "ψ_m (kPa)", color = "Depth (m)", x=NULL) + theme(legend.position = c(0.87, 0.25))
  
  # Output the plots
  ggsave(filename = fnamePNG,
         # width = 2200, height = 1800, units = "px", dpi = 220,
         # plot = grid.arrange(pl0_aboveground + ggtitle("Aboveground LWFBrook90R"),
         #                     pl1_belowground_theta + ggtitle("Belowground LWFBrook90R"),
         #                     pl1_belowground_psimi + ggtitle("Belowground LWFBrook90R"),
         #                     ncol=1,
         #                     layout_matrix = matrix(c(c(1,2),
         #                                              c(1,3)),
         #                                            ncol=2, byrow = T)))
         width = 2200, height = 1800, units = "px", dpi = 220,
         plot = grid.arrange(pl0_aboveground + ggtitle("Aboveground LWFBrook90R"),
                             pl1_belowground_theta + ggtitle("Belowground LWFBrook90R"),
                             pl1_belowground_psimi + ggtitle("Belowground LWFBrook90R"),
                             pl2_belowground_theta + ggtitle("Belowground LWFBrook90R"),
                             pl2_belowground_psi + ggtitle("Belowground LWFBrook90R"),
                             ncol=1,
                             layout_matrix = matrix(c(c(1,2,4),
                                                      c(1,3,5)),
                                                    ncol=3, byrow = T)))
}
