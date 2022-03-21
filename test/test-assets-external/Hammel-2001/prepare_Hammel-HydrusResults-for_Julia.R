# F.Bernhard 2022-02-02

rm(list=ls())

library(tidyr)
library(dplyr)
library(LWFBrook90R)
library(ggplot2)
library(gridExtra)

# Prepare Hydrus results for treatment with Julia
conc2delta18O <- function(conc){
  #round(1000*(conc/(2005.2e-6*(1-conc))-1),2)
  if (all(is.character(conc))) {
    return(NA_real)
  } else {
    return(conc-15)
  }
}
conc2delta2H <- function(conc){
  #round(delta <- 1000*(conc/(155.76e-6*(1-conc))-1),2)
  # return(conc-120)
  if (all(is.character(conc))) {
    return(NA_real)
  } else {
    return(conc-120)
  }
}

# Hammel_select <- "Hammel_loam"
for (Hammel_select in c("Hammel_loam", "Hammel_sand")) {
  # hydrus_folder <- if_else(Hammel_select == "Hammel_loam",
  #                          "Hammel-IntegrationTests-Hydrus/Hammel_Test_Loam",
  #                          "Hammel-IntegrationTests-Hydrus/Hammel_Test_Sand")
  hydrus_folder <- if_else(Hammel_select == "Hammel_loam",
                           "Hammel-IntegrationTests-Hydrus_ISO/Hammel_Test_Loam_ISO2",
                           "Hammel-IntegrationTests-Hydrus_ISO/Hammel_Test_Sand_ISO2")
  hydrus_file_NodInf <- file.path(hydrus_folder, "Nod_Inf.out")
  hydrus_file_Atmos <- file.path(hydrus_folder, "ATMOSPH.IN")
  hydrus_file_Selector <- file.path(hydrus_folder, "SELECTOR.IN")

  #----------- gather information from Selector.in----------------------------------
  rawSelector <- readLines(hydrus_file_Selector)
  units <- data.frame(transpose(list(
    rawSelector[grep("LUnit",rawSelector)+(1:3)]))) %>%
    rename("LengthUnit"=1,"TimeUnit"=2,"ConcentrationMassUnit"=3)
  unit2cm <- ifelse(units$LengthUnit=="mm",0.1,1)

  #----------- gather information from ATMOSPH.IN----------------------------------
  rawAtmos <- readLines(hydrus_file_Atmos)

  nCycleLine <- grep("Number of Cycles",rawAtmos)+1
  nCycles <- ifelse(length(nCycleLine) == 0,
                    1,
                    as.numeric(rawAtmos[nCycleLine]))

  startLine <- grep("Prec +rSoil",rawAtmos) + 1 # Header Line
  endLine   <- grep("end",rawAtmos) - 1         # End line
  ATMOSPH_DAT <- readr::read_table(
    rawAtmos[startLine:endLine],
    col_names = c("timeStep","prcp","Epot","Tpot","hCritA","rB","hB",
                  "PressureHead_top","Temp_top","Temp_bot","Ampl",
                  "conc18O_top","conc18O_bot","conc2H_top","conc2H_bot")) %>%
    mutate(delta18O_top = conc2delta18O(conc18O_top),
           delta2H_top  = conc2delta2H(conc2H_top))

  timeSteps <- NULL
  for(n in 1:nCycles){
    timeSteps <- c(timeSteps, ATMOSPH_DAT$timeStep+n-1)
  }

  ATMOSPH_DAT <- as.data.frame(lapply(ATMOSPH_DAT, rep, nCycles))

  ATMOSPH_DAT <- ATMOSPH_DAT %>%
    mutate(timeStepStart = c(0,timeSteps[-length(timeSteps)]),
           timeStepEnd   = timeSteps)

  #----------- post process Nod_Inf.out --------------------------------------------
  raw <- paste(readLines(hydrus_file_NodInf),collapse='\n')
  raw <- unlist(strsplit(raw, "Time:"))
  raw <- raw[-c(1,2)] # remove first two elements

  NOD_INF_DAT_list <- list()
  for(n in 1:length(raw)){ # cycle over nTimeSteps = length(raw)
    y <- readLines(textConnection(raw[n]))
    # each raw element `y` has the form:
    # [1] "        0.0000"
    # [2] ""
    # [3] ""
    # [4] " Node      Depth      Head Moisture       K          C         Flux        Sink         Kappa   v/KsTop   Temp   Conc(1..NS) Sorb(1...NS)"
    # [5] "           [L]        [L]    [-]        [L/T]      [1/L]      [L/T]        [1/T]         [-]      [-]      [C]      [M/L*3]"
    # [6] ""
    # [7] "   1     0.0000    -164.554 0.0603   0.8336E-03  0.9905E-04 -0.8336E-03  0.0000E+00      -1  -0.130E-05   20.00  0.5000E+01  0.5000E+02  0.0000E+00  0.0000E+00"
    # [8] "   2    -2.0000    -164.554 0.0603   0.8336E-03  0.9905E-04 -0.8336E-03  0.0000E+00      -1  -0.130E-05   20.00  0.5000E+01  0.5000E+02  0.0000E+00  0.0000E+00"
    # ...
    # [106] " 100  -198.0000    -164.554 0.0603   0.8336E-03  0.9905E-04 -0.8336E-03  0.0000E+00      -1  -0.130E-05   20.00  0.5000E+01  0.5000E+02  0.0000E+00  0.0000E+00"
    # [107] " 101  -200.0000    -164.554 0.0603   0.8336E-03  0.9905E-04 -0.8336E-03  0.0000E+00      -1  -0.130E-05   20.00  0.5000E+01  0.5000E+02  0.0000E+00  0.0000E+00"
    # [108] "end"
    # [109] ""
    # [110] ""
    # [111] " "

    curr_timeStep <- as.numeric(y[1])
    curr_header   <- y[4] # unused, as it is 13 elements long, but data can be 15 or more...
    curr_data     <- y[7:(grep("end",y)-1)]

    NOD_INF_DAT_list[[n]] <-
      readr::read_table(curr_data,
                        col_names = c("node","depth","head","theta", "K", "C",
                                      "Flux", "RWU", "Kappa", "vKsTop", "Temp",
                                      "Conc18O", "Sorb18O","Conc2H", "Sorb2H")) %>%
      mutate(delta18O = conc2delta18O(Conc18O),
             delta2H  = conc2delta2H(Conc2H),
             timeStep = curr_timeStep)
  }
  NOD_INF_DAT <- dplyr::bind_rows(NOD_INF_DAT_list) %>%
    # bring depth and pressure head to cm
    mutate(depth_cm = depth*unit2cm,
           head_cm  = head*unit2cm)


  write.csv(x = NOD_INF_DAT,
            file = file.path(dirname(hydrus_file_NodInf),
                             "Nod_inf_processed.csv"))

  # resultParSpecs <-
  #   data.frame(name        = c("head",           "theta",                    "delta18O",             "RWU",               "K",                      "Flux"),
  #              description = c("Hydraulic Head", "Volumetric Soil Moisture", "Soil Water delta 18O", "Root Water Uptake", "Hydraulic Conductivity", "Water Flux"),
  #              unit        = c("cm",             "-",                        "\u2030",               "1/d",               "cm/d",                   "cm/d"))

  pl_theta <- ggplot(NOD_INF_DAT,
         aes(x=timeStep, y=depth_cm, fill = theta)) +
    geom_tile() + scale_x_continuous(breaks = -100:100) + theme_bw()
  pl_18O <- ggplot(NOD_INF_DAT,
         aes(x=timeStep, y=depth_cm, fill = delta18O)) +
    geom_tile() + scale_x_continuous(breaks = -100:100) + theme_bw()
  pl_2H <- ggplot(NOD_INF_DAT,
         aes(x=timeStep, y=depth_cm, fill = delta2H)) +
    geom_tile() + scale_x_continuous(breaks = -100:100) + theme_bw()

  pl_theta_depths <- ggplot(filter(NOD_INF_DAT, depth_cm %in% c(-10,-100,-150, -190)),
         aes(x=timeStep, y=theta, color = as.factor(depth_cm), group = depth_cm)) +
    geom_line() + scale_x_continuous(breaks = -100:100) + theme_bw()
  pl_18O_depths <- ggplot(filter(NOD_INF_DAT, depth_cm %in% c(-10,-100,-150, -190)),
                          aes(x=timeStep, y=delta18O, color = as.factor(depth_cm), group = depth_cm)) +
    geom_line() + scale_x_continuous(breaks = -100:100) + theme_bw()

  pl_BC_18O <- ggplot(tibble(ATMOSPH_DAT),
                      aes(y=prcp,x=timeStepEnd,
                          xmin = timeStepStart, xmax=timeStepEnd, ymin=0, ymax=prcp)) +
    geom_rect(aes(fill=delta18O_top)) +
    scale_x_continuous(breaks = -100:100) + theme_bw()


  pl <- egg::ggarrange(
    ncol=1,
    pl_BC_18O,
    pl_theta_depths,
    pl_18O_depths,
    pl_18O)
  ggsave(plot = pl, file.path(dirname(hydrus_file_NodInf), "Nod_inf_processed.png"))
}






