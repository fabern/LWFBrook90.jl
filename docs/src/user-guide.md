# User Guide



## Installation and usage

### Installation
First, download and install Julia: [https://julialang.org/downloads/](https://julialang.org/downloads/). If you like you can also install [julia-vscode.org](https://www.julia-vscode.org) to get a complete IDE. Both these steps are nicely summarized by https://techytok.com/julia-vscode/. ALternatively, you can also follow the official documentation of the VS Code Julia: [https://www.julia-vscode.org/docs/dev/gettingstarted/#Installation-and-Configuration-1](https://www.julia-vscode.org/docs/dev/gettingstarted/#Installation-and-Configuration-1).

To install LWFBrook90.jl open a Julia REPL, enter the Pkg REPL by pressing `]` and add the package:
```
(@v1.5) pkg> add LWFBrook90
(@v1.5) pkg> status
```
Dependencies of LWFBrookJulia.jl should automatically be installed. After that hit `Ctrl-C` to quit the Pkg REPL and return to the default Julia REPL.

### Usage
Check out a step-by-step in guide for a simulation in section [Example](@ref)

The steps in a typical simulation script are:
- load the package `using LWFBrook90`
- load the package dependency `using OrdinaryDiffEq: solve, Tsit5`
- read input data
- set up model options
- set up an `ODE` problem (`u0`,`tspan`, `p`) and solve it with DifferentialEquations.jl
- plot and/or postprocess simulation results

LWFBrook90.jl can output additional quantities in daily resolution derived during simulation. Monthly or yearly quantities can be derived in the post processing.
For performance reasons, e.g. for Bayesian parameter estimation, computation of these additional quantities can also be deactivated during simulation, and they could be calculated in a post-processing step from the state vector (not yet implemented).



## Input data

### Overview of input data
To run a simulation following input data are needed
- `soil_horizons.csv` - containing the hydraulic parameters of the different soil horizons
- `meteoveg.csv` - containing daily values of meteorologic variables and stand properties
- `meteo_storm_durations.csv` - containing parameters of sub-daily storm/precipitation event patterns
- `initial_conditions.csv` - containing initial conditions of scalar state variables
- `soil_discretization.csv` - containing the initial conditions of the soil water status in the form of the soil matric potential (kPa) (vector state variable); continuously defined parameters of relative root density distributions; as well as the definition of the numerical discretization of the soil domain (nodes with upper and lower limits).
- `param.csv` - containing further scalar model parameters

The structure of the input data is illustrated by the example input data set `BEA2016-*` located in the folder `examples/`, as well as below in this documentation. Please follow these examples closely when generating your own input files, including the exact column names and header lines containing the units.

For convenience, input files can be generated from a script that sets up a simulation with the R package [LWFBrook90R (v0.4.3)](https://github.com/pschmidtwalter/LWFBrook90R#usage). Instead of running the simulation with `run_LWFB90()`, the same arguments can be used to generate the input files for LWFBrook90.jl using the R function provided in the file `generate_LWFBrook90jl_Input.R`.

To load load input data and prepare a simulation follow the instructions in section [Example](@ref) or alternatively use the sample script `main.jl`.

In case you're unfamiliar to Julia, there are various ways to run a script such as `main.jl`: One possibility is to open the Julia REPL and run the script using `include(“main.jl”)`. Alternatively, the editor VS Code in combination with the Julia extension ([julia-vscode.org](https://www.julia-vscode.org)), provides a complete IDE for programming in Julia.

### Structure of input data

`soil_horizons.csv`: (when using Mualem-van Genuchten parametrization)

| HorizonNr | Upper_m | Lower_m | ths_volFrac         | thr_volFrac         | alpha_perMeter | npar_   | ksat_mmDay | tort_   | gravel_volFrac |
| --------- | ------- | ------- | ------------------- | ------------------- | -------------- | ------- | ---------- | ------- | -------------- |
| -         | m       | m       | volume fraction (-) | volume fraction (-) | perMeter       | -       | mm per day | -       | volFrac        |
| 1         | 0       | -0.04   | 0.7149              | 0.069               | 1147.88919     | 1.05123 | 24864.6784 | 4.67037 | 0.01           |
| 2         | -0.04   | -0.08   | 0.66873             | 0.069               | 1274.88602     | 1.05105 | 12881.4864 | 4.47828 | 0.175          |
| 3         | -0.08   | -0.2    | 0.6561              | 0.069               | 1215.92721     | 1.05106 | 10516.6166 | 4.50162 | 0.175          |
| 4         | -0.2    | -0.45   | 0.60159             | 0.069               | 795.20401      | 1.05099 | 3438.06275 | 4.34489 | 0.175          |
| 5         | -0.45   | -0.75   | 0.52331             | 0.069               | 352.36825      | 1.05097 | 450.35802  | 4.29122 | 0.01           |
| 6         | -0.75   | -1.1    | 0.56472             | 0.069               | 570.68168      | 1.05111 | 1488.20958 | 5.14818 | 0.01           |
| 7         | -1.1    | -1.2    | 0.46743             | 0.069               | 164.564        | 1.05103 | 67.97846   | 0.01    | 0.95           |


`soil_horizons.csv`: (when using Clapp-Hornberger parametrization) (NOT COMPLETELY IMPLEMENTED!)

| HorizonNr  | Upper_m | Lower_m | thsat_volFrac       | thetaf_volFrac      | psif_kPa    | bexp_  | kf_mmDay     | wtinf_ | gravel_volFrac      |
| ---------- | ------- | ------- | ------------------- | ------------------- | ----------- | ------ | ------------ | ------ | ------------------- |
| -          | m       | m       | volume fraction (-) | volume fraction (-) | kPa         | -      | mm per day   | -      | volume fraction (-) |
| 1          | 0       | -0.04   | NA                  | NA                  | NA          | NA     | NA           | NA     | NA                  |


`meteoveg.csv`: contains time dependent parameters (meterologic variables and vegetation parameters):

| dates      | globrad_MJDayM2 | tmax_degC | tmin_degC | vappres_kPa | windspeed_ms | prec_mmDay | densef | height | lai   | sai  | age       |
| ---------- | --------------- | --------- | --------- | ----------- | ------------ | ---------- | ------ | ------ | ----- | ---- | --------- |
| YYYY-MM-DD | MJ/Day/m2       | degree C  | degree C  | kPa         | m per s      | mm per day | -      | m      | -     | -    | years     |
| 01.01.10   | 4.75            | -0.68     | -5.2      | 0.4         | 1.27         | 0.6        | 1      | 25     | 2.274 | 1    | 240       |
| 02.01.10   | 4.64            | -3.95     | -14.57    | 0.22        | 3.32         | 0          | 1      | 25     | 2.274 | 1    | 240.00275 |
| 03.01.10   | 6.58            | -8.28     | -16.07    | 0.13        | 1.59         | 0          | 1      | 25     | 2.274 | 1    | 240.00549 |
| ...        | ...             | ...       | ...       | ...         | ...          | ...        | ...    | ...    | ...   | ...  | ...       |
| 21.06.10   | 10.94           | 5.18      | 2.03      | 0.72        | 2.29         | 0.9        | 1      | 25     | 3.79  | 1    | 240.46978 |
| 22.06.10   | 8.17            | 5.62      | 2.85      | 0.77        | 2.17         | 0.1        | 1      | 25     | 3.79  | 1    | 240.47253 |
| 23.06.10   | 30.39           | 13.53     | 1.68      | 0.81        | 2.66         | 0          | 1      | 25     | 3.79  | 1    | 240.47527 |
| 24.06.10   | 30.57           | 16.13     | 3.95      | 0.88        | 2.28         | 0          | 1      | 25     | 3.79  | 1    | 240.47802 |

`meteo_storm_durations.csv`:

| month                                                        | average_storm_duration_h |
| ------------------------------------------------------------ | ------------------------ |
| ### Typical average durations of a single  storm event for each month ------- | NA                       |
| January                                                      | 4                        |
| Februrary                                                    | 4                        |
| March                                                        | 4                        |
| April                                                        | 4                        |
| May                                                          | 4                        |
| June                                                         | 4                        |
| July                                                         | 4                        |
| August                                                       | 4                        |
| September                                                    | 4                        |
| October                                                      | 4                        |
| November                                                     | 4                        |
| December                                                     | 4                        |
|                                                              |                          |

`initial_conditions.csv`:

| param_id                                                            | x    |
| ------------------------------------------------------------------- | ---- |
| ### Initial conditions (except for  depth-dependent u_aux_PSIM) --- | NA   |
| u_GWAT_init_mm                                                      | 0    |
| u_INTS_init_mm                                                      | 0    |
| u_INTR_init_mm                                                      | 0    |
| u_SNOW_init_mm                                                      | 0    |
| u_CC_init_MJ_per_m2                                                 | 0    |
| u_SNOWLQ_init_mm                                                    | 0    |

`soil_discretization.csv` contains the initial conditions of the soil water status, root density distributions, and the definition of the numerical discretization of the soil domain (nodes with upper and lower limits):

| Upper_m | Lower_m | Rootden_ | uAux_PSIM_init_kPa | u_delta18O_init_mUr | u_delta2H_init_mUr |
| ------- | ------- | -------- | ------------------ | ------------------- | ------------------ |
| m       | m       | -        | kPa                | mUr                 | mUr                |
| 0       | -0.04   | 0.02868  | -6.3               | NA                  | NA                 |
| -0.04   | -0.08   | 0.02539  | -6.3               | NA                  | NA                 |
| -0.08   | -0.2    | 0.02     | -6.3               | NA                  | NA                 |
| -0.2    | -0.45   | 0.01159  | -6.3               | NA                  | NA                 |
| -0.45   | -0.75   | 0.00507  | -6.3               | NA                  | NA                 |
| -0.75   | -1.1    | 0.00191  | -6.3               | NA                  | NA                 |
| -1.1    | -1.2    | 0.00092  | -6.3               | NA                  | NA                 |


`param.csv` contains scalar model parameters:

| param_id                                                     | x        |
| ------------------------------------------------------------ | -------- |
| **### Meteorologic site parameters -------**                 | NA       |
| LAT_DEG                                                      | 46.70052 |
| ESLOPE_DEG                                                   | 18.26    |
| ASPECT_DEG                                                   | 225      |
| ALB                                                          | 0.2      |
| ALBSN                                                        | 0.5      |
| C1                                                           | 0.25     |
| C2                                                           | 0.5      |
| C3                                                           | 0.2      |
| WNDRAT                                                       | 0.3      |
| FETCH                                                        | 5000     |
| Z0W                                                          | 0.005    |
| ZW                                                           | 2        |
| **### Canopy parameters -------**                            | NA       |
| LWIDTH                                                       | 0.1      |
| Z0G                                                          | 0.00325  |
| Z0S                                                          | 0.001    |
| LPC                                                          | 4        |
| CS                                                           | 0.035    |
| CZS                                                          | 0.13     |
| CZR                                                          | 0.05     |
| HS                                                           | 1        |
| HR                                                           | 10       |
| ZMINH                                                        | 2        |
| RHOTP                                                        | 2        |
| NN                                                           | 2.5      |
| **### Interception parameters -------**                      | NA       |
| FRINTLAI                                                     | 0.06     |
| FSINTLAI                                                     | 0.04     |
| FRINTSAI                                                     | 0.06     |
| FSINTSAI                                                     | 0.04     |
| CINTRL                                                       | 0.15     |
| CINTRS                                                       | 0.15     |
| CINTSL                                                       | 0.6      |
| CINTSS                                                       | 0.6      |
| RSTEMP                                                       | -0.5     |
| **### Snowpack parameters -------**                          | NA       |
| MELFAC                                                       | 1.5      |
| CCFAC                                                        | 0.3      |
| LAIMLT                                                       | 0.2      |
| SAIMLT                                                       | 0.5      |
| GRDMLT                                                       | 0.35     |
| MAXLQF                                                       | 0.05     |
| KSNVP                                                        | 0.3      |
| SNODEN                                                       | 0.3      |
| **### Leaf evaporation parameters  (affecting PE) -------**  | NA       |
| GLMAX                                                        | 0.0053   |
| GLMIN                                                        | 0.0003   |
| CR                                                           | 0.5      |
| RM                                                           | 1000     |
| R5                                                           | 100      |
| CVPD                                                         | 2        |
| TL                                                           | 0        |
| T1                                                           | 10       |
| T2                                                           | 30       |
| TH                                                           | 40       |
| **### Plant parameters (affecting  soil-water supply) -------** | NA       |
| MXKPL                                                        | 8        |
| MXRTLN                                                       | 3000     |
| INITRLEN                                                     | 12       |
| INITRDEP                                                     | 0.25     |
| RGRORATE                                                     | 0.03     |
| RGROPER                                                      | 30       |
| FXYLEM                                                       | 0.5      |
| PSICR                                                        | -2       |
| RTRAD                                                        | 0.35     |
| NOOUTF                                                       | 1        |
| **### Soil parameters -------**                              | NA       |
| IDEPTH_m                                                     | 0.040    |
| QDEPTH_m                                                     | 0.0      |
| RSSA                                                         | 100      |
| RSSB                                                         | 1        |
| INFEXP                                                       | 0        |
| BYPAR                                                        | 0        |
| QFPAR                                                        | 1        |
| QFFC                                                         | 0        |
| IMPERV                                                       | 0        |
| DSLOPE                                                       | 0        |
| LENGTH_SLOPE                                                 | 200      |
| DRAIN                                                        | 1        |
| GSC                                                          | 0        |
| GSP                                                          | 0        |
| **### Numerical solver parameters -------**                  | NA       |
| DTIMAX                                                       | 0.5      |
| DSWMAX                                                       | 0.05     |
| DPSIMAX                                                      | 0.0005   |

## Calibration data (calibration not yet implemented)
Roadmap to include calibration data intends to allow inclusion of:
- Throughfall amounts (for parametrisation of interception)
- Soil moisture (volumetric water content, θ)
- Soil matric potential (ψ)
