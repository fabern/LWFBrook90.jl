# User Guide



## Installation and usage

### Installation
First, download and install Julia: [https://julialang.org/downloads/](https://julialang.org/downloads/). If you like you can also install [julia-vscode.org](https://www.julia-vscode.org) to get a complete IDE. Both these steps are nicely summarized by https://techytok.com/julia-vscode/. ALternatively, you can also follow the official documentation of the VS Code Julia: [https://www.julia-vscode.org/docs/dev/gettingstarted/#Installation-and-Configuration-1](https://www.julia-vscode.org/docs/dev/gettingstarted/#Installation-and-Configuration-1).

To install LWFBrook90.jl open a Julia REPL, enter the Pkg REPL by pressing `]` and add the package:
```
(@v1.7) pkg> add LWFBrook90
(@v1.7) pkg> status
```
Dependencies of LWFBrook90.jl should automatically be installed. After that hit `Ctrl-C` to quit the Pkg REPL and return to the default Julia REPL. Now the installation can be tested with a single line `using LWFBrook90; run_example()`.

### Usage
Check out a step-by-step guide for a simulation in section [Example Script 01](@ref)

The steps in a typical simulation script are:
- load the package `using LWFBrook90`
- read input data `model = loadSPAC()`
    - optional arguments in `loadSPAC()` can be used to define model parameters (-> type `?loadSPAC` to see documentation)
- pre-process model for simulation `simulation = setup(model)`
- compute simulation `simulate!(simulation)`
- post-process and plot simulation results

LWFBrook90.jl outputs quantities in daily resolution.
Monthly or yearly quantities can be derived in the post processing.
~~For performance reasons, e.g. for Bayesian parameter estimation, computation of these additional quantities can also be deactivated during simulation, and they could be calculated in a post-processing step from the state vector.~~

## Input data

### Overview of input data
To run a simulation following input files are expected in a single folder
- `soil_horizons.csv` - containing the hydraulic parameters of the different soil horizons
- `param.csv` - containing scalar model parameters
- `meteoveg.csv` - containing daily values of meteorologic variables [and stand properties]
- `meteoiso.csv` - containing isotopic signatures of aggregate precipitation samples
- `initial_conditions.csv` - containing initial conditions of scalar state variables
- `soil_discretization.csv` - containing the initial conditions of the soil water status in the form of the soil matric potential (kPa) and the initial isotopic signatures (vector state variables); continuously defined parameters of relative root density distributions; as well as the definition of the numerical discretization of the soil domain (nodes with upper and lower limits).
- `meteo_storm_durations.csv` - containing parameters of sub-daily storm/precipitation event patterns for each month

Of these only the CSV files for the parameters `soil_horizons.csv` and `param.csv` and the forcings `meteoveg.csv`, `meteoiso.csv` are absolutely needed. The remaing can be provided by the user as arguments to `loadSPAC()` in the Julia script:
- `soil_discretizations.csv`: not needed if `Δz_thickness_m`, `root_distribution`, and `IC_soil` provided.
- `initial_conditions.csv`: not needed if `IC_scalar` provided.
- `meteo_storm_durations.csv`: not needed if `storm_durations_h` provided.
- `meteoveg.csv`: columns with vegetation parameters are not needed if `loadSPAC(canopy_evolution = ...)` is provided.

The structure of the input CSV's is illustrated by the example input data sets `isoBEA2010-18-*` or `DAV2020-bare-minimum` or `DAV2020-full` located in the folder `examples/`, as well as below in this documentation. Please follow these examples closely when generating your own input files, including the exact column names and header lines containing the units.

For convenience, input CSV files can be generated from a script that sets up a simulation with the R package [LWFBrook90R (v0.4.3)](https://github.com/pschmidtwalter/LWFBrook90R#usage). Instead of running the simulation with `run_LWFB90()`, the same arguments can be used to generate the input files for LWFBrook90.jl using the R function provided in the file `generate_LWFBrook90jl_Input.R`. Note that the input file `meteoiso.csv` needs to be generated separately and the files containing the initial conditions (`initial_conditions.csv` and `soil_discretization.csv`) also need to be extended manually with the isotope values (see structure of these input files below).

To load input data and prepare a simulation follow the instructions in section [Example Script 01](@ref) or alternatively use the sample script `main_with_isotopes.jl`. NOTE: these will be replaced with Jupyter-notebooks generated with Literate.jl.

In case you're unfamiliar with Julia, there are various ways to run a script such as `main.jl`: One possibility is to open the Julia REPL and run the script using `include(“main.jl”)`. Alternatively, the editor VS Code in combination with the Julia extension ([julia-vscode.org](https://www.julia-vscode.org)), provides a complete IDE for programming in Julia.

### Structure of input data
```@setup read_csv_inputs
# from: https://juliadocs.github.io/Documenter.jl/stable/man/syntax/#@eval-block
using CSV
using DataFrames
using Latexify

example_soil_horizons         = CSV.read("../../examples/isoBEA2010-18-reset-FALSE-input/isoBEA2010-18-reset-FALSE_soil_horizons.csv",          DataFrame)
example_meteoveg              = CSV.read("../../examples/isoBEA2010-18-reset-FALSE-input/isoBEA2010-18-reset-FALSE_meteoveg.csv",               DataFrame)
example_meteoiso              = CSV.read("../../examples/isoBEA2010-18-reset-FALSE-input/isoBEA2010-18-reset-FALSE_meteoiso.csv",               DataFrame)
example_meteo_storm_durations = CSV.read("../../examples/isoBEA2010-18-reset-FALSE-input/isoBEA2010-18-reset-FALSE_meteo_storm_durations.csv",  DataFrame)
example_initial_conditions    = CSV.read("../../examples/isoBEA2010-18-reset-FALSE-input/isoBEA2010-18-reset-FALSE_initial_conditions.csv",     DataFrame)
example_soil_discretization   = CSV.read("../../examples/isoBEA2010-18-reset-FALSE-input/isoBEA2010-18-reset-FALSE_soil_discretization.csv",    DataFrame)
example_param                 = CSV.read("../../examples/isoBEA2010-18-reset-FALSE-input/isoBEA2010-18-reset-FALSE_param.csv",                  DataFrame)

# Fix underscores
# Fix row title underscores
rename!(example_soil_horizons         , replace.(names(example_soil_horizons         ), "_" => "\\_"))
rename!(example_meteoveg              , replace.(names(example_meteoveg              ), "_" => "\\_"))
rename!(example_meteoiso              , replace.(names(example_meteoiso              ), "_" => "\\_"))
rename!(example_meteo_storm_durations , replace.(names(example_meteo_storm_durations ), "_" => "\\_"))
rename!(example_initial_conditions    , replace.(names(example_initial_conditions    ), "_" => "\\_"))
rename!(example_soil_discretization   , replace.(names(example_soil_discretization   ), "_" => "\\_"))
rename!(example_param                 , replace.(names(example_param                 ), "_" => "\\_"))
# Fix table content underscores
# Fix table content underscores
example_soil_horizons          = combine(example_soil_horizons         , All() .=> (x)->replace.(x, "_" => "\\_"), renamecols=false)
example_meteoveg               = combine(example_meteoveg              , All() .=> (x)->replace.(x, "_" => "\\_"), renamecols=false)
example_meteoiso               = combine(example_meteoiso              , All() .=> (x)->replace.(x, "_" => "\\_"), renamecols=false)
example_meteo_storm_durations  = combine(example_meteo_storm_durations , All() .=> (x)->replace.(x, "_" => "\\_"), renamecols=false)
example_initial_conditions     = combine(example_initial_conditions    , All() .=> (x)->replace.(x, "_" => "\\_"), renamecols=false)
example_soil_discretization    = combine(example_soil_discretization   , All() .=> (x)->replace.(x, "_" => "\\_"), renamecols=false)
example_param                  = combine(example_param                 , All() .=> (x)->replace.(x, "_" => "\\_"), renamecols=false)

n_half = 3  # hide
example_meteoveg_to_print = [push!(copy(first(example_meteoveg, n_half+1)),
    fill("...",size(example_meteoveg,2)))
    last(example_meteoveg, n_half)]  # hide
example_meteoiso_to_print = [push!(copy(first(example_meteoiso, n_half+1)),
    fill("...",size(example_meteoiso,2)))
    last(example_meteoiso, n_half)]  # hide
```

`soil_horizons.csv`: (when using Mualem-VanGenuchten parametrization of the soil retention curve):
```@example read_csv_inputs
mdtable(example_soil_horizons,latex=false) # hide
```


`meteoveg.csv`: contains time dependent parameters (meterologic variables and vegetation parameters):
```@example read_csv_inputs
mdtable(example_meteoveg_to_print , latex=false)  # hide
```

If meteoveg.csv does not contain columns for DENSEF, HEIGHT, LAI, SAI, the parametrization
of these should be provided to the function `loadSPAC()`.

`meteoiso.csv`: contains time dependent isotopic signatures of precipiation.
(Note that the `dates` contain the end dates of the collection interval of cumulative
isotope samples. The values are backward interpolated with piecewise constants. The first
row containng the `NA` is assumed to be the start date of the first collection interval.):

```@example read_csv_inputs
mdtable(example_meteoiso_to_print , latex=false)  # hide
```


`meteo_storm_durations.csv`:
```@example read_csv_inputs
mdtable(example_meteo_storm_durations , latex=false) # hide
```


`initial_conditions.csv`:
```@example read_csv_inputs
mdtable(example_initial_conditions    , latex=false) # hide
```


`soil_discretization.csv` contains the initial conditions of the soil water status, root density distributions, and the definition of the numerical discretization of the soil domain (nodes with upper and lower limits):
```@example read_csv_inputs
mdtable(example_soil_discretization   , latex=false) # hide
```


`param.csv` contains scalar model parameters:
```@example read_csv_inputs
mdtable(example_param                 , latex=false) # hide
```

## Calibration data (calibration not yet implemented)
Roadmap to include calibration data intends include of:
- Throughfall amounts (for parametrisation of interception)
- Soil moisture (volumetric water content, θ)
- Soil matric potential (ψ)

as well as

- Isotopic composition of soil water (δ_soil)
- Isotopic composition of xylem water (δ_xylem)

Below the example structure of data sets for
`psi.csv` contains the soil matric potential (ψ) of different sensor series or as site average:

| dates      | depth_cm | depth_nominal_cm | psi_kPa | series  |
| ---------- | -------- | ---------------- | ------- | ------- |
| YYYY-MM-DD | cm       | cm               | kPa     | SITEAVG |
| 2022-08-30 | 15       | 15               | -5.77   | SITEAVG |
| 2022-08-30 | 200      | 150              | -2.87   | SITEAVG |
| 2022-08-30 | 300      | 150              | -0.10   | SITEAVG |
| 2022-08-30 | 50       | 50               | -1.98   | SITEAVG |
| 2022-08-31 | 80       | 80               | -0.13   | SITEAVG |
| 2022-08-31 | 15       | 15               | -7.53   | SITEAVG |
| 2022-08-31 | 200      | 150              | -2.95   | SITEAVG |
| 2022-08-31 | 300      | 150              | -0.10   | SITEAVG |
| 2021-06-17 | 50       | 50               | -3.00   | SITEAVG |
| 2021-06-17 | 80       | 80               | -0.17   | SITEAVG |
| ...        | ...      | ...              | ...     | ...     |
|            |          |                  |         |         |

`theta.csv` contains the soil moisture (volumetric water content, θ):

| dates      | depth_m | theta_m3m3 |
| ---------- | ------- | ---------- |
| YYYY-MM-DD | m       | m3/m3      |
| 2020-12-23 | 0.15    | 0.24       |
| 2020-12-23 | 0.50    | 0.32       |
| 2020-12-23 | 0.80    | 0.30       |
| 2020-12-24 | 0.15    | 0.24       |
| 2020-12-24 | 0.50    | 0.32       |
| ...        | ...     | ...        |


`delta_soil.csv` contains the isotopic signature of soil water:

| dates      | depth_m | delta18O_permil | delta2h_permil |
| ---------- | ------- | -----------     | -------        |
| YYYY-MM-DD | cm      | permil          | permil         |
| 2020-07-07 | 0       | -5.89           | -35.96         |
| 2020-07-07 | 15      | -8.75           | -57.99         |
| 2020-07-07 | 50      | -8.81           | -59.12         |
| 2020-07-07 | 80      | -9.72           | -65.41         |
| 2020-07-21 | 0       | -5.83           | -29.54         |
| 2020-07-21 | 15      | -8.73           | -53.58         |
| 2020-07-21 | 50      | -9.31           | -58.47         |
| 2020-07-21 | 80      | -10.14          | -64.71         |
| ...        | ...     | ...             | ...            |

`delta_xylem.csv` contains the isotopic signature of xylem water:

| dates      | species | treeID | delta18O_permil | delta2H_permil |
| ---------- | ------- | ------ | --------------- | -------------- |
| YYYY-MM-DD | -       | -      | permil          |permil          |
| 2021-08-03 | Beech   | 271    | -8.416          | -72.562        |
| 2021-08-03 | Beech   | 3518   | -8.946          | -76.543        |
| 2021-08-03 | Beech   | 3519   | -8.484          | -72.523        |
| 2021-08-17 | Beech   | 271    | -8.378          | -76.044        |
| 2021-08-17 | Beech   | 3518   | -8.769          | -76.036        |
| 2021-08-17 | Beech   | 3519   | -8.135          | -70.475        |
| ...        | ...     | ...    | ...             |...             |
