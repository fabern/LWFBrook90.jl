# User Guide



## Installation and usage

### Installation
First, download and install Julia: [https://julialang.org/downloads/](https://julialang.org/downloads/). If you like you can also install [julia-vscode.org](https://www.julia-vscode.org) to get a complete IDE. Both these steps are nicely summarized by https://techytok.com/julia-vscode/. ALternatively, you can also follow the official documentation of the VS Code Julia: [https://www.julia-vscode.org/docs/dev/gettingstarted/#Installation-and-Configuration-1](https://www.julia-vscode.org/docs/dev/gettingstarted/#Installation-and-Configuration-1).

To install LWFBrook90.jl open a Julia REPL, enter the Pkg REPL by pressing `]` and add the package:
```
(@v1.7) pkg> add LWFBrook90
(@v1.7) pkg> status
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
For performance reasons, e.g. for Bayesian parameter estimation, computation of these additional quantities can also be deactivated during simulation, and they could be calculated in a post-processing step from the state vector.



## Input data

### Overview of input data
To run a simulation following input data are needed
- `soil_horizons.csv` - containing the hydraulic parameters of the different soil horizons
- `meteoveg.csv` - containing daily values of meteorologic variables and stand properties
- `meteoiso.csv` - containing isotopic signatures of aggregate precipitation samples
- `meteo_storm_durations.csv` - containing parameters of sub-daily storm/precipitation event patterns for each month
- `initial_conditions.csv` - containing initial conditions of scalar state variables
- `soil_discretization.csv` - containing the initial conditions of the soil water status in the form of the soil matric potential (kPa) and the initial isotopic signatures (vector state variables); continuously defined parameters of relative root density distributions; as well as the definition of the numerical discretization of the soil domain (nodes with upper and lower limits).
- `param.csv` - containing further scalar model parameters

The structure of the input data is illustrated by the example input data set `isoBEA2010-18-*` located in the folder `examples/`, as well as below in this documentation. Please follow these examples closely when generating your own input files, including the exact column names and header lines containing the units.

For convenience, input files can be generated from a script that sets up a simulation with the R package [LWFBrook90R (v0.4.3)](https://github.com/pschmidtwalter/LWFBrook90R#usage). Instead of running the simulation with `run_LWFB90()`, the same arguments can be used to generate the input files for LWFBrook90.jl using the R function provided in the file `generate_LWFBrook90jl_Input.R`. Note that the input file `meteoiso.csv` needs to be generated separately and the files containing the initial conditions (`initial_conditions.csv` and `soil_discretization.csv`) also need to be extended manually with the isotope values (see structure of these input files below).

To load load input data and prepare a simulation follow the instructions in section [Example](@ref) or alternatively use the sample script `main_with_isotopes.jl`. NOTE: these will shortly be replace with Jupyter-notebooks generated with Literate.jl (TODO).

In case you're unfamiliar to Julia, there are various ways to run a script such as `main.jl`: One possibility is to open the Julia REPL and run the script using `include(“main.jl”)`. Alternatively, the editor VS Code in combination with the Julia extension ([julia-vscode.org](https://www.julia-vscode.org)), provides a complete IDE for programming in Julia.

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

- isotopic composition of soil water (δ_soil)
- isotopic composition of xylem water (δ_xylem)

Below the example structure of data sets for
`psi.csv` contains the soil matric potential (ψ):

| dates      | depth_m | psi_kPa |
| ---------- | ------- | ------- |
| YYYY-MM-DD | m       | kPa     |
| 2019-12-23 | 0.00    | -0.1    |
| 2019-12-23 | 0.20    | -0.93   |
| 2019-12-23 | 0.40    | -0.1    |
| 2019-12-23 | 0.80    | -0.27   |
| 2019-12-23 | 1.60    | -0.1    |
| 2019-12-24 | 0.00    | -0.1    |
| ...        | ...     | ...     |
|            |         |         |

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

TODO...

`delta_xylem.csv` contains the isotopic signature of xylem water:

TODO...
