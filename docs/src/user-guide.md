# User Guide



## Installation and usage

### Installation
First, download and install Julia: [https://julialang.org/downloads/](https://julialang.org/downloads/).

To install LWFBrook90.jl open a Julia REPL and run:
```Julia
using Pkg
Pkg.add("LWFBrook90.jl")
```
!!! note

    Currently the package is not yet registerd. So above command will not work.

To install a current development version of the package you can install it directly from GitHub using:
```Julia
using Pkg
Pkg.add(PackageSpec(name = "Omniscape", rev = "main"))
```

Dependencies of LWFBrookJulia.jl should automatically be installed.

### Usage
Check out a step by step in guide for installation in [Example](@ref)

The steps in a typical simulation script are:
- load the package `using LWFBrook90`
- load the package dependency `using DifferentialEquations`
- read input data
- set up model options
- set up an `ODE` problem (`u0`,`tspan`, `p`) and solve it with DifferentialEquations.jl
- plot and/or postprocess output



## Input data
To run a simulation following input data are needed
- siteparam:
- meteo:
- precdat (currently unused):
- pdur: (currently unused??)
- soil_materials:
- soil_nodes:
- param:

Which contain problem-specific settings such as site paramters, meteorological drivers, and soil parameters as well as implementation-specific control settings in `input_param`.

The input files can be specified in two ways:
1. `read_inputData` (TODO: not yet implemented) or
2. [`read_LWFBrook90R_inputData`](@ref) (implemented).

The first input data is not yet implemented.
The second type of input data is more detailed and should be generated automatically by the user. The structure of the second type of input data is illustrated by the example input data `BEA2016-*` set located in the folder `example/`.

An easy way to generate the second type of input data is by setting up a simulation with the R package [LWFBrook90R (v0.4.3)](https://github.com/pschmidtwalter/LWFBrook90R#usage). Instead of running the simulation with `run_LWFB90()`, the same arguments can be used to generate the input files for LWFBrook90.jl using the R function provided in the file `generate_LWFBrook90jl_Input.R`.

To load load input data and prepare a simulation follow the instructions in section [Example](@ref) or alternatively use the sample script `main.jl`.

!!! note

    TODO(bernhard): include units descriptions in tables below. And also discuss rounding.

### Needed time dependent parameters (daily time step): meteo data and stand properties

Time dependent parameters (climate and vegetation) are provided in the following form:

`BEA2016-reset-FALSE_climveg.csv`:

| dates      | globrad | tmax | tmin  | vappres | windspeed | prec | mesfl | densef | height | lai   | sai  | age       |
| ---------- | ------- | ---- | ----- | ------- | --------- | ---- | ----- | ------ | ------ | ----- | ---- | --------- |
| 2016-01-01 | 4.08    | 5.5  | -1.1  | 0.53    | 1.22      | 3.2  | 0     | 1      | 23     | 1.752 | 1    | 200       |
| 2016-01-02 | 1.61    | 3.36 | -2.08 | 0.49    | 0.89      | 0.2  | 0     | 1      | 23     | 1.752 | 1    | 200.00274 |



### Needed constant meteo data
`pdur.csv`:
```
"x"
4
4
4
4
4
4
4
4
4
4
4
4
```

`param.csv`:

```
"x","param_id"
366,"ndays"
0,"0_heat"
18.26,"eslope"
225,"aspect"
0.2,"alb"
0.5,"albsn"
0.25,"c1"
0.5,"c2"
0.2,"c3"
0.3,"wndrat"
5000,"fetch"
0.005,"z0w"
2,"zw"
0.1,"lwidth"
0.00325,"obsheight_x_czs"
0.001,"z0s"
4,"lpc"
0.035,"cs"
0.13,"czs"
0.05,"czr"
...
0.5,"dtimax"
0.05,"dswmax"
5e-04,"dpsimax"
```

`precdat.csv``

- TODO(bernhard): support for `precdat.csv` is currently not implemented


### Needed soil properties
`soil_nodes.csv`:

| layer | midpoint | thick | mat  | psiini | rootden |
| ----- | -------- | ----- | ---- | ------ | ------- |
| 1     | -0.02    | 40    | 1    | -6.3   | 0.029   |
| 2     | -0.06    | 40    | 2    | -6.3   | 0.025   |
| 3     | -0.14    | 120   | 3    | -6.3   | 0.020   |
| 4     | -0.325   | 250   | 4    | -6.3   | 0.012   |
| 5     | -0.6     | 300   | 5    | -6.3   | 0.005   |
| 6     | -0.925   | 350   | 6    | -6.3   | 0.002   |
| 7     | -1.15    | 100   | 7    | -6.3   | 0.001   |

`soil_materials.csv`:

| mat  | ths   | thr   | alpha    | npar  | ksat      | tort  | gravel |
| ---- | ----- | ----- | -------- | ----- | --------- | ----- | ------ |
| 1    | 0.715 | 0.069 | 1147.889 | 1.051 | 24864.678 | 4.670 | 0.010  |
| 2    | 0.669 | 0.069 | 1274.886 | 1.051 | 12881.486 | 4.478 | 0.175  |
| 3    | 0.656 | 0.069 | 1215.927 | 1.051 | 10516.617 | 4.502 | 0.175  |
| 4    | 0.602 | 0.069 | 795.204  | 1.051 | 3438.063  | 4.345 | 0.175  |
| 5    | 0.523 | 0.069 | 352.368  | 1.051 | 450.358   | 4.291 | 0.010  |
| 6    | 0.565 | 0.069 | 570.682  | 1.051 | 1488.210  | 5.148 | 0.010  |
| 7    | 0.467 | 0.069 | 164.564  | 1.051 | 67.978    | 0.010 | 0.950  |



!!! note

    TODO(bernhard): clean up documentation from here on...


## Calibration data (calibration not yet implemented)
Possible calibration data could in the future could be:
- Throughfall amounts (for parametrisation of interception)
- Soil moisture (θ)
- Soil matric potential (ψ)



!!! note

    TODO(bernhard): after impelementation of function `read_inputData`
    Below is an illustration of the initial data sets that could be used.

#### Needed time dependent parameters (daily time step): meteo data and stand properties

| site_id | dates      | tmin   | tmax  | tmean  | prec | relhum | globrad | wind | vappres |
| ------- | ---------- | ------ | ----- | ------ | ---- | ------ | ------- | ---- | ------- |
|         |            | °C     | °C    | °C     | mm   | %      | MJ/m2   | m/s  | kPa     |
| BEA     | 01.01.2010 | -6.39  | 1.85  | -1.78  | 0.2  | 79.39  | 3.45    | 2.87 | 0.43    |
| BEA     | 02.01.2010 | -14.09 | -6.98 | -10.83 | 0    | 79.25  | 7.26    | 3.6  | 0.21    |


| site     | X      | Y      | long_wgs84 | lat_wgs84   | tree height (m) | max root depth (m) | LAI  | tree age (yr) | % of deciduous plants |
| -------- | ------ | ------ | ---------- | ----------- | --------------- | ------------------ | ---- | ------------- | --------------------- |
| BEA      | 827262 | 165790 | 10.40555   | 46.60469346 | 3.5             | 1.1                | 1.91 | 80            | 45                    |
| BEA      | 827620 | 165710 | 10.41018   | 46.60385246 | 8.5             | 1.0                | 1.91 | 80            | 50                    |



#### Needed soil properties
LWF sites (TODO: format as LWFBrook90R):
- TODO(Bernhard): remove this once an example data set is here.

| site_id             | Total  soil depth | horizon | texture                                   | upper | lower | bd       | gravel   | sand     | silt     | clay     | c_org    |
| ------------------- | ---------------- | ------- | ----------------------------------------- | ----- | ----- | -------- | -------- | -------- | -------- | -------- | -------- |
| unit                | m                | KA5     | GSTCS* | (m)   | (m)   | (g cm-3) | fraction | (mass-%) | (mass-%) | (mass-%) | (mass-%) |
| Münstertal_pit1 | 1.2              | Ah      | Sl3                                       | 0     | -0.05 | 1.32     | 0.1      | 70       | 20       | 10       | 9        |
| -                   | -                | Bv      | …                                         | -0.05 | -0.4  |          |          |          |          |          |          |
| Münstertal_pit2 | 1.3              | Ah      |                                           | ...   |       |          |          |          |          |          |          |
| -                   | -                | Bv      |                                           | ...   |       |          |          |          |          |          |          |
\* GSTCS refers to: German soil texture classification system


