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

In case you're unfamiliar to Julia, there are various ways to run a script such as `main.jl`: One possibility is to open the Julia REPL and run the script using `include(“main.jl”)`. Alternatively, the editor VS Code in combination with the Julia extension ([julia-vscode.org](https://www.julia-vscode.org)), provides a complete IDE for programming in Julia.

!!! note

    TODO(bernhard): include units descriptions in tables below. And also discuss rounding.

### Needed time dependent parameters (daily time step): meteo data and stand properties

Time dependent parameters (climate and vegetation) are provided in the following form:

`BEA2016-reset-FALSE_meteoveg.csv`:

| dates      | globrad | tmax    | tmin    | vappres | windspeed | prec     | mesfl    | densef | height | lai   | sai   | age       |
| ---------- | ------- | ----    | -----   | ------- | --------- | ----     | -----    | ------ | ------ | ----- | ----  | -------   |
|            | weather | weather | weather | weather | weather   | weather  | stream   | stand  | stand  | stand | stand | stand     |
|            | (MJ/m2) | (°C)    | (°C)    | (kPa)   | (m/s)     | (mm/day) | (unused) | (-)    | (m)    | (-)   | (-)   | (years)   |
| 2016-01-01 | 4.08    | 5.5     | -1.1    | 0.53    | 1.22      | 3.2      | 0        | 1      | 23     | 1.752 | 1     | 200       |
| 2016-01-02 | 1.61    | 3.36    | -2.08   | 0.49    | 0.89      | 0.2      | 0        | 1      | 23     | 1.752 | 1     | 200.00274 |

Note that the second and third rows containing description and unit headers is not contained in the input dataset.

`precdat.csv``

- TODO(bernhard): support for `precdat.csv` is currently not implemented


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
"param_id","x"
"0_heat",0
"eslope",18.26
"aspect",225
"alb",0.2
"albsn",0.5
"c1",0.25
"c2",0.5
"c3",0.2
"wndrat",0.3
"fetch",5000
"z0w",0.005
"zw",2
"lwidth",0.1
"obsheight_x_czs",0.00325
"z0s",0.001
"lpc",4
"cs",0.035
"czs",0.13
"czr",0.05
...
"dtimax",0.5
"dswmax",0.05
"dpsimax",5e-04
```



### Needed soil properties
`soil_nodes.csv`:

| layer | midpoint | thick | mat  | psiini | rootden |
| ----- | -------- | ----- | ---- | ------ | ------- |
| (-)   | (m)      | (mm)  | (#)  | (kPa)  | (-)     |
| 1     | -0.02    | 40    | 1    | -6.3   | 0.029   |
| 2     | -0.06    | 40    | 2    | -6.3   | 0.025   |
| 3     | -0.14    | 120   | 3    | -6.3   | 0.020   |
| 4     | -0.325   | 250   | 4    | -6.3   | 0.012   |
| 5     | -0.6     | 300   | 5    | -6.3   | 0.005   |
| 6     | -0.925   | 350   | 6    | -6.3   | 0.002   |
| 7     | -1.15    | 100   | 7    | -6.3   | 0.001   |

Note that the second row containing units is not contained in the input dataset.

`soil_materials.csv`: (when using Mualem-van Genuchten parmetrization)

| mat  | ths   | thr   | alpha    | npar  | ksat      | tort  | gravel |
| ---- | ----- | ----- | -------- | ----- | --------- | ----- | ------ |
| (#)  | (-)   | (-)   | (m-1)    | (-)   | (mm/day)  | (-)   | (-)    |
| 1    | 0.715 | 0.069 | 1147.889 | 1.051 | 24864.678 | 4.670 | 0.010  |
| 2    | 0.669 | 0.069 | 1274.886 | 1.051 | 12881.486 | 4.478 | 0.175  |
| 3    | 0.656 | 0.069 | 1215.927 | 1.051 | 10516.617 | 4.502 | 0.175  |
| 4    | 0.602 | 0.069 | 795.204  | 1.051 | 3438.063  | 4.345 | 0.175  |
| 5    | 0.523 | 0.069 | 352.368  | 1.051 | 450.358   | 4.291 | 0.010  |
| 6    | 0.565 | 0.069 | 570.682  | 1.051 | 1488.210  | 5.148 | 0.010  |
| 7    | 0.467 | 0.069 | 164.564  | 1.051 | 67.978    | 0.010 | 0.950  |

Note that the second row containing units is not contained in the input dataset.

`soil_materials.csv`: (when using Clapp-Hornberger parametrization) (NOT IMPLEMENTED!)

| mat  | thsat | thetaf| psif    | bexp  | kf       | wtinf | gravel |
| ---- | ----- | ----- | ------- | ----- | -------  | ----- | ------ |
| (#)  | (-)   | (-)   | (kPa)   | (-)   | (mm/day) | (-)   | (-)    |
| NA   | NA    | NA    | NA      | NA    | NA       | NA    | NA     |



## Calibration data (calibration not yet implemented)
Roadmap to include calibration data intends to allow inclusion of:
- Throughfall amounts (for parametrisation of interception)
- Soil moisture (volumetric water content, θ)
- Soil matric potential (ψ)






!!! note

    USER GUIDE ENDS HERE.

    TODO(bernhard): after impelementation of function `read_inputData`
    Below is an illustration of the initial data sets that could be used.

#### Needed time dependent parameters (daily time step): meteo data and stand properties

| site_id | dates      | tmin   | tmax  | tmean  | prec | relhum | globrad | wind | vappres |
| ------- | ---------- | ------ | ----- | ------ | ---- | ------ | ------- | ---- | ------- |
|         |            | °C     | °C    | °C     | mm   | %      | MJ/m2   | m/s  | kPa     |
| BEA     | 01.01.2010 | -6.39  | 1.85  | -1.78  | 0.2  | 79.39  | 3.45    | 2.87 | 0.43    |
| BEA     | 02.01.2010 | -14.09 | -6.98 | -10.83 | 0    | 79.25  | 7.26    | 3.6  | 0.21    |

Note that the second row containing units is not contained in the input dataset.

| site     | X      | Y      | long_wgs84 | lat_wgs84   | tree height  | max root depth | LAI  | tree age  | % of deciduous plants |
| -------- | ------ | ------ | ---------- | ----------- | ------------ | -------------- | ---- | --------- | --------------------- |
| (-)      | (m)    | (-)    | (°E WGS84) | (°N WGS84)  | (m)          | (m)            | (-)  | (years)   | (%)                   |
| BEA      | 827262 | 165790 | 10.40555   | 46.60469346 | 3.5          | 1.1            | 1.91 | 80        | 45                    |
| BEA      | 827620 | 165710 | 10.41018   | 46.60385246 | 8.5          | 1.0            | 1.91 | 80        | 50                    |
TODO: remove (X, Y)

Note that the second row containing units is not contained in the input dataset.


#### Needed soil properties
LWF sites (TODO: format as LWFBrook90R):
- TODO(Bernhard): remove this once an example data set is here.

| site_id         | Total soil depth | horizon | texture | upper | lower | bd       | gravel   | sand     | silt     | clay     | c_org    |
| --------------- | ---------------- | ------- | ------- | ----- | ----- | -------- | -------- | -------- | -------- | -------- | -------- |
| unit            | m                | KA5     | GSTCS*  | (m)   | (m)   | (g cm-3) | fraction | (mass-%) | (mass-%) | (mass-%) | (mass-%) |
| Münstertal_pit1 | 1.2              | Ah      | Sl3     | 0     | -0.05 | 1.32     | 0.1      | 70       | 20       | 10       | 9        |
| -               | -                | Bv      | …       | -0.05 | -0.4  |          |          |          |          |          |          |
| Münstertal_pit2 | 1.3              | Ah      |         | ...   |       |          |          |          |          |          |          |
| -               | -                | Bv      |         | ...   |       |          |          |          |          |          |          |
\* GSTCS refers to: German soil texture classification system

Note that the second row containing units is not contained in the input dataset.

