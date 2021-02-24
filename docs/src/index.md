```@meta
CurrentModule = LWFBrook90Julia
```

# LWFBrook90Julia

LWFBrook90Julia is a 1D soil vegetation atmosphere transport model for forested areas based on LWFBrook90R.

- TODO(bernhard): describe use of model here.
- TODO(bernhard): Remove paragraph on implementations: "LWFBrook90R is implemented in Fortran and has an R package as interface. LWFBrook90Julia is 100% implemented in Julia. It defines the dynamical system in terms of ordinary differential equations and corresponding (time-varying) parameters and makes use of the package DifferentialEquations.jl for solving this dynamical system for a specified time `tspan` and intial conditions `u0`."
- TODO(bernhard): Refer the documentation of original Brook90: http://www.ecoshift.net/brook/b90doc.html

## Input data
To run a simulation following input data are needed
- siteparam
- meteo
- precdat (currently unused)
- pdur
- soil_materials
- soil_nodes
- param

Which contain problem-specific settings such as site paramters, meteorological drivers, and soil parameters as well as implementation-specific control settings in `input_param`.

These files need to be provided as CSVs in a single folder and can be loaded by `read_inputData(folder)`.

The structure of these input data is illustrated by the example input data set located in "".
They are inspired by the definition of the input arguments to the function `run_LWFB90()` in the R package `LWFBrook90R`.
[They are acutally equal with the arguments that are handed internally to the function `r_lwfbrook90()`]

### Needed time dependent parameters (daily time step): meteo data and stand properties

Time dependent parameters (climate and vegetation) are provided in the following form:

KAUFENRING (LWFBrook90R):
- TODO(Bernhard): remove this once an example data set is prepared.
- TODO(bernhard): include site_id
- TODO(bernhard): include tmean
- TODO(bernhard): remove MESFL (measured streamflow)
- TODO(bernhard): switch to dates (YYYY.MM.DD) instead of 3 columns
- TODO(bernhard): first row (ID, ...) and third row (UNITS) are not part of csv

```
climveg.csv
```

| yr   | mo | da | globrad | tmax  | tmin   | vappres | wind | prec | mesfl | densef | height | lai | sai   | age |
| ---- | ---| ---| ------- | ----- | ------ | ------- | ---- | ---- | ----- | ------ | ------ | --- | ----- | --- |
| ID   | ID | ID | METEO | METEO  | METEO   | METEO   | METEO| METEO| STREAM| VEG    | VEG    | VEG | VEG   | VEG |
|      |    |    | MJ/m2   | °C    | °C     | kPa     | m/s  | mm   |       |        |        |     |       |     |
| 1980 | 1  | 1  | 3.45    | 1.85  | -6.39  | 0.5     | 2.87 | 0.2  | 0     | 1      | 25     | 4.8 | 0.875 | 100 |



LWF sites (TODO: format as LWFBrook90R):
- TODO(Bernhard): remove this once an example data set is here.

| site_id | dates      | tmin   | tmax  | tmean  | prec | relhum | globrad | wind | vappres |
| ------- | ---------- | ------ | ----- | ------ | ---- | ------ | ------- | ---- | ------- |
|         |            | °C     | °C    | °C     | mm   | %      | MJ/m2   | m/s  | kPa     |
| BEA     | 01.01.2010 | -6.39  | 1.85  | -1.78  | 0.2  | 79.39  | 3.45    | 2.87 | 0.43    |
| BEA     | 02.01.2010 | -14.09 | -6.98 | -10.83 | 0    | 79.25  | 7.26    | 3.6  | 0.21    |


| site     | X      | Y      | long_wgs84 | lat_wgs84   | tree height (m) | max root depth (m) | LAI  | tree age (yr) | % of deciduous plants |
| -------- | ------ | ------ | ---------- | ----------- | --------------- | ------------------ | ---- | ------------- | --------------------- |
| BEA      | 827262 | 165790 | 10.40555   | 46.60469346 | 3.5             | 1.1                | 1.91 | 80            | 45                    |
| BEA      | 827620 | 165710 | 10.41018   | 46.60385246 | 8.5             | 1.0                | 1.91 | 80            | 50                    |



### Needed constant meteo data

```
pdur.csv
```

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

`param.csv``

```

"x"
11323
0
0
0
0.14
0.14
0.25
0.5
0.2
0.3
5000
0.005
2
0.004
0.00325
0.001
4
0.035
0.13
...

```

`precdat.csv``

- TODO(bernhard): currently not implemented





### Needed soil properties

KAUFENRING:
- TODO(Bernhard): remove this once an example data set is here.

`soil_materials.csv`:

| mat | ths  | thr   | alpha | npar | ksat  | tort   | gravel |
| --- | ---  | ---   | ---   | ---  | ---   | ---    | ---    |
| 1   | 0.44 | 0.069 | 2.6   | 1.18 | 196.0 | 1.03   | 0.016  |
| 2   | 0.39 | 0.069 | 3.86  | 1.14 | 30.6  | -0.87  | 0.017  |


`soil_nodes.csv`:

| layer | midpoint | thick | mat | psiini | rootden|
| ---   | ----     | ----- | --- | ------ | ------ |
| 1     | -0.01    | 20    | 1   | -6.3   | 0.024  |
| 2     | -0.03    | 20    | 1   | -6.3   | 0.023  |
| 3     | -0.05    | 20    | 1   | -6.3   | 0.022  |
| 4     | -0.07    | 20    | 1   | -6.3   | 0.020  |
| 5     | -0.09    | 20    | 1   | -6.3   | 0.020  |
| 6     | -0.125   | 50    | 1   | -6.3   | 0.018  |
| 7     | -0.175   | 50    | 1   | -6.3   | 0.016  |
| 8     | -0.225   | 50    | 1   | -6.3   | 0.014  |
| 9     | -0.275   | 50    | 1   | -6.3   | 0.012  |
| 10    | -0.35    | 100   | 2   | -6.3   | 0.010  |

LWF sites (TODO: format as LWFBrook90R):
- TODO(Bernhard): remove this once an example data set is here.

| site_id             | Total soilde pth | horizon | texture                                   | upper | lower | bd       | gravel   | sand     | silt     | clay     | c_org    |
| ------------------- | ---------------- | ------- | ----------------------------------------- | ----- | ----- | -------- | -------- | -------- | -------- | -------- | -------- |
| unit                | m                | KA5     | German soil texture classification system | (m)   | (m)   | (g cm-3) | fraction | (mass-%) | (mass-%) | (mass-%) | (mass-%) |
| ***\*Münstertal_pit1\**** | 1.2              | Ah      | Sl3                                       | 0     | -0.05 | 1.32     | 0.1      | 70       | 20       | 10       | 9        |
| -                   | -                | Bv      | …                                         | -0.05 | -0.4  |          |          |          |          |          |          |
| ***\*Münstertal_pit2\**** | 1.3              | Ah      |                                           | ...   |       |          |          |          |          |          |          |
| -                   | -                | Bv      |                                           | ...   |       |          |          |          |          |          |          |



## Calibration data (calibration not yet implemented)
Possible calibration data could in the future could be:
- Throughfall amounts (for parametrisation of interception)
- Soil moisture (θ)
- Soil matric potential (ψ)



# Package content

```@index
```
