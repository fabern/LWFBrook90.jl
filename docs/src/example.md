# Example

Example data from Beatenberg is located in subfolder `example/`. WSL is acknowledged for providing the input data (see section [Acknowledgments](@ref)).

## Step-by-step instructions
To run the example simulation simulation simply call LWFBrook90.run_example(). For more control either run the script `main.jl`or follow the step-by-step instructiosn below.

Load packages:
```Julia
using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5
```

 Define and read in input data
 ```Julia
 # 1a) Read in input data
input_prefix = "BEA2016-reset-FALSE"
input_path = "example/"*input_prefix*"-input/"

(input_meteoveg,
    input_param,
    input_siteparam,
    input_precdat,    #TODO(bernhard): input_precdat is unused
    input_pdur,
    input_soil_materials,
    input_soil_nodes,
    input_reference_date) = read_LWFBrook90R_inputData(input_path, input_prefix)
```

Now the user has the possibility to modify these input objects:
```Julia
# 1b) Here posibility to modify dataframes input_[...] manually
```

After that the objects are further parsed and simulation parameters are defined. The user can select whether intermediate quantities such as e.g. evaporation fluxes should be stored during simulation (`compute_intermediate_quantities = true`) or whether processing steps during simulation should be kept to a minimum for performance reasons (`compute_intermediate_quantities = false`).

```Julia
# 1c) Parse loaded/redefined input files
(pfile_meteoveg, pfile_param, pfile_siteparam, pfile_precdat, pfile_pdur, pfile_soil) =
    derive_params_from_inputData(input_meteoveg,
                                 input_param,
                                 input_siteparam,
                                 input_precdat,
                                 input_pdur,
                                 input_soil_materials,
                                 input_soil_nodes,
                                 input_reference_date)

```
```Julia
####################
# Define simulation
# Soil hydraulic model
IMODEL = pfile_param["IMODEL"] # 0 for Clapp-Hornberger; 1 for Mualem-van Genuchten
NLAYER = pfile_param["NLAYER"]

# Define solver options
NOOUTF    = 1 == pfile_param["NOOUTF"] # 1 if no outflow allowed from roots, otherwise 0
Reset     = 0 # currently only Reset = 0 implemented

constant_dt_solver = 1 # [days]

# Flag whether ODE containes additional quantities than only states
compute_intermediate_quantities = true
####################
```


Then functions from the package are used to define the problem that will be handed to the solver from DifferentialEquations.jl. That is we need to define `p`, `u0`, and `tspan`

`p`:
```Julia
####################
# Define parameters for differential equation
p = define_LWFB90_p(NLAYER, IMODEL, constant_dt_solver,
                    NOOUTF, Reset, compute_intermediate_quantities,
                    pfile_meteoveg,
                    pfile_siteparam,
                    pfile_param,
                    pfile_soil,
                    pfile_pdur)
####################
```
`u0`:
```Julia
####################
# Define initial states of differential equation
# state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI
u_GWAT_init = pfile_siteparam["u_GWAT_init"]
u_SNOW_init = pfile_siteparam["u_SNOW_init"]
u_INTS_init = pfile_param["INTS_init"];
u_INTR_init = pfile_param["INTR_init"];
u_CC_init     = 0; # any initial snow has zero liquid water and cold content
u_SNOWLQ_init = 0; # any initial snow has zero liquid water and cold content

u_aux_PSIM_init = pfile_soil["PSIM_init"]

######
# Transform initial value of auxiliary state u_aux_PSIM_init into state u_SWATIinit:
if any( u_aux_PSIM_init.> 0)
    error("Initial matrix psi must be negative or zero")
end
p_soil = p[1][1][6] # TODO(bernhard): this hardcoded index is dangerous in case definition of p vector changes

u_aux_WETNESinit = LWFBrook90.KPT.FWETNES(u_aux_PSIM_init, p_soil)
u_SWATIinit      = p_soil.p_SWATMX ./ p_soil.p_THSAT .* LWFBrook90.KPT.FTheta(u_aux_WETNESinit, p_soil)
######


# Create u0 for DiffEq.jl
u0 = define_LWFB90_u0(u_GWAT_init,
                      u_INTS_init,
                      u_INTR_init,
                      u_SNOW_init,
                      u_CC_init,
                      u_SNOWLQ_init,
                      u_SWATIinit,
                      compute_intermediate_quantities)
####################
```

`tspan`:
```Julia
####################
# Define ODE problem which consists of
#   - definition of right-hand-side (RHS) function f
#   - definition of callback function cb
#   - u0:     initial condition of states
#   - tspan:  definition of simulation time span
#   - p:      parameters

# Define simulation time span:
tspan = (minimum(input_meteoveg[:,"days"]),
         maximum(input_meteoveg[:,"days"])) # simulate all available days
# Define ODE:
ode_LWFBrook90 = define_LWFB90_ODE(u0, tspan, p)

# Alternative definitions of tspan:
# tspan = (0.,  5.) # simulate 5 days
# tspan = (0.,  100.) # simulate 100 days
# simulates specific period:
# tspan = (LWFBrook90.jl.DateTime2RelativeDaysFloat(DateTime(1980,1,1), reference_date),
#          LWFBrook90.jl.DateTime2RelativeDaysFloat(DateTime(1985,1,1), reference_date))
####################
```


Then the ODE problem can be solved:
```Julia
####################
## Solve ODE:
sol_LWFBrook90 = solve(ode_LWFBrook90, Tsit5();
    progress = true,
    saveat = tspan[1]:tspan[2], dt=1e-6, adaptive = true); # dt will be overwritten
####################
```


Note to use the progress bar during the solve the optional package should additionally be installed and loaded: `using ProgressLogging`.



The generated solution can be plotted using the plotting recipes of DifferntialEquations.jl for Plots.jl ([see instructions](https://diffeq.sciml.ai/stable/basics/plot/))

```Julia
####################
## Plotting
using Plots
sol_LWFBrook90_Dates =
    LWFBrook90.jl.RelativeDaysFloat2DateTime.(
        sol_LWFBrook90.t,
        input_reference_date)

# Plot scalar quantities
# Using dates (but not interpolated)
plot(sol_LWFBrook90_Dates,
    sol_LWFBrook90[[1,2,3,4,5,6],:]',
    label=["GWAT" "INTS" "INTR" "SNOW" "CC" "SNOWLQ"])

# Using simple plot recipe that interpolates, but without dates
plot(sol_LWFBrook90;
    vars = [1, 2, 3, 4, 5, 6],
    label=["GWAT" "INTS" "INTR" "SNOW" "CC" "SNOWLQ"])

# Plot vector quantities
# http://docs.juliaplots.org/latest/generated/gr/#gr-ref43
x = sol_LWFBrook90_Dates
y = cumsum(pfile_soil["THICK"])
z = sol_LWFBrook90[7 .+ (0:example["NLAYER"]-1), :]./pfile_soil["THICK"]
heatmap(x, y, z,
    yflip = true,
    xlabel = "Date",
    ylabel = "Depth",
    colorbar_title = "θ")
```

## Plotting results
Following plots illustrate results of the provided data set. The scalar state variables and depth-depenedent (vector) state variables can be plotted:
```@raw html
<p align="center">
<img src="../assets/git-hash-b3f7183/2021-02-24_16h56_LWFBrook90Julia_plot_u_scalar.png" width="400"><br>
<br><em><b>Figure 2</b>: Example simulation: scalar results</em><br>
<p>
```
```@raw html
<p align="center">
<img src="../assets/git-hash-b3f7183/2021-02-24_16h56_LWFBrook90Julia_plot_u_vector.png" width="400"><br>
<br><em><b>Figure 3</b>: Example simulation: vector results soil water</em><br>
<p>
```


## Comparison with LWFBrook90R
Tests are run to assert agreement with results from LWFBrook90R. Visualizations are reported below. Note that minor discrepancies ```@raw htmlare still present linked to the adaptive time stepping and intermediate updates of state variables.
```@raw html
<p align="center">
<img src="../assets/git-hash-b3f7183/2021-02-24_16h56_R-vs-Julia_comparison_DailyRawValues.png" width="400"><br>
<br><em><b>Figure 4</b>: Comparing daily outputs of LWFBrook90R and LWFBrook90.jl for example data set over a year</em><br>
<p>
```
```@raw html
<p align="center">
<img src="../assets/git-hash-b3f7183/2021-02-24_16h56_R-vs-Julia_comparison_DailyRawValues_first3months.png" width="400"><br>
<br><em><b>Figure 5</b>: Comparing daily outputs of LWFBrook90R and LWFBrook90.jl for example data set over 3 months</em><br>
<p>
```

Note that some features of LWFBrook90R are not implemented in the main version of LWFBrook90.jl. The time step adaptivity and `Reset==1` are major ones that require some code refactoring that is not how the library for ODEs DiffEq.jl is intended to be used. Because of that implementation of these features is currently in a feature branch here on git `feature 005`. Below are some of the results of that code:

```@raw html
<p align="center">
<img src="../assets/git-hash-2-55ca42d-feature005/2021-02-24_19h10_R-vs-Julia_comparison_DailyRawValues.png" width="400"><br>
<br><em><b>Figure 6</b>: Comparing daily outputs of LWFBrook90R and experimental LWFBrook90.jl:feature-005 for example data set over a year</em><br>
<p>
```

```@raw html
<p align="center">
<img src="../assets/git-hash-2-55ca42d-feature005/2021-02-24_19h10_R-vs-Julia_comparison_DailyRawValues_first3months.png" width="400"><br>
<br><em><b>Figure 7</b>: Comparing daily outputs of LWFBrook90R and experimental LWFBrook90.jl:feature-005 for example data set over 3 months</em><br>
<p>
```