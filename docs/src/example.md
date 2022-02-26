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
# Read in input data
input_prefix = "BEA2016-reset-FALSE"
input_path = "example/"*input_prefix*"-input/"

####################
(input_meteoveg,
    input_meteoveg_reference_date,
    input_param,
    input_storm_durations,
    input_initial_conditions,
    input_soil_horizons,
    simOption_FLAG_MualVanGen) = read_inputData(input_path, input_prefix)

input_soil_discretization = discretize_soil(input_path, input_prefix)
####################
```

Then simulation parameters are defined. The user has the possibility to modify the input variables before continuing with the simulation. Further the user can select whether intermediate quantities such as e.g. evaporation fluxes should be stored during simulation (`compute_intermediate_quantities = true`) or whether processing steps during simulation should be kept to a minimum for performance reasons (`compute_intermediate_quantities = false`).

```Julia
####################
# Define solver options
Reset = false                          # currently only Reset = 0 implemented
compute_intermediate_quantities = true # Flag whether ODE containes additional quantities than only states

# Override input file settings
# Here possibility to check and override dataframes input_[...] manually
    # # E.g:
    # # Soil hydraulic model
    # input_param[1,"NOOUTF"] = true # `true` if outflow from roots prevented, `false` if allowed
####################

```

Then functions from the package are used to define the problem that will be handed to the solver from DifferentialEquations.jl. That is we need to define `p`, `u0`, and `tspan`

`p`:
```Julia
####################
# Define parameters for differential equation
p = define_LWFB90_p(
    input_meteoveg,
    input_meteoveg_reference_date,
    input_param,
    input_storm_durations,
    input_soil_horizons,
    input_soil_discretization,
    simOption_FLAG_MualVanGen;
    Reset = Reset,
    compute_intermediate_quantities = compute_intermediate_quantities)
####################
```
`u0`:
```Julia
####################
# Define initial states of differential equation
# state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI

######
# Transform initial value of auxiliary state u_aux_PSIM_init into state u_SWATIinit:
u_aux_PSIM_init = input_soil_discretization[:,"uAux_PSIM_init_kPa"]
if any( u_aux_PSIM_init.> 0)
    error("Initial matrix psi must be negative or zero")
end
######

# Create u0 for DiffEq.jl
u0 = define_LWFB90_u0(p, input_initial_conditions,
                      uSoil_initial,
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
# tspan = (0.,  5.) # simulate days 0 to 5 (in the reference frame of the input data)
# Simulation specific period (provided it is within the input data):
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
    saveat = tspan[1]:tspan[2], dt=1e-6, adaptive = true); # dt is initial dt, but adaptive
####################
```


Note to use a progress bar indicating advancement of the solving is possible. It is sufficient to load the package `using ProgressLogging`.

The generated solution can be plotted using the plotting recipes of DifferntialEquations.jl for Plots.jl ([see instructions](https://diffeq.sciml.ai/stable/basics/plot/)). An example is provided below:

```Julia
####################
## Plotting
using Plots
sol_LWFBrook90_Dates =
    LWFBrook90.jl.RelativeDaysFloat2DateTime.(
        sol_LWFBrook90.t,
        input_meteoveg_reference_date)

# Plot scalar quantities
# Using dates (but not interpolated)
plot(sol_LWFBrook90_Dates,
    sol_LWFBrook90[[1,2,3,4,5,6],:]',
    label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

# Using simple plot recipe that interpolates, but without dates
plot(sol_LWFBrook90;
    vars = [1, 2, 3, 4, 5, 6],
    label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

# Plot vector quantities
# http://docs.juliaplots.org/latest/generated/gr/#gr-ref43
x = sol_LWFBrook90_Dates
y = cumsum(pfile_soil["THICK"])
z = sol_LWFBrook90[7 .+ (0:example["NLAYER"]-1), :]./pfile_soil["THICK"]
heatmap(x, y, z,
    yflip = true,
    xlabel = "Date",
    ylabel = "Depth [mm]",
    colorbar_title = "θ [-]")
```

## Plotting results
Following plots illustrate results of the provided data set. The scalar state variables and depth-depenedent (vector) state variables can be plotted:
```@raw html
<p align="center">
<img src="../assets/example-results/2021-03-16_11h27-Reset0-git+6e6946d+clean_LWFBrook90Julia_plot_u_scalar.png" width="400"><br>
<br><em><b>Figure 2</b>: Example simulation: scalar results</em><br>
<p>
```
```@raw html
<p align="center">
<img src="../assets/example-results/2021-03-16_11h27-Reset0-git+6e6946d+clean_LWFBrook90Julia_plot_u_vector.png" width="400"><br>
<br><em><b>Figure 3</b>: Example simulation: vector results soil water</em><br>
<p>
```


## Comparison with LWFBrook90R
Tests are run to assert agreement with results from LWFBrook90R. Visualizations are reported below. Note that minor discrepancies ```@raw htmlare still present linked to the adaptive time stepping and intermediate updates of state variables.
```@raw html
<p align="center">
<img src="../assets/example-results/2021-03-16_11h27-Reset0-git+6e6946d+clean_R-vs-Julia_Daily_first12M.png" width="400"><br>
<br><em><b>Figure 4</b>: Comparing daily outputs of LWFBrook90R and LWFBrook90.jl for example data set over a year</em><br>
<p>
```
```@raw html
<p align="center">
<img src="../assets/example-results/2021-03-16_11h27-Reset0-git+6e6946d+clean_R-vs-Julia_Daily_first2M.png" width="400"><br>
<br><em><b>Figure 5</b>: Comparing daily outputs of LWFBrook90R and LWFBrook90.jl for example data set over 2 months</em><br>
<p>
```

Note that some features of LWFBrook90R are not implemented in the main version of LWFBrook90.jl. The time step adaptivity and `Reset==1` are major ones that require some code refactoring that is not how the library for ODEs DiffEq.jl is intended to be used. Because of that implementation of these features is currently in a feature branch here on git `feature 005`.