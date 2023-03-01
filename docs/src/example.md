# Example manual

Example data from Beatenberg is located in subfolder `examples/`. WSL is acknowledged for providing the input data (see section [Acknowledgments](@ref)).

## Step-by-step instructions
To run your own simulations: follow the step-by-step instructions below.

Load packages:
```Julia
using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5
```

 Define and read in input data
 ```Julia
# Read in input data
input_prefix = "isoBEA2010-18-reset-FALSE"
input_path = "examples/"*input_prefix*"-input/"

####################
# Define simulation model by reading in system definition and input data
model = SPAC(input_path, input_prefix;
                simulate_isotopes = contains(input_prefix, "iso"));
####################


####################
# Prepare simulation by discretizing spatial domain
simulation = LWFBrook90.discretize(model; tspan = (0,100));

# Solve system of ODEs:
LWFBrook90.simulate!(simulation)

sol_LWFBrook90 = simulation.ODESolution
####################
```

TODO: rewrite this documentation: Then simulation parameters are defined. The user has the possibility to modify the input variables before continuing with the simulation. Further the user can select whether intermediate quantities such as e.g. evaporation fluxes should be stored during simulation (`compute_intermediate_quantities = true`) or whether processing steps during simulation should be kept to a minimum for performance reasons (`compute_intermediate_quantities = false`). TODO: define tspan... define Δz as additional input to SPAC()
TODO: check if still true: Note to use a progress bar indicating advancement of the solving is possible. It is sufficient to load the package `using ProgressLogging`.

The generated solution can be plotted using the in-built function `LWFBrook90.ISO.plotisotopes()`. Alternatively plots can be constructed using the plotting recipes of DifferentialEquations.jl for Plots.jl ([see instructions](https://diffeq.sciml.ai/stable/basics/plot/)). Examples of both approaches are provided below:

```Julia
####################
## Plotting
using Plots, Measures

###
# A) Use in-built plotting function
optim_ticks = (x1, x2) -> Plots.optimize_ticks(x1, x2; k_min = 4)
pl_inbuilt = LWFBrook90.ISO.plotisotopes(
    sol_LWFBrook90, optim_ticks;
    layout = grid(4, 1, heights=[0.1 ,0.4, 0.1, 0.4]),
    size=(1000,1400), dpi = 300, leftmargin = 15mm);
savefig(pl_inbuilt, "Isotopeplots_pl_inbuilt.png")

###
# B) Construct plots yourself using the solution object
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