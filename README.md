# LWFBrook90.jl
Implementation of the LWF-BROOK90 hydrological model in Julia

<p align="center">User info:
    <a href="https://www.repostatus.org/#wip"><img title="Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public." src="https://www.repostatus.org/badges/latest/wip.svg"/></a>
<a href="https://www.gnu.org/licenses/gpl-3.0"><img title="License: GPL v3" src="https://img.shields.io/badge/License-GPLv3-blue.svg"></a>
<a href="https://fabern.github.io/LWFBrook90.jl/stable"><img title="Stable" src="https://img.shields.io/badge/docs-stable-blue.svg"></a>
</p>
<p align="center">Developer info:
<a href="https://github.com/fabern/LWFBrook90.jl/actions"><img title="Build Status" src="https://github.com/fabern/LWFBrook90.jl/workflows/CI/badge.svg"/></a>
<a href="https://codecov.io/gh/fabern/LWFBrook90.jl"><img title="Coverage" src="https://codecov.io/gh/fabern/LWFBrook90.jl/branch/master/graph/badge.svg"/></a>
<a href="https://github.com/invenia/BlueStyle"><img title="Code Style: Blue" src="https://img.shields.io/badge/code%20style-blue-4495d1.svg"/></a>
<a href="https://fabern.github.io/LWFBrook90.jl/dev"><img title="Dev" src="https://img.shields.io/badge/docs-dev-blue.svg"></a>
<a href="https://discourse.julialang.org/t/good-practices-for-package-development-in-the-julia-ecosystem/8175/2">Git branch model</a>
</p>



## What is LWFBrook90.jl?
The model LWF-BROOK90 is a 1D Soil Vegetation Atmosphere Transport (SVAT) model, calculating the soil water balance in forest soil. Modelled processes include vertical soil water movement, soil and plant evapotranspiration and temporary storages in snowpack or interception layer.
LWFBrook90.jl is an implementation of the existing LWF-BROOK90 rewritten entirely in the Julia programming language.

Processes and state variables of the model in LWF-BROOK90 are summarised visually in Figure 1 below.
LWFBrook90.jl is developed with the following objectives in mind:
- [ ] support for stable isotopes (δ¹⁸O and δ²H) by including transport equation and fractionation processes
- [ ] efficient parameter estimation (optimizing computational costs)
- [ ] model flexibility for alternative processes parametrizations (possibly resulting in a [flexible model framework](https://presentations.copernicus.org/EGU2020/EGU2020-17975_presentation.pdf))



<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/LWFBrook90jl_overview_v1.4-juliacolors.png?raw=true" width="400"><p>
<p align="center">Figure 1: Summary of processes and state variables used in LWFBrook90.jl<p>

Further model description can be found in the [documentation](https://fabern.github.io/LWFBrook90.jl/stable).



## Acknowledgements
Development of LWFBrook90.jl builds on the following works:
- BROOK90 (v4.8) by C. Anthony Federer, Licensed under CC0 1.0, http://www.ecoshift.net/brook/brook90.htm
- LWFBrook90R (v0.4.3) by Paul Schmidt-Walter, Volodymyr Trotsiuk, Klaus Hammel, Martin Kennel, Anthony Federer, Robert Nuske, Licensed under GPL-3.0

Furthermore, Matthias Häni, Katrin Meusburger, Peter Waldner, Lorenz Walthert, Stephan Zimmermann of [WSL](http://www.wsl.ch) and the Long-term Forest Ecosystem Research (LWF) project of WSL are gratefully acknowledged for providing example data files located in `example/BEA2016*`.

For the license of LWFBrook90.jl see the file `LICENSE`.



## For Users: getting started
An example data set `BEA2016*` was generated using the R package LWFBrook90R and is located inthe folder `example/`. See outputs in the next section.

To run this example simulation simulation simply call LWFBrook90.run_example(). Note, that the first time it takes some time to load and compile the package. Another possibility is to follow the script `main.jl`, which sets up a simulation and shows some simple plotting commands.


Further documentation of [LWFBrook90.jl](https://fabern.github.io/LWFBrook90.jl/stable), and [BROOK90](http://www.ecoshift.net/brook/b90doc.html) are available.

### Example data set:
Following plots illustrate results of the provided data set. The scalar state variables and depth-depenedent (vector) state variables can be plotted:

<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/git-hash-b3f7183/2021-02-24_16h56_LWFBrook90Julia_plot_u_scalar.png?raw=true" width="400"><p>
<p align="center">Figure 2: Example simulation: scalar results<p>

<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/git-hash-b3f7183/2021-02-24_16h56_LWFBrook90Julia_plot_u_vector.png?raw=true" width="400"><p>
<p align="center">Figure 3: Example simulation: vector results soil water<p>

Tests are run to assert agreement with results from LWFBrook90R. Visualizations are reported below. Note that minor discrepancies are still present linked to the adaptive time stepping and intermediate updates of state variables.
<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/git-hash-b3f7183/2021-02-24_16h56_R-vs-Julia_comparison_DailyRawValues.png?raw=true" width="400"><p>
<p align="center">Figure 4: Comparing daily outputs of LWFBrook90R and LWFBrook90.jl for example data set over a year<p>

<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/git-hash-b3f7183/2021-02-24_16h56_R-vs-Julia_comparison_DailyRawValues_first3months.png?raw=true" width="400"><p>
<p align="center">Figure 5: Comparing daily outputs of LWFBrook90R and LWFBrook90.jl for example data set over 3 months<p>



## Limitations of LWFBrook90.jl
- Preprocessing steps to define a convenient format of input data sets are not yet designed. Currently, LWFBrook90.jl makes use of the preprocessing steps provded in the the R package LWFBrook90R.
- LWFBrook90R is based on LWF-BROOK90. LWF-BROOK90 was itself based an older version of BROOK90 (v3.1F). Developments in BROOK90 up to v4.8 (such as intercell averages of hydraulic conductivity `KKMEAN`) are therefore only partially included, but not activated by default in LWFBrook90.jl.
- Currently a part of LWFBrook90R Reset==1, is not implemented in LWFBrook90.jl. (A test implementation of this is available in the branch `005-simplify-time-step-control`.)

### Improve agreement with LWFBrook90R
Note that some features of LWFBrook90R are not implemented in the main version of LWFBrook90.jl. The time step adaptivity and Reset==1 are major ones that require some code refactoring that is not how the library for ODEs DiffEq.jl is intended to be used. Because of that implementation of these features is currently in a feature branch `005` here on git. Below are some of the results of that branch:

<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/git-hash-2-55ca42d-feature005/2021-02-24_19h10_R-vs-Julia_comparison_DailyRawValues.png?raw=true" width="400"><p>
<p align="center">Figure 6: Comparing daily outputs of LWFBrook90R and experimental LWFBrook90.jl:feature-005 for example data set over a year<p>

<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/git-hash-2-55ca42d-feature005/2021-02-24_19h10_R-vs-Julia_comparison_DailyRawValues_first3months.png?raw=true" width="400"><p>
<p align="center">Figure 7: Comparing daily outputs of LWFBrook90R and experimental LWFBrook90.jl:feature-005 for example data set over 3 months<p>



## For Developers:
Any help in form of discussions, pull requests, example data sets, or otherwise is very welcome. Please don't hesitate to contact us.

LWFBrook90.jl makes use of [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) to solve the system of Ordinary Differential Equations (ODE). Each state variable (`u`) has a corresponding ODE. The ODEs are defined by their right hand side defined in the function `f` (which sets `du` that is the change in `u`). The right hand side `f(u,p,t)` depends on time `t`, parameters `p` (time-dependent or constant), and the state `u`.

Variable naming generally follows the convention by BROOK90 and LWFBrook90R, but uses additionally a prefix to indicate their dependencies (`p_*`, `p_fT_*`, `p_fu_*` for constant, time dependent or state dependent parameters, `u_*`, `u_aux_*` for elementary and auxiliary state variables respectively, and `aux_du_*` for auxiliary rate of changes of state variables).

Note that some state variables (rain and snow interception storage, `u_INTR`, `u_INTS` and snow storage and energy `u_SNOW`, `u_CC`, `u_SNOWLQ`) are updated once per simulation day and other state variables (groundwater and soil water storages `u_GWAT`,`u_SWATI`) are solved on a higher resolved time discretization set by the ODE solver, resulting in a scheme based on operator splitting.