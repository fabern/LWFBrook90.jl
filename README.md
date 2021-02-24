# LWFBrook90Julia
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.com/fabern/LWFBrook90Julia.jl.svg?token=Wmy6jUbNaUsJTRx8zJVf&branch=main)](https://travis-ci.com/fabern/LWFBrook90Julia.jl)
[![Coverage](https://codecov.io/gh/fabern/LWFBrook90Julia.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/fabern/LWFBrook90Julia.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fabern.github.io/LWFBrook90Julia.jl/stable)
<!--- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fabern.github.io/LWFBrook90Julia.jl/dev)-->

## What is LWFBrook90Julia?
Implementation of the LWF-BROOK90 hydrological model in Julia

## Tell me more!

The hydrological model LWF-BROOK90 is a 1D Soil Vegetation Atmosphere Transport (SVAT) model, that includes vertical soil water movement; soil and plant evapotranspiration; and snowpack, interception and groundwater storages.
LWFBrook90Julia is an implementation of LWF-BROOK90 rewritten entirely in Julia.

- Processes and state variables in LWF-BROOK90 are summarised visually in Figure 1 below.
- LWFBrook90Julia can output additional quantities in daily, monthly or yearly resolution derived during simulation.
- For performance reasons computation of these additional quantities can also be deactivated during simulation. (The quantities can be computed post-simulation from the state vector (not yet implemented).)

LWFBrook90Julia is developed with the following objectives in mind:
- [ ] efficient parameter estimation
- [ ] model flexibility allowing to reparametrize processes (possibly resulting in a [flexible model framework](https://presentations.copernicus.org/EGU2020/EGU2020-17975_presentation.pdf))
- [ ] support for stable isotopes (δ18O and δ2H) by including transport equation and fractionation processes

<p align="center"><img src="https://github.com/fabern/LWFBrook90Julia.jl/blob/main/docs/src/figs/LWFBrook90Julia_overview_v1.4.png?raw=true" width="400"><p>
<p align="center">Figure 1: Summary of processes and state variables used in LWFBrook90Julia<p>


## Acknowledgement
Development of LWFBrook90Julia builds on the following works:
- Brook90 (v4.8) by C. Anthony Federer, Licensed under CC0 1.0, http://www.ecoshift.net/brook/brook90.htm
- LWFBrook90R (v0.4.3) by Paul Schmidt-Walter, Volodymyr Trotsiuk, Klaus Hammel, Martin Kennel, Anthony Federer, Robert Nuske, Licensed under GPL-3.0

See the file `LICENSE` for the license of LWFBrook90Julia.

Furthermore, Matthias Häni, Katrin Meusburger, Peter Waldner, Lorenz Walthert, Stephan Zimmermann of [WSL](www.wsl.ch) and the Long-term Forest Ecosystem Research (LWF) project of WSL are gratefully acknowledged for providing example data files located in `example/BEA2016*`.

## For Users: getting started
The file `main.jl` contains a setup and simple plotting command to get started with a first simulation.
An example data set `BEA2016*` was generated using the R package LWFBrook90R and is located inthe folder `example/`. See output at the bottom of this README.

Further documentation of [LWFBrook90Julia](https://fabern.github.io/LWFBrook90Julia.jl/stable)), and [BROOK90](http://www.ecoshift.net/brook/b90doc.html) are available by following the links.

## Limitations of LWFBrook90Julia
- Preprocessing steps to define a conevnient format of onput data sets are not yet designed. Currently, LWFBrook90Julia makes use of the preprocessing steps provded in the the R package LWFBrook90R.
- LWFBrook90R is based on LWF-BROOK90. LWF-BROOK90 was itself based an older version of Brook90 (v3.1F). Developments in Brook90 up to v4.8 (such as intercell averages of hydraulic conductivity) are therefore only partially included, but not activated by default in LWFBrook90Julia.
- Currently a part of LWFBrook90R Reset==1, is not implemented in LWFBrook90Julia. (A test implementation of this is available in the branch `005-simplify-time-step-control`.)

## For Developers:
Any help in form of discussions, pull requests, example data sets, or otherwise is greatly welcome. Please don't hesitate to contact us.

LWFBrook90Julia makes use of DifferentialEquations.jl to solve the system of Ordinary Differential Equations (ODE). Each state variable (`u`) has a corresponding ODE. The ODEs are defined by their right hand side defined in the function `f` (which sets `du` that is the change in `u`). The right hand side `f(u,p,t)` depends on time `t`, parameters `p` (time-dependent or constant), and the state `u`.

Variable naming generally follows the convention by BROOK90 and LWFBrook90R, but uses additionally a prefix to indicate their dependencies (`p_*`, `p_fT_*`, `p_fu_*` for constant, time dependent or state dependent parameters, `u_*`, `u_aux_*` for elementary and auxiliary state variables respectively, and `aux_du_*` for auxiliary rate of changes of state variables).

Note that some state variables (rain and snow interception storage, `u_INTR`, `u_INTS` and snow storage and energy `u_SNOW`, `u_CC`, `u_SNOWLQ`) are updated once per simulation day and other state variables (groundwater and soil water storages `u_GWAT`,`u_SWATI`) are solved on a higher resolved time discretization set by the ODE solver, resulting in a scheme based on operator splitting.

## Example data set:
Following plots illustrate results of the provided data set. The scalar state variables and depth-depenedent (vector) state variables can be plotted:

<p align="center"><img src="https://github.com/fabern/LWFBrook90Julia.jl/blob/main/docs/src/figs/git-hash-b3f7183/2021-02-24_16h56_LWFBrook90Julia_plot_u_scalar.png?raw=true" width="400"><p>
<p align="center">Figure 2: Example simulation: scalar results<p>

<p align="center"><img src="https://github.com/fabern/LWFBrook90Julia.jl/blob/main/docs/src/figs/git-hash-b3f7183/2021-02-24_16h56_LWFBrook90Julia_plot_u_vector.png?raw=true" width="400"><p>
<p align="center">Figure 3: Example simulation: vector results soil water<p>

Tests are run to assert agreement with results from LWFBrook90R. Currently minor discrepancies are still present (linked among other to the deactivation of "Reset==1" in LWFBrook90Julia and different time stepping algorithms).
<p align="center"><img src="https://github.com/fabern/LWFBrook90Julia.jl/blob/main/docs/src/figs/git-hash-b3f7183/2021-02-24_16h56_R-vs-Julia_comparison_DailyRawValues.png?raw=true" width="400"><p>
<p align="center">Figure 4: Comparing daily outputs of LWFBrook90R and LWFBrook90Julia for example data set over a year<p>

<p align="center"><img src="https://github.com/fabern/LWFBrook90Julia.jl/blob/main/docs/src/figs/git-hash-b3f7183/2021-02-24_16h56_R-vs-Julia_comparison_DailyRawValues_first3months.png?raw=true" width="400"><p>
<p align="center">Figure 5: Comparing daily outputs of LWFBrook90R and LWFBrook90Julia for example data set over 3 months<p>


### Improve agreement
Note that some features of LWFBrook90R are not implemented in the main version of LWFBrook90Julia. The time step adaptivity and Reset==1 are major ones that require some code refactoring that is not how the library for ODEs DiffEq.jl is intended to be used. Because of that implementation of these features is currently in a feature branch here on git `feature 005`. Below are some of the results of that code:

<p align="center"><img src="https://github.com/fabern/LWFBrook90Julia.jl/blob/main/docs/src/figs/git-hash-2-55ca42d-feature005/2021-02-24_19h10_R-vs-Julia_comparison_DailyRawValues.png?raw=true" width="400"><p>
<p align="center">Figure 6: Comparing daily outputs of LWFBrook90R and experimental LWFBrook90Julia:feature-005 for example data set over a year<p>

<p align="center"><img src="https://github.com/fabern/LWFBrook90Julia.jl/blob/main/docs/src/figs/git-hash-2-55ca42d-feature005/2021-02-24_19h10_R-vs-Julia_comparison_DailyRawValues_first3months.png?raw=true" width="400"><p>
<p align="center">Figure 7: Comparing daily outputs of LWFBrook90R and experimental LWFBrook90Julia:feature-005 for example data set over 3 months<p>
