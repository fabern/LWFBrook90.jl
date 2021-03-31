# LWFBrook90.jl
Implementation of the LWF-BROOK90 hydrological model in Julia

<p align="center">User info:
<a href="https://www.repostatus.org/#wip"><img title="Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public." src="https://www.repostatus.org/badges/latest/wip.svg"/></a>
<a href="https://www.gnu.org/licenses/gpl-3.0"><img title="License: GPL v3" src="https://img.shields.io/badge/License-GPLv3-blue.svg"></a>
<a href="https://gitter.im/fabern/LWFBrook90.jl"><img title="Join the chat at https://gitter.im/fabern/LWFBrook90.jl" src="https://badges.gitter.im/fabern/LWFBrook90.jl.svg"></a>
<a href="https://fabern.github.io/LWFBrook90.jl/stable"><img title="Stable" src="https://img.shields.io/badge/docs-stable-blue.svg"></a>
</p>
<p align="center">Developer info:
<a href="https://github.com/fabern/LWFBrook90.jl/actions"><img title="Build Status" src="https://github.com/fabern/LWFBrook90.jl/workflows/CI/badge.svg"/></a>
<a href="https://codecov.io/gh/fabern/LWFBrook90.jl"><img title="Coverage" src="https://codecov.io/gh/fabern/LWFBrook90.jl/branch/develop/graph/badge.svg?token=6F1BUJ8UR9"/>
<a href="https://github.com/invenia/BlueStyle"><img title="Code Style: Blue" src="https://img.shields.io/badge/code%20style-blue-4495d1.svg"/></a>
<a href="https://fabern.github.io/LWFBrook90.jl/dev"><img title="Dev" src="https://img.shields.io/badge/docs-dev-blue.svg"></a>
<a href="https://discourse.julialang.org/t/good-practices-for-package-development-in-the-julia-ecosystem/8175/2">Git branch model</a>
</p>

## What is LWFBrook90.jl?
<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/logo.png?raw=true" width="150"><p>
The model LWF-BROOK90 is a 1D Soil-Vegetation-Atmosphere Transfer (SVAT) model, calculating the soil water balance in forest soil. Modelled processes include vertical soil water movement, soil and plant evapotranspiration and temporary storages in snowpack or interception layer.
LWFBrook90.jl is an implementation of the existing LWF-BROOK90 rewritten entirely in the Julia programming language.

Processes and state variables of the model in LWF-BROOK90 are summarised visually in Figure 1 below.
LWFBrook90.jl is developed with the following objectives in mind:
- [ ] support for stable isotopes (δ¹⁸O and δ²H) by including transport equation and fractionation processes
- [ ] efficient parameter estimation (optimizing computational costs)
- [ ] model flexibility for alternative processes parametrizations (possibly resulting in a [flexible model framework](https://presentations.copernicus.org/EGU2020/EGU2020-17975_presentation.pdf))



<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/LWFBrook90jl_overview_v1.4-juliacolors.png?raw=true" width="400"><p>
<p align="center">Figure 1: Summary of processes and state variables used in LWFBrook90.jl<p>

Further model description can be found in the [documentation](https://fabern.github.io/LWFBrook90.jl/stable).


## License
LWFBrook90.jl is licensed under GPL-3.0 (for details see the file `LICENSE`).

## Acknowledgements
Development of LWFBrook90.jl builds on the following works:
- BROOK90 (v4.8) by C. Anthony Federer, Licensed under CC0 1.0, http://www.ecoshift.net/brook/brook90.htm
- LWFBrook90R (v0.4.3) by Paul Schmidt-Walter, Volodymyr Trotsiuk, Klaus Hammel, Martin Kennel, Anthony Federer, Robert Nuske, Licensed under GPL-3.0, https://github.com/pschmidtwalter/LWFBrook90R

Note that LWFBrook90R itself uses LWF-BROOK90 by Hammel and Kennel, 2001:
Hammel, K., & Kennel, M. (2001). Charakterisierung und Analyse der Wasserverfügbarkeit und des Wasserhaushalts von Waldstandorten in Bayern mit dem Simulationsmodell BROOK90 (No. 185; *Forstliche Forschungsberichte München*, p. 135). Technische Uni München Wissenschaftszentrum Weihenstephan. ISBN 3-933506-16-6

All literature references are reported in the section [References](https://fabern.github.io/LWFBrook90.jl/stable/#References) in the documentation.

Matthias Häni, Katrin Meusburger, Peter Waldner, Lorenz Walthert, Stephan Zimmermann of [WSL](http://www.wsl.ch) and its Long-term Forest Ecosystem Research (LWF) project gratefully acknowledged for providing example data files located in `example/BEA2016*`.

## For Users: getting started
To get started with Julia: see the section [Installation](https://fabern.github.io/LWFBrook90.jl/stable/user-guide/#Installation) in the documentation.

An example data set `BEA2016*` was generated using the R package LWFBrook90R and is located in the folder `example/`. See outputs in the next section.

To run this example simulation simply call `LWFBrook90.run_example()`. Note, that the first time might take some time to load and compile the package. Another possibility is to follow the script `main.jl`, which sets up a simulation and shows some simple plotting commands. See documentation section [Example](https://fabern.github.io/LWFBrook90.jl/stable/example/).

Further documentation for LWFBrook90.jl is available [here](https://fabern.github.io/LWFBrook90.jl/stable) and for BROOK90 [here](http://www.ecoshift.net/brook/b90doc.html).

### Example data set:
Following plots illustrate results of the provided data set. The scalar state variables and depth-depenedent (vector) state variables can be plotted:

<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/example-results/2021-03-16_11h27-Reset0-git+6e6946d+clean_LWFBrook90Julia_plot_u_scalar.png?raw=true" width="400"><p>
<p align="center">Figure 2: Example simulation: scalar results<p>

<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/example-results/2021-03-16_11h27-Reset0-git+6e6946d+clean_LWFBrook90Julia_plot_u_vector.png?raw=true" width="400"><p>
<p align="center">Figure 3: Example simulation: vector results soil water<p>

Tests are run to assert agreement with results from LWFBrook90R. Visualizations are reported below. Note that minor discrepancies are still present linked to the adaptive time stepping and intermediate updates of state variables.
<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/example-results/2021-03-16_11h27-Reset0-git+6e6946d+clean_R-vs-Julia_Daily_first12M.png?raw=true" width="400"><p>
<p align="center">Figure 4: Comparing daily outputs of LWFBrook90R and LWFBrook90.jl for example data set over a year<p>

<p align="center"><img src="https://github.com/fabern/LWFBrook90.jl/blob/develop/docs/src/assets/example-results/2021-03-16_11h27-Reset0-git+6e6946d+clean_R-vs-Julia_Daily_first2M.png?raw=true" width="400"><p>
<p align="center">Figure 5: Comparing daily outputs of LWFBrook90R and LWFBrook90.jl for example data set over 2 months<p>



## Limitations of LWFBrook90.jl
- Preprocessing steps to define a convenient format of input data sets are not yet designed. Currently, LWFBrook90.jl makes use of an input data set that is best prepared using the preprocessing steps provded in the the R package LWFBrook90R.
- Newest developments of BROOK90 are not included. LWFBrook90R is based on LWF-BROOK90, which was itself forked from an older vresion of BROOK90 (v3.1F). Developments in BROOK90 up to v4.8 (such as `KKMEAN` the intercell averages of hydraulic conductivity) are partially included, but not activated by default in LWFBrook90.jl.
- Currently a part of LWFBrook90R activated when `Reset==1`, is not implemented in LWFBrook90.jl. (A test implementation of this is available in the branch `005-simplify-time-step-control`.)

### Improve agreement with LWFBrook90R
Note that some features of LWFBrook90R are not implemented in the main version of LWFBrook90.jl. The time step adaptivity and `Reset==1` of LWFBrook90R would require code refactoring, which goes slightly against the intended use of the library for ODEs DifferentialEquations.jl. Because of this, implementation of these features is left away from the main version. However, an attempt at their implementation resides currently in a feature branch `005`.



## For Developers:
Any help in form of discussions, pull requests, example data sets, or otherwise is very welcome. Please don't hesitate to contact us.

For implementation details: see documentation section [Implementation](https://fabern.github.io/LWFBrook90.jl/stable/model/#Implementation).

If you're new to scientific computing in Julia. There are many useful tutorials around. For a quick start you can check out:
- Getting Started with Julia (for Experienced Programmers), by C.Rackauckas, 35min: https://www.youtube.com/watch?v=-lJK92bEKow
- Noteworthy Differences from other Languages: https://docs.julialang.org/en/v1/manual/noteworthy-differences/

For more in depth treatments of topics related to scientific computing and machine learning in Julia be sure to check out the course MIT 18.337J by C. Rackauckas at https://github.com/mitmath/18337. This course lists further introduction material under lecture 1.1.