# LWFBrook90.jl
Package repository: [https://github.com/fabern/LWFBrook90.jl](https://github.com/fabern/LWFBrook90.jl)

```@contents
Pages = ["index.md",
         "model.md",
         "user-guide.md",
         "example.md",
         "code-lst.md",
         "function-docs.md"]
Depth = 2
```



## About LWFBrook90.jl
LWFBrook90.jl implements a 1D soil-vegetation-atmosphere transfer model.
Intended use cases of this impelemntation are:
- support of stable isotopes (δ¹⁸O and δ²H)
- efficient model calibration for data analysis
- increased model flexibility

To read about the model structure, see section [SVAT Model](@ref).
For a quick start refer to the step-by-step guide in section [Example](@ref)
For further details read through the section [User Guide](@ref) or refer to sections [Code Listing](@ref) and [Function Documentations](@ref) for further technical intricacies e.g. for development of the package.



## Citing LWFBrook90.jl
When using LWFBrook90.jl please cite <!-- [TODO.et.al (2021)](TODO) -->
>TODO: generate citation: (Journal of Open Source Software? Zenodo? Journal of Open Research Software (DifferentialEquations.jl)? Alternatives?)



## Acknowledgments
Development of LWFBrook90.jl builds on the following works:
- BROOK90 (v4.8) by C. Anthony Federer, Licensed under CC0 1.0, http://www.ecoshift.net/brook/brook90.htm
- LWFBrook90R (v0.4.3) by Paul Schmidt-Walter, Volodymyr Trotsiuk, Klaus Hammel, Martin Kennel, Anthony Federer, Robert Nuske, Licensed under GPL-3.0, https://github.com/pschmidtwalter/LWFBrook90R

Matthias Häni, Katrin Meusburger, Peter Waldner, Lorenz Walthert, Stephan Zimmermann of [WSL](http://www.wsl.ch) and its Long-term Forest Ecosystem Research (LWF) project gratefully acknowledged for providing example data files located in `example/BEA2016*`.


## References
Federer, C. A., Vörösmarty, C., & Fekete, B. (2003). Sensitivity of Annual Evaporation to Soil and Root Properties in Two Models of Contrasting Complexity. *Journal of Hydrometeorology*, 4(6), 1276–1290. [https://doi.org/10.1175/1525-7541(2003)004<1276:SOAETS>2.0.CO;2](https://doi.org/10.1175/1525-7541(2003)004<1276:SOAETS>2.0.CO;2)

Federer, C. A. (2002). BROOK 90: A simulation model output for evaporation, soil water, and streamflow. http://www.ecoshift.net/brook/brook90.htm

Hammel, K., & Kennel, M. (2001). Charakterisierung und Analyse der Wasserverfügbarkeit und des Wasserhaushalts von Waldstandorten in Bayern mit dem Simulationsmodell BROOK90 (No. 185; *Forstliche Forschungsberichte München*, p. 135). Technische Uni München Wissenschaftszentrum Weihenstephan. ISBN 3-933506-16-6

Schmidt-Walter, P., Trotsiuk, V., Meusburger, K., Zacios, M., & Meesenburg, H. (2020). Advancing simulations of water fluxes, soil moisture and drought stress by using the LWF-Brook90 hydrological model in R. *Agricultural and Forest Meteorology*, 291, 108023. https://doi.org/10.1016/j.agrformet.2020.108023