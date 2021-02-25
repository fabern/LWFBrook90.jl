# LWFBrook90Julia.jl
Package repository: [https://github.com/fabern/LWFBrook90Julia.jl](https://github.com/fabern/LWFBrook90Julia.jl)

```@contents
Pages = ["index.md",
         "model.md",
         "user-guide.md",
         "example.md",
         "function-refs.md",
         "function-docs.md"]
Depth = 2
```



## About LWFBrook90Julia.jl
LWFBrook90Julia.jl implements a 1D soil vegetation atmosphere transport model.
Intended use cases of this impelemntation are:
- efficient model calibration for data analysis
- support of stable isotopes (δ¹⁸O and δ²H)
- increased flexibility for reparametrizations

To read about the model structure, see section [SVAT Model](@ref).
For a quick start refer to the step-by-step guide in section [Example](@ref)
For further details read through the [User Guide](@ref) or refer to sections [Function References](@ref) and [Function Documentations](@ref) for further technical intricacies e.g. for development of the package.



## Citing LWFBrook90Julia.jl
When using LWFBrook90Julia.jl please cite <!-- [TODO.et.al (2021)](TODO) -->
>TODO: generate citation: (Journal of Open Source Software? Zenodo? Alternatives?)



## Acknowledgments
- Brook90 (v4.8) by C. Anthony Federer
- LWFBrook90R (v0.4.3) by Paul Schmidt-Walter, Volodymyr Trotsiuk, Klaus Hammel, Martin Kennel, Anthony Federer, Robert Nuske
- Matthias Häni, Katrin Meusburger, Peter Waldner, Lorenz Walthert, Stephan Zimmermann of [WSL](https://www.wsl.ch) and the Long-term Forest Ecosystem Research (LWF) project of WSL are gratefully acknowledged for providing example data files located in `example/BEA2016*`.



## References
Federer, C. A., Vörösmarty, C., & Fekete, B. (2003). Sensitivity of Annual Evaporation to Soil and Root Properties in Two Models of Contrasting Complexity. *Journal of Hydrometeorology*, 4(6), 1276–1290. https://doi.org/10.1175/1525-7541(2003)004<1276:SOAETS>2.0.CO;2

Federer, C. A. (2002). BROOK 90: A simulation model output for evaporation, soil water, and streamflow. http://www.ecoshift.net/brook/brook90.htm

Hammel, K., & Kennel, M. (2001). Charakterisierung und Analyse der Wasserverfügbarkeit und des Wasserhaushalts von Waldstandorten in Bayern mit dem Simulationsmodell BROOK90 (No. 185; *Forstliche Forschungsberichte München*, p. 135). Technische Uni München Wissenschaftszentrum Weihenstephan. ISBN 3-933506-16-6

Hammel, K., & Kennel, M. (2001). Charakterisierung und Analyse der Wasserverfügbarkeit und des Wasserhaushalts von Waldstandorten in Bayern mit dem Simulationsmodell BROOK90 (No. 185; *Forstliche Forschungsberichte München*, p. 135). Technische Uni München Wissenschaftszentrum Weihenstephan. ISBN 3-933506-16-6

Schmidt-Walter, P., Trotsiuk, V., Meusburger, K., Zacios, M., & Meesenburg, H. (2020). Advancing simulations of water fluxes, soil moisture and drought stress by using the LWF-Brook90 hydrological model in R. *Agricultural and Forest Meteorology*, 291, 108023. https://doi.org/10.1016/j.agrformet.2020.108023