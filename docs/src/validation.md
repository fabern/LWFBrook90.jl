# Validation

## Comparison with LWFBrook90R
Tests are run to assert agreement with results from LWFBrook90R. Visualizations are reported below. Note that minor discrepancies are still present linked to the adaptive time stepping and intermediate updates of state variables.
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

Note that some features of LWFBrook90R are not implemented in the main version of LWFBrook90.jl. The time step adaptivity and `Reset==1` are major ones that require some code refactoring in order to use them together with the library for ODEs DiffEq.jl. Because of that implementation of these features is currently unsupported.