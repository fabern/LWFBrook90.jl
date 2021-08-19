using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5

# Read in input data
input_prefix = "BEA2016-reset-FALSE"
input_path = "example/"*input_prefix*"-input/"

####################
(input_meteoveg,
    input_meteoveg_reference_date,
    input_param,
    input_pdur,
    input_initial_conditions,
    input_soil_materials,
    input_soil_nodes) = read_LWFBrook90R_inputData(input_path, input_prefix)
####################

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

####################
# Define parameters for differential equation
p = define_LWFB90_p(
    input_meteoveg,
    input_meteoveg_reference_date,
    input_param,
    input_pdur,
    input_soil_materials,
    input_soil_nodes;
    Reset = Reset,
    compute_intermediate_quantities = compute_intermediate_quantities)
####################

####################
# Define initial states of differential equation
# state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI

######
# Transform initial value of auxiliary state u_aux_PSIM_init into state u_SWATIinit:
u_aux_PSIM_init = input_soil_nodes[:,"psiini_kPa"]
if any( u_aux_PSIM_init.> 0)
    error("Initial matrix psi must be negative or zero")
end
p_soil = p[1][1]
u_aux_WETNESinit = LWFBrook90.KPT.FWETNES(u_aux_PSIM_init, p_soil)
u_SWATIinit      = p_soil.p_SWATMAX ./ p_soil.p_THSAT .* LWFBrook90.KPT.FTheta(u_aux_WETNESinit, p_soil)
######

# Create u0 for DiffEq.jl
u0 = define_LWFB90_u0(input_initial_conditions[1,"u_GWAT_init"],
                      input_initial_conditions[1,"u_INTS_init"],
                      input_initial_conditions[1,"u_INTR_init"],
                      input_initial_conditions[1,"u_SNOW_init"],
                      input_initial_conditions[1,"u_CC_init"],
                      input_initial_conditions[1,"u_SNOWLQ_init"],
                      u_SWATIinit,
                      compute_intermediate_quantities)
####################

####################
# Define ODE problem which consists of
#   - definition of right-hand-side (RHS) function f
#   - definition of callback function cb
#   - u0:     initial condition of states
#   - tspan:  definition of simulation time span
#   - p:      parameters

# Define simulation time span:
# tspan = (0.,  5.) # simulate 5 days
# tspan = (0.,  100.) # simulate 100 days # NOTE: KAU bugs in "branch 005-" when at least 3*365
tspan = (minimum(input_meteoveg[:,"days"]),
         maximum(input_meteoveg[:,"days"])) # simulate all available days
# tspan = (LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1980,1,1), reference_date),
#          LWFBrook90.DateTime2RelativeDaysFloat(DateTime(1985,1,1), reference_date)) # simulates selected period

# Define ODE:
ode_LWFBrook90 = define_LWFB90_ODE(u0, tspan, p)
####################

####################
## Solve ODE:
sol_LWFBrook90 = solve(ode_LWFBrook90, Tsit5();
    progress = true,
    saveat = tspan[1]:tspan[2], dt=1e-6, adaptive = true); # dt is initial dt, but adaptive
####################

####################
## Benchmarking
# @time sol_LWFBrook90 = solve(ode_LWFBrook90, progress = true;
#     saveat = tspan[1]:tspan[2], dt=1e-1, adaptive = false);
# @time sol_LWFBrook90 = solve(ode_LWFBrook90, progress = true, Euler(); # Note: Euler sometimes hangs
#     saveat = tspan[1]:tspan[2], dt=1e-1, adaptive = false);
# using BenchmarkTools # for benchmarking
# sol_LWFBrook90 = @btime solve(ode_LWFBrook90; dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity
# sol_LWFBrook90 = @btime solve(ode_LWFBrook90; saveat = tspan[1]:tspan[2], dt=1.0e-1, adaptive = false); # dt will be overwritten, adaptive deacives DiffEq.jl adaptivity
####################


####################
## Plotting
# 1) very basic
# using Plots # Install plot package at first use with `]` and then `add Plots`
# # Plot 1
# plot(sol_LWFBrook90; vars = [1, 2, 3, 4, 5, 6],label=["GWAT" "INTS" "INTR" "SNOW" "CC" "SNOWLQ"])
# # Plot 2
# # http://docs.juliaplots.org/latest/generated/gr/#gr-ref43
# x = LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_meteoveg_reference_date)
# y = cumsum(p_soil.p_THICK)
# z = sol_LWFBrook90[7 .+ (0:p_soil.NLAYER-1), :]./p_soil.p_THICK
# heatmap(x, y, z, yflip = true,
#         xlabel = "Date",
#         ylabel = "Depth [mm]",
#         colorbar_title = "θ [-]")
####################
