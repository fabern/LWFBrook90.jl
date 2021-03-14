using LWFBrook90
using OrdinaryDiffEq: solve, Tsit5

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

# 1b) Here posibility to modify dataframes input_[...] manually
# TODO

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

####################
# Define simulation
# Soil hydraulic model
IMODEL = pfile_param["IMODEL"] # 0 for Clapp-Hornberger; 1 for Mualem-van Genuchten
NLAYER = pfile_param["NLAYER"]

# Define solver options
NOOUTF    = 1 == pfile_param["NOOUTF"] # 1 if no outflow allowed from roots, otherwise 0
Reset     = 0 # currently only Reset = 0 implemented

constant_dt_solver = 1 # [days]

compute_intermediate_quantities = true # Flag whether ODE containes additional quantities than only states
####################

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
    saveat = tspan[1]:tspan[2], dt=1e-6, adaptive = true); # dt will be overwritten
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
# x = LWFBrook90.RelativeDaysFloat2DateTime.(sol_LWFBrook90.t, input_reference_date)
# y = cumsum(pfile_soil["THICK"])
# z = sol_LWFBrook90[6 .+ (1:NLAYER), :]./pfile_soil["THICK"]
# heatmap(x, y, z, yflip = true,
#         xlabel = "Date",
#         ylabel = "Depth",
#         colorbar_title = "θ")
####################
