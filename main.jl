# fabian.bernhard@wsl.ch, 2021-01-02
# using Infiltrator
using DifferentialEquations
using LWFBrook90Julia
greet()


# 1a) Read in input data
(input_meteo,
    input_param,
    input_siteparam,
    input_precdat,    #TODO(benhard): input_precdat is unused
    input_pdur,
    input_soil_materials,
    input_soil_nodes,

    input_reference_date) = read_KAUFENRING_inputData("../Brook90_Julia/Input_data_KAU/evergreen/")

# 1b) Here posibility to modify dataframes
# TODO

# 1c) Parse loaded input files
(pfile_meteo, pfile_param, pfile_siteparam, pfile_precdat, pfile_pdur, pfile_soil) =
    LWFBrook90Julia.GLBLDECL.derive_params_from_input(input_meteo,
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
####################

####################
# Define parameters for differential equation
compute_intermediate_quantities = false # Flag whether ODE containes additional quantities than only states
p = LWFBrook90Julia.define_DiffEq_parameters(NLAYER, IMODEL, constant_dt_solver, NOOUTF, Reset, compute_intermediate_quantities,
                                             pfile_meteo,
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
u_SWATIinit  = fill(NaN, NLAYER)
if any( u_aux_PSIM_init.> 0)
    error("Initial matrix psi must be negative or zero")
end
if IMODEL == 0
    error("IMODEL==0 is not implemented to get initial SWATI.")
    # TODO(bernhard): implement this.
elseif IMODEL == 1
    # TODO(bernhard): are these store somewhere else than p[1][1][3] and pfile_soil
    p_SWATMX = p[1][1][6] #p_cst_1[3] # TODO(bernhard): this hardcoded index is dangerous in case definition of p vector changes
    p_MvGα   = pfile_soil["PAR"][!,"α"]
    p_MvGn   = pfile_soil["PAR"][!,"n"]
    p_THSAT  = pfile_soil["PAR"][!,"θs"]
    p_θr     = pfile_soil["PAR"][!,"θr"]
    for i = 1:NLAYER
        # Define initial u_SWATI based on input parameter
        u_aux_WETNESinit_i = LWFBrook90Julia.KPT.FWETNES_MvG(u_aux_PSIM_init[i], p_MvGα[i], p_MvGn[i])
        u_SWATIinit[i]     = LWFBrook90Julia.KPT.FTheta_MvG(u_aux_WETNESinit_i, p_THSAT[i], p_θr[i]) * p_SWATMX[i]/p_THSAT[i]
    end
end
######

# Create u0 for DiffEq.jl
u0 = LWFBrook90Julia.define_DiffEq_u0(u_GWAT_init,
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
tspan = (0.,  100.) # simulate 100 days
tspan = (0.,  5.) # simulate 5 days
# tspan = (minimum(input_meteo[:,"days"]),
#          maximum(input_meteo[:,"days"])) # simulate all available days
# tspan = (LWFBrook90Julia.DateTime2RelativeDaysFloat(DateTime(1980,1,1), reference_date),
#          LWFBrook90Julia.DateTime2RelativeDaysFloat(DateTime(1980,1,1), reference_date)) # simulates selected period

ode_LWFBrook90Julia = LWFBrook90Julia.define_DiffEq_ODE(u0, tspan, p)
sol_LWFBrook90Julia = solve(ode_LWFBrook90Julia, progress = true)

# using Plots
# plot(sol_LWFBrook90Julia; vars = [1, 2, 3, 4, 5, 6],label=["GWAT" "INTS" "INTR" "SNOW" "CC" "SNOWLQ"])
# plot(sol_LWFBrook90Julia; vars = 6+1:NLAYER,label="SWATI " .* string.(1:NLAYER))


# @time solve(ode_LWFBrook90Julia, progress = true)
####################