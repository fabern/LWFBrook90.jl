# fabian.bernhard@wsl.ch, 2021-01-02
# using Infiltrator
using LWFBrook90Julia
greet()


# 1a) Read in input data
(input_meteo,
    input_param,
    input_siteparam,
    input_precdat,    #TODO(benhard): input_precdat is unused
    input_pdur,
    input_soil_materials,
    input_soil_nodes) = read_KAUFENRING_inputData("../Brook90_Julia/Input_data_KAU/evergreen/")

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
                                                      input_soil_nodes)

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
p = LWFBrook90Julia.define_DiffEq_parameters(NLAYER, IMODEL, constant_dt_solver, NOOUTF, Reset,
                                             pfile_meteo, 
                                             pfile_siteparam,
                                             pfile_param, 
                                             pfile_soil, 
                                             pfile_pdur)
####################


####################
# Define initial states of differential equation
# state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI
# TODO(bernhard): are these store somewhere else than input_siteparam
u_GWAT_init = input_siteparam[1,5];
u_SNOW_init = input_siteparam[1,4];
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
    p_SWATMX = p[1][1][3] #p_cst_1[3]
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
global compute_intermediate_quantities = true # Flag whether ODE containes additional quantities than only states
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
#   - initial condition of states
#   - definition of simulation time span
#   - parameters
ode_LWFBrook90Julia = LWFBrook90Julia.define_DiffEq_ODE(u0, tspan, p, Reset)

####################
