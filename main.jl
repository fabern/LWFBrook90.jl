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

# TODO(bernhard): clean up the four above functions...

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


