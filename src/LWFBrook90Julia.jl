# fabian.bernhard@wsl.ch, 2021-02-01

module LWFBrook90Julia

include("module_CONSTANTS.jl"); #using CONSTANTS
include("module_GLBLDECL.jl");   #using GLBLDECL
include("module_KPT.jl");       #using KPT
include("module_WAT.jl");       #using WAT
include("module_SUN.jl");       #using SUN

export greet, read_KAUFENRING_inputData, derive_params_from_input_meteo
export my_f

include("func_input_definition.jl")
include("func_DiffEq_definition.jl")

include("extra_file.jl")        #TODO(bernhard): remove this placeholder function.
greet() = print("Hello World!") #TODO(bernhard): remove this placeholder function.

end # module
