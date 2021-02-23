# fabian.bernhard@wsl.ch, 2021-02-01

module LWFBrook90Julia

using DifferentialEquations
#using Infiltrator

# on modules: https://discourse.julialang.org/t/large-programs-structuring-modules-include-such-that-to-increase-performance-and-readability/29102/5
include("module_CONSTANTS.jl");  # to bring into scope: using .CONSTANTS
include("module_GLBLDECL.jl");   # to bring into scope: using .GLBLDECL
include("module_KPT.jl");        # to bring into scope: using .KPT
include("module_WAT.jl");        # to bring into scope: using .WAT
include("module_SUN.jl");        # to bring into scope: using .SUN
include("module_PET.jl");        # to bring into scope: using .PET
include("module_SNO.jl");        # to bring into scope: using .SNO
include("module_EVP.jl");        # to bring into scope: using .SNO

export greet, read_LWFBrook90R_inputData, derive_params_from_input_meteo
export my_f

include("func_input_definition.jl")
include("func_DiffEq_definition_u0.jl")
include("func_DiffEq_definition_p.jl")
include("func_DiffEq_definition_cb.jl")
include("func_DiffEq_definition_f.jl")
include("func_DiffEq_definition_ode.jl")
include("func_MSB_functions.jl")

include("extra_file.jl")        #TODO(bernhard): remove this placeholder function once testing of own functions is implemented.
greet() = print("Hello World!") #TODO(bernhard): remove this placeholder function once testing of own functions is implemented.

end # module
