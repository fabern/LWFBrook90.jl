# fabian.bernhard@wsl.ch, 2021-02-01

module LWFBrook90Julia

using DifferentialEquations

# on modules: https://discourse.julialang.org/t/large-programs-structuring-modules-include-such-that-to-increase-performance-and-readability/29102/5
include("module_CONSTANTS.jl");  # to bring into scope: using .CONSTANTS
include("module_GLBLDECL.jl");   # to bring into scope: using .GLBLDECL
include("module_KPT.jl");        # to bring into scope: using .KPT
include("module_WAT.jl");        # to bring into scope: using .WAT
include("module_SUN.jl");        # to bring into scope: using .SUN
include("module_PET.jl");        # to bring into scope: using .PET
include("module_SNO.jl");        # to bring into scope: using .SNO
include("module_EVP.jl");        # to bring into scope: using .SNO

export greet, read_KAUFENRING_inputData, derive_params_from_input_meteo
export my_f

include("func_input_definition.jl")
include("func_DiffEq_definition.jl")
include("func_MSB_functions.jl")

include("extra_file.jl")        #TODO(bernhard): remove this placeholder function.
greet() = print("Hello World!") #TODO(bernhard): remove this placeholder function.

end # module
