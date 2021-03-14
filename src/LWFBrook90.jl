module LWFBrook90

using OrdinaryDiffEq  # instead of loading the full DifferentialEquations
using DiffEqCallbacks # instead of loading the full DifferentialEquations

export read_LWFBrook90R_inputData, derive_params_from_inputData
export define_LWFB90_p, define_LWFB90_u0, define_LWFB90_ODE
export KPT_SOILPAR_Mvg1d, KPT_SOILPAR_Ch1d

# on modules: https://discourse.julialang.org/t/large-programs-structuring-modules-include-such-that-to-increase-performance-and-readability/29102/5
include("module_CONSTANTS.jl");  # to bring into scope: using .CONSTANTS
include("module_KPT.jl");        using .KPT # using to bring exports into scope
include("module_WAT.jl");        using .WAT # using to bring exports into scope
include("module_SUN.jl");        # to bring into scope: using .SUN
include("module_PET.jl");        # to bring into scope: using .PET
include("module_SNO.jl");        # to bring into scope: using .SNO
include("module_EVP.jl");        # to bring into scope: using .SNO
include("module_GLBLDECL.jl");   using .GLBLDECL # using to bring exports into scope

include("func_input_definition.jl")
include("func_DiffEq_definition_u0.jl")
include("func_DiffEq_definition_p.jl")
include("func_DiffEq_definition_cb.jl")
include("func_DiffEq_definition_f.jl")
include("func_DiffEq_definition_ode.jl")
include("func_MSB_functions.jl")

include("../example/BEA2016-reset-FALSE-input/func_run_example.jl")

end # module
