```@meta
CurrentModule = LWFBrook90
```

# Function Documentations

## Functions from file `examples/BEA2016...-input/func_run_example.jl`
```@autodocs
Modules = [LWFBrook90]
Pages   = ["examples/func_run_example.jl"]
```

## Functions from file `LWFBrook90.jl`
```@autodocs
Modules = [LWFBrook90]
Pages   = ["LWFBrook90.jl"]
```

## Functions from files `func_read_inputData.jl`, `func_discretize_soil_domain.jl`, or `func_postprocess.jl`
```@autodocs
Modules = [LWFBrook90]
Pages   = ["func_read_inputData.jl",
           "func_discretize_soil_domain.jl",
           "func_postprocess.jl",
           "generate_LWFBrook90R_Input.jl"]
```

## Functions defining the DiffEq.jl system of ODE (p, u0, f, callbacks, ...)
```@autodocs
Modules = [LWFBrook90]
Pages   = ["func_DiffEq_definition_u0.jl",
           "func_DiffEq_definition_p.jl",
           "func_DiffEq_definition_cb.jl",
           "func_DiffEq_definition_f.jl",
           "func_DiffEq_definition_ode.jl",
           "func_MSB_functions.jl"
           ]
```

## Functions from the different modules defining LWFBrook90
```@autodocs
Modules = [LWFBrook90.CONSTANTS,
           LWFBrook90.KPT,
           LWFBrook90.WAT,
           LWFBrook90.SUN,
           LWFBrook90.PET,
           LWFBrook90.SNO,
           LWFBrook90.EVP,
           LWFBrook90.ISO,
           ]
Pages   = ["module_CONSTANTS.jl",
           "module_KPT.jl",
           "module_WAT.jl",
           "module_SUN.jl",
           "module_PET.jl",
           "module_SNO.jl",
           "module_EVP.jl",
           "module_ISO.jl"
           ]
Order = [:constant, :type, :function, :macro]
```