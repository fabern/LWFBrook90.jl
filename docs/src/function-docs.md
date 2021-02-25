```@meta
CurrentModule = LWFBrook90
```

# Function Documentations

## Functions from file `LWFBrook90.jl`
```@autodocs
Modules = [LWFBrook90]
Pages   = ["LWFBrook90.jl"]
```

## Functions from file `func_input_definition.jl`
```@autodocs
Modules = [LWFBrook90]
Pages   = ["func_input_definition.jl"]
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

## Functions from the different modules defining LWF-BROOK90
```@autodocs
Modules = [LWFBrook90.CONSTANTS,
           LWFBrook90.GLBLDECL,
           LWFBrook90.KPT,
           LWFBrook90.WAT,
           LWFBrook90.SUN,
           LWFBrook90.PET,
           LWFBrook90.SNO,
           LWFBrook90.EVP
           ]
Pages   = ["module_CONSTANTS.jl",
           "module_GLBLDECL.jl",
           "module_KPT.jl",
           "module_WAT.jl",
           "module_SUN.jl",
           "module_PET.jl",
           "module_SNO.jl",
           "module_EVP.jl"
           ]
Order = [:constant, :type, :function, :macro]
```