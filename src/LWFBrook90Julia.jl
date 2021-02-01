# fabian.bernhard@wsl.ch, 2021-02-01

module LWFBrook90Julia

export greet, input_definitions, my_f

include("extra_file.jl")


greet() = print("Hello World!")

function input_definitions(folder, prefix)
    return("Define here input definitions...")
end
end # module
