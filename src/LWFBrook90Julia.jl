# fabian.bernhard@wsl.ch, 2021-02-01

module LWFBrook90Julia

export greet, input_definitions, my_f
include("extra_file.jl")

"""
    greet()
Return "Hello World!".
"""
greet() = print("Hello World!")

"""
input_definitions(folder::String, prefix::String)
Return "Define here input definitions...".
"""
function input_definitions(folder::String, prefix::String)
    return("Define here input definitions...")
end

end # module
