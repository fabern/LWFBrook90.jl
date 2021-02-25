using LWFBrook90Julia
using Documenter

makedocs(;
    modules=[LWFBrook90Julia],
    authors="Fabian Bernhard",
    repo="https://github.com/fabern/LWFBrook90Julia.jl/blob/{commit}{path}#L{line}",
    sitename="LWFBrook90Julia.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fabern.github.io/LWFBrook90Julia.jl",
        assets=String[],
    ),
    pages=[
        "About"        => "index.md",
        "SVAT Model"   => "model.md",
        "User Guide"   => "user-guide.md",
        "Examples"     => "example.md",
        "Function References" => "function-refs.md",
        "Function Documentation" => "function-docs.md"
    ],
)

deploydocs(;
    repo="github.com/fabern/LWFBrook90Julia.jl",
    devbranch = "main"
)
