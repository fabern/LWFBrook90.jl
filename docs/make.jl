using LWFBrook90
using Documenter

makedocs(;
    modules=[LWFBrook90],
    authors="Fabian Bernhard",
    repo="https://github.com/fabern/LWFBrook90.jl/blob/{commit}{path}#L{line}",
    sitename="LWFBrook90.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fabern.github.io/LWFBrook90.jl",
        assets=String[],
    ),
    pages=[
        "About"        => "index.md",
        "SVAT Model"   => "model.md",
        "User Guide"   => "user-guide.md",
        "Example"      => "example.md",
        "Code Listing" => "code-lst.md",
        "Function Documentations" => "function-docs.md"
    ],
)

deploydocs(;
    repo="github.com/fabern/LWFBrook90.jl",
    devbranch = "main",
    devurl    = "dev",
    branch    = "gh-pages"
)
