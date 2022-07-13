using LWFBrook90
using Documenter
using Literate

#### code pertaining to Literate.jl
### From Literate.jl (https://github.com/fredrikekre/Literate.jl/blob/master/docs/make.jl):
if haskey(ENV, "GITHUB_ACTIONS")
    ENV["JULIA_DEBUG"] = "Documenter"
end

deployconfig = Documenter.auto_detect_deploy_system()
Documenter.post_status(deployconfig; type="pending", repo="github.com/fabern/LWFBrook90.jl")

# generate examples
EXAMPLE = joinpath(@__DIR__, "..", "examples", "scripts", "example.jl")
EXAMPLE = joinpath(@__DIR__, "..", "examples", "scripts", "test_literate.jl")
OUTPUT = joinpath(@__DIR__, "src/generated")

Literate.markdown(EXAMPLE, OUTPUT)
#### end code for Literate.jl


### From Documenter.jl:
DocMeta.setdocmeta!(LWFBrook90, :DocTestSetup, :(using LWFBrook90); recursive=true)

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
        "Function Documentations" => "function-docs.md",
        "Autogenerated Example" => "generated/example.md"
    ],
)

deploydocs(;
    repo="github.com/fabern/LWFBrook90.jl",
    devbranch = "develop",
    devurl    = "dev",
    branch    = "gh-pages",
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
)
