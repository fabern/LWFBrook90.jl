using LWFBrook90
using Documenter
using Literate

### From Literate.jl (https://github.com/fredrikekre/Literate.jl/blob/master/docs/make.jl):
if haskey(ENV, "GITHUB_ACTIONS")
    ENV["JULIA_DEBUG"] = "Documenter"
end

deployconfig = Documenter.auto_detect_deploy_system()
Documenter.post_status(deployconfig; type="pending", repo="github.com/fabern/LWFBrook90.jl")

# generate examples
EXAMPLE = joinpath(@__DIR__, "..", "examples", "example.jl")
OUTPUT = joinpath(@__DIR__, "src/generated")
function preprocess(str)
    str = replace(str, "x = 123" => "y = 321"; count=1)
    return str
end

Literate.markdown(EXAMPLE, OUTPUT, preprocess = preprocess)
Literate.notebook(EXAMPLE, OUTPUT, preprocess = preprocess)
Literate.script(EXAMPLE, OUTPUT, preprocess = preprocess)

# Replace the link in outputformats.md
# since that page is not "literated"
if haskey(ENV, "GITHUB_ACTIONS")
    folder = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        Documenter.deploy_folder(
            deployconfig;
            repo="github.com/fabern/LWFBrook90.jl",
            devbranch = "develop",
            push_preview = true,
            devurl = "dev",
        ).subfolder
    end
    url = "https://nbviewer.jupyter.org/github/fabern/LWFBrook90.jl/blob/gh-pages/$(folder)/"
    str = read(joinpath(@__DIR__, "src/outputformats.md"), String)
    str = replace(str, "[notebook.ipynb](generated/notebook.ipynb)." => "[notebook.ipynb]($(url)generated/notebook.ipynb).")
    write(joinpath(@__DIR__, "src/outputformats.md"), str)
end



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
        "Function Documentations" => "function-docs.md"
    ],
)

deploydocs(;
    repo="github.com/fabern/LWFBrook90.jl",
    devbranch = "develop",
    devurl    = "dev",
    branch    = "gh-pages",
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
)
