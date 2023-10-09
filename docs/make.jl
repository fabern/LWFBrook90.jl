# To build the documentation:
#    - julia --project=docs docs/make.jl
#    - empty!(ARGS); include("make.jl")
# To build the documentation without running the tests:
#    - julia --project="." make.jl preview
#    - push!(ARGS,"preview"); include("make.jl")

using
    Documenter,
    Literate,
    LWFBrook90,
    Dates


#### code pertaining to Literate.jl
if haskey(ENV, "GITHUB_ACTIONS")
    ### From Literate.jl (https://github.com/fredrikekre/Literate.jl/blob/master/docs/make.jl):
    ENV["JULIA_DEBUG"] = "Documenter"
    ENV["GKSwstype"] = "100" # following https://github.com/JuliaPlots/Plots.jl/issues/1182
end

deployconfig = Documenter.auto_detect_deploy_system()
Documenter.post_status(deployconfig; type="pending", repo="github.com/fabern/LWFBrook90.jl")

# a) define preprocessor for Literate.jl
function update_date(content)
    content = replace(content, "DATEOFTODAY" => Date(now()))
    return content
end

# b) generate examples
OUTPUT_markdown_for_Docs = joinpath(@__DIR__, "src/generated")
OUTPUT_scripts_for_examp = joinpath(@__DIR__, "src/generated")

# For Literate.jl use following formatting in the files:
    # # A markdown H1 title
    # A non-code markdown normal line
    ## A comment within the code chunk

            # Literate.markdown(joinpath(@__DIR__, "../examples/scripts/main.jl"), OUTPUT_markdown_for_Docs; name = "example", documenter = true, preprocess = update_date)
            # Literate.script(  joinpath(@__DIR__, "../examples/scripts/main.jl"), OUTPUT_scripts_for_examp; name = "example", documenter = true, preprocess = update_date)
            # Literate.notebook(joinpath(@__DIR__, "../examples/scripts/main.jl"), OUTPUT_folder; name = "example", documenter = true, preprocess = update_date)

#TODO: fix this:
Literate.markdown(joinpath(@__DIR__, "../examples/scripts/example-script-01.jl"), OUTPUT_markdown_for_Docs; name = "example-script-01", documenter = true, preprocess = update_date)
#TODO: fix this: Literate.script(  joinpath(@__DIR__, "../examples/scripts/example-script-01.jl"), OUTPUT_scripts_for_examp; name = "example-script-01", documenter = true, preprocess = update_date)
            # Literate.notebook(joinpath(@__DIR__, "../examples/scripts/example-script-01.jl"), OUTPUT_folder; name = "example-script-01", documenter = true, preprocess = update_date)
#### end code for Literate.jl


### code pertaining to Documenter.jl:
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
        "Example-01 (autogenerated)"  => "generated/example-script-01.md",
        "Validation"   => "example.md",
        "Code Listing" => "code-lst.md",
        "Function Documentations" => "function-docs.md",
    ],
)

deploydocs(;
    repo="github.com/fabern/LWFBrook90.jl",
    devbranch = "develop",
    devurl    = "dev",
    branch    = "gh-pages",
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
)
#### end code for Documenter.jl
