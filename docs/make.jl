##
using Documenter,RTI

makedocs(
    modules=[RTI],
    authors="Yi Zhang",
    sitename="RTI",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    pages=[
        "Home" => "index.md",
        "Forward problems"=>"Forward problems.md",
        "Inverse problems" => Any[
        "A toy model" =>"A toy model.md"
        "Crati"=>"Crati.md"],
        "Utilities"=>"Utilities.md",
        "List of functions" =>"list_of_functions.md"
        ],
)

deploydocs(;
    repo="github.com/deconvolution/RTI.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "main",
    devurl = "dev",
    forcepush=true,
    push_preview =true
)
