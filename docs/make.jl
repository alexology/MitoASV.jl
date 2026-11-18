using Documenter
using MitoASV

DocMeta.setdocmeta!(MitoASV, :DocTestSetup, :(using MitoASV); recursive=true)

makedocs(;
    modules=[MitoASV],
    sitename="MitoASV.jl",
    authors="Alex Laini",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alexology.github.io/MitoASV.jl",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "Analyses" => "analyses.md",
        "Output tables and flags" => "output.md",
        "Settings and complexity" => "settings.md",
        "Methods" => "methods.md",
        "Limitations and testing" => "limitations.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(; repo="github.com/alexology/MitoASV.jl", devbranch="main")
