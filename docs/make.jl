using PfaffianSystems
using Documenter

DocMeta.setdocmeta!(PfaffianSystems, :DocTestSetup, :(using PfaffianSystems); recursive=true)

makedocs(;
    modules=[PfaffianSystems],
    authors="t0rny",
    repo="https://github.com/t0rny/PfaffianSystems.jl/blob/{commit}{path}#{line}",
    sitename="PfaffianSystems.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://t0rny.github.io/PfaffianSystems.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/t0rny/PfaffianSystems.jl",
    devbranch="main",
)
