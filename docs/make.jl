ENV["GKSwstype"] = "100"
using RandomKrylov
using Documenter

DocMeta.setdocmeta!(RandomKrylov, :DocTestSetup, :(using RandomKrylov); recursive = true)

makedocs(;
    modules=[RandomKrylov],
    authors="Tangi Migot tangi.migot@gmail.com",
    repo="https://github.com/tmigot/RandomKrylov.jl/blob/{commit}{path}#{line}",
    sitename="RandomKrylov.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tmigot.github.io/RandomKrylov.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Benchmark overdetermined" => "rectangular.md",
        "Benchmark sdp matrices" => "sdp.md",
        "Benchmark MatrixDepot" => "matrixdepot.md",
    ],
)

deploydocs(; repo = "github.com/tmigot/RandomKrylov.jl", devbranch = "main")
