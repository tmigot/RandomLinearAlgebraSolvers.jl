ENV["GKSwstype"] = "100"
using RandomLinearAlgebraSolvers
using Documenter

DocMeta.setdocmeta!(RandomLinearAlgebraSolvers, :DocTestSetup, :(using RandomLinearAlgebraSolvers); recursive = true)

makedocs(;
    modules=[RandomLinearAlgebraSolvers],
    authors="Tangi Migot tangi.migot@gmail.com",
    repo="https://github.com/tmigot/RandomLinearAlgebraSolvers.jl/blob/{commit}{path}#{line}",
    sitename="RandomLinearAlgebraSolvers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tmigot.github.io/RandomLinearAlgebraSolvers.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Benchmark overdetermined" => "rectangular.md",
        "Benchmark sdp matrices" => "sdp.md",
        "Benchmark MatrixDepot" => "matrixdepot.md",
    ],
)

deploydocs(; repo = "github.com/tmigot/RandomLinearAlgebraSolvers.jl", devbranch = "main")
