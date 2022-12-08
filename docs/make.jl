using Documenter, NumericalMethods

makedocs(
    sitename = "NumericalMethods.jl",
    authors  = "Mitchell Mirano",
    modules=[NumericalMethods],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
    pages = [
        "Home" => "index.md",
        "Interpolations" => "interpolations.md",
        "Integrals" => "integrals.md",
        "Differential Equations" => "differential-equations.md"
    ]
)

deploydocs(
    repo = "github.com/Mitchell-Mirano/NumericalMethods.jl.git"
)