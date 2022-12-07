using Documenter, NumericalMethods

makedocs(
    modules = [NumericalMethods],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "NumericalMethods.jl",
    authors  = "Mitchell Mirano",
    pages = [
        "Home" => "index.md",
        "Interpolations" => "interpolations.md",
        "Integrals" => "integrals.md",
        "Differential Equations" => "differential-equations.md"
    ]
)

deploydocs(
    repo = "github.com/Mitchell-Mirano/NumericalMethods.jl.git",
)