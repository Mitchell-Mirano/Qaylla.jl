using Documenter, Yaqa

makedocs(
    sitename = "Yaqa.jl",
    authors  = "Mitchell Mirano",
    modules=[Yaqa],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
    pages = [
        "Home" => "index.md",
        "Interpolations" => "interpolations.md",
        "Integrals" => "integrals.md",
        "Differential Equations" => "differential-equations.md"
    ]
)

deploydocs(
    repo = "github.com/Mitchell-Mirano/Yaqa.jl.git"
)