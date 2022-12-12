using Documenter, Qaylla

makedocs(
    sitename = "Qaylla.jl",
    authors  = "Mitchell Mirano",
    modules=[Qaylla],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
    pages = [
        "Home" => "index.md",
        "Interpolations" => "interpolations.md",
        "Integrals" => "integrals.md",
        "Differential Equations" => "differential-equations.md"
    ]
)

deploydocs(
    repo = "github.com/Mitchell-Mirano/Qaylla.jl.git"
)