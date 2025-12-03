using Lattices
using Documenter

# ドキュメントを生成する
makedocs(
    sitename = "Lattices.jl",
    modules  = [Lattices],
    pages    = [
        "Home" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/sotashimozono/Lattices.jl.git",
)