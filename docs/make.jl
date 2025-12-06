using Lattice2D
using Documenter

# ドキュメントを生成する
makedocs(
    sitename = "Lattice2D.jl",
    modules  = [Lattice2D],
    pages    = [
        "Home" => "index.md",
        "Gallery" => "Gallery.md"
    ]
)

deploydocs(
    repo = "github.com/sotashimozono/Lattice2D.jl.git",
)