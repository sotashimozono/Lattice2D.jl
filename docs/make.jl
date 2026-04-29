using Lattice2D
using Documenter

# Run the feature figure generator script automatically before building docs
@info "Generating geometric feature figures..."
include(joinpath(@__DIR__, "src", "generate_features.jl"))

# ドキュメントを生成する
makedocs(;
    sitename="Lattice2D.jl",
    modules=[Lattice2D],
    pages=["Home" => "index.md", "Gallery" => "Gallery.md"],
)

deploydocs(; repo="github.com/sotashimozono/Lattice2D.jl.git")
