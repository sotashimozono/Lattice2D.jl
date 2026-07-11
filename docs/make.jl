using Lattice2D
using LatticeCore
using Documenter

# Run the feature figure generator script automatically before building docs
@info "Generating geometric feature figures..."
include(joinpath(@__DIR__, "src", "generate_features.jl"))

# ドキュメントを生成する
makedocs(;
    sitename="Lattice2D.jl",
    # Only Lattice2D is in `modules` (drives the missing-docs check). LatticeCore
    # symbols are surfaced via the `@autodocs Modules = [LatticeCore]` block in
    # latticecore.md purely so their `@ref` cross-references resolve (issue #68);
    # we intentionally do not require *every* LatticeCore docstring to be present.
    modules=[Lattice2D],
    pages=[
        "Home" => "index.md",
        "Gallery" => "Gallery.md",
        "LatticeCore API" => "latticecore.md",
    ],
)

deploydocs(; repo="github.com/sotashimozono/Lattice2D.jl.git")
