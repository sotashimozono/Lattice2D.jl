ENV["GKSwstype"] = "100"

using Lattices, Test, Plots
using LinearAlgebra

const FIG_BASE = joinpath(pkgdir(Lattices), "docs", "src", "assets", "figures")
const FIG_LAT = joinpath(FIG_BASE, "lattice")
const PATHS = Dict(
    :geometry => joinpath(FIG_LAT, "geometry")
)
mkpath.(values(PATHS))

const dirs = ["core", "lattices", "utils"]

@testset "tests" begin
    test_args = copy(ARGS)
    println("Passed arguments ARGS = $(test_args) to tests.")
    @time for dir in dirs
        dirpath = joinpath(@__DIR__, dir)
        println("\nTest $(dirpath)")
        # Find all files named test_*.jl in the directory and include them.
        files = sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(dirpath)))
        if isempty(files)
            println("  No test files found in $(dirpath).")
        else
            for f in files
                filepath = joinpath(dirpath, f)
                @time begin
                    println("  Including $(filepath)")
                    include(filepath)
                end
            end
        end
    end
end
