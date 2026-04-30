ENV["GKSwstype"] = "100"

using Lattice2D, Test
using LinearAlgebra
using Random
using StaticArrays
using Aqua

const FIG_BASE = joinpath(pkgdir(Lattice2D), "docs", "src", "assets", "figures")
const FIG_LAT = joinpath(FIG_BASE, "lattice")
const PATHS = Dict(:geometry => joinpath(FIG_LAT, "geometry"))
mkpath.(values(PATHS))

const dirs = ["core", "lattices", "utils", "api", "disorder"]

@testset "tests" begin
    test_args = copy(ARGS)
    println("Passed arguments ARGS = $(test_args) to tests.")
    @time for dir in dirs
        dirpath = joinpath(@__DIR__, dir)
        println("\nTest $(dirpath)")
        # Find all files named test_*.jl in the directory and include them.
        files = sort(
            filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(dirpath))
        )
        if isempty(files)
            println("  No test files found in $(dirpath).")
            @test true
        else
            for f in files
                filepath = joinpath(dirpath, f)
                println("  Including $(filepath)")
                include(filepath)
            end
        end
    end
end

@testset "Aqua" begin
    # Restrict ambiguity detection to Lattice2D itself so that method
    # ambiguities living in upstream packages (e.g. LatticeCore, Plots)
    # do not turn the suite red.
    #
    # `Plots` is currently a regular `[deps]` entry for backward
    # compatibility (see README "Plots dependency"); plotting is in
    # practice provided through `LatticeCorePlotsExt`. Migrating it to
    # `[weakdeps]` is tracked separately and will be a breaking release,
    # so `stale_deps` is configured to ignore `Plots`.
    Aqua.test_all(
        Lattice2D;
        ambiguities=(; broken=false, recursive=false),
        piracies=true,
        stale_deps=(; ignore=[:Plots]),
        deps_compat=true,
        project_extras=true,
        unbound_args=true,
        undefined_exports=true,
    )
end
