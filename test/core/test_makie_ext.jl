using CairoMakie          # loads Makie -> triggers Lattice2DMakieExt
using LatticeCore
using Test

# `runtests.jl` sets ENV["GKSwstype"] = "100"; CairoMakie renders headlessly
# to disk regardless. We build the figures, render each to a temp PNG (a real
# end-to-end check that the draw calls are valid), and check the returned type.
const MAKIE_EXT = Base.get_extension(Lattice2D, :Lattice2DMakieExt)

renders(fig) =
    mktemp() do path, io
        close(io)
        p = path * ".png"
        CairoMakie.save(p, fig)
        ok = isfile(p) && filesize(p) > 0
        rm(p; force=true)
        return ok
    end

@testset "Lattice2DMakieExt" begin
    @testset "extension loads with Makie" begin
        @test MAKIE_EXT !== nothing
    end

    @testset "makie_lattice returns a renderable Figure" begin
        for lat in (build_lattice(Square, 4, 4), build_lattice(Honeycomb, 3, 3))
            fig = makie_lattice(lat)
            @test fig isa Makie.Figure
            @test renders(fig)
        end
        # sublattice colouring + bond highlight paths
        lat = build_lattice(Honeycomb, 3, 3)
        fig = makie_lattice(lat; colorby=:sublattice, highlight_bonds=[1, 2, 3])
        @test fig isa Makie.Figure
        @test renders(fig)
        # single-colour path
        @test makie_lattice(lat; colorby=:none, show_sites=true) isa Makie.Figure
    end

    @testset "makie_state returns a renderable Figure" begin
        lat = build_lattice(Square, 5, 5)
        N = num_sites(lat)
        state = Float64.(1:N)
        fig = makie_state(lat, state)
        @test fig isa Makie.Figure
        @test renders(fig)
        # XY-angle arrows path
        angles = [2π * i / N for i in 1:N]
        @test makie_state(lat, angles; arrows=true) isa Makie.Figure
        # length mismatch is rejected
        @test_throws DimensionMismatch makie_state(lat, Float64.(1:(N - 1)))
    end

    @testset "makie_structure_factor: figure + S(k) numerics" begin
        lat = build_lattice(Square, 6, 6)
        N = num_sites(lat)
        # smoke: builds and renders at low resolution
        fig = makie_structure_factor(lat, ones(N); resolution=16)
        @test fig isa Makie.Figure
        @test renders(fig)

        # Independent numerics via the internal grid helper. Odd resolution puts
        # k = 0 at the centre; a uniform state must light up the Γ peak:
        #   S(0) = |Σ_j 1|² / N = N.
        R = 21
        ks, S = MAKIE_EXT._structure_factor_grid(
            lat, ones(N); k_range=(-π, π), resolution=R
        )
        c = (R + 1) ÷ 2
        @test abs(ks[c]) < 1e-12                      # centre is k = 0
        @test S[c, c] ≈ float(N) atol = 1e-8
        # S(k) = S(-k) for a real state (inversion symmetry of |·|²)
        @test S ≈ reverse(S) rtol = 1e-10
        # S ≥ 0 everywhere
        @test all(S .≥ -1e-12)

        # validation
        @test_throws DimensionMismatch makie_structure_factor(lat, ones(N - 1))
        @test_throws ArgumentError makie_structure_factor(lat, ones(N); resolution=0)
    end
end
