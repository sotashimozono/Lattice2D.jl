using Plots

# `runtests.jl` already sets `ENV["GKSwstype"] = "100"` so GR runs
# headlessly. We only build the figure objects and inspect their
# series -- nothing is written to disk -- but the env var still
# applies in case any kwarg path triggers a backend init.

@testset "Lattice2DPlotsExt -- brillouin_zone" begin
    @testset "Extension is loaded once Plots is in scope" begin
        ext = Base.get_extension(Lattice2D, :Lattice2DPlotsExt)
        @test ext !== nothing
    end

    @testset "brillouin_zone throws under any open axis" begin
        lat_obc = build_lattice(Square, 4, 4; boundary=OpenAxis())
        @test_throws ArgumentError brillouin_zone(lat_obc)

        cyl = build_lattice(
            Square, 4, 4; boundary=LatticeBoundary((PeriodicAxis(), OpenAxis()))
        )
        @test_throws ArgumentError brillouin_zone(cyl)
    end

    @testset "Square BZ is a square of half-diagonal sqrt(2)*pi centred at Gamma" begin
        lat = build_lattice(Square, 4, 4)
        verts = brillouin_zone(lat)
        @test length(verts) == 4
        # All vertices at distance sqrt(2)*pi from Gamma.
        for v in verts
            @test norm(v) ≈ sqrt(2) * π atol = 1e-10
        end
        # Opposite-vertex pairs cancel through the origin.
        c = sum(verts) / length(verts)
        @test norm(c) <= 1e-10
        # Vertices walk CCW: cross-products of consecutive edges are
        # all positive.
        for k in 1:length(verts)
            v_a = verts[k]
            v_b = verts[mod1(k + 1, length(verts))]
            v_c = verts[mod1(k + 2, length(verts))]
            cross =
                (v_b[1] - v_a[1]) * (v_c[2] - v_b[2]) -
                (v_b[2] - v_a[2]) * (v_c[1] - v_b[1])
            @test cross > 0
        end
    end

    @testset "Triangular BZ is a regular hexagon" begin
        lat = build_lattice(Triangular, 4, 4)
        verts = brillouin_zone(lat)
        @test length(verts) == 6
        # Equal radii (regular hexagon).
        rs = [norm(v) for v in verts]
        @test maximum(rs) - minimum(rs) <= 1e-10
        # Edge lengths equal too.
        edges = [norm(verts[mod1(k + 1, 6)] - verts[k]) for k in 1:6]
        @test maximum(edges) - minimum(edges) <= 1e-10
    end

    @testset "Honeycomb BZ is a hexagon (same Bravais as Triangular)" begin
        # Honeycomb shares the triangular Bravais primitives (the
        # 2-site basis lives inside the unit cell, not in the basis
        # itself), so its BZ has the same hexagonal shape.
        lat_h = build_lattice(Honeycomb, 4, 4)
        verts_h = brillouin_zone(lat_h)
        @test length(verts_h) == 6
        rs = [norm(v) for v in verts_h]
        @test maximum(rs) - minimum(rs) <= 1e-10
    end

    @testset "high_symmetry_points reports the expected sets" begin
        @test haskey(high_symmetry_points(build_lattice(Square, 4, 4)), :X)
        @test haskey(high_symmetry_points(build_lattice(Triangular, 4, 4)), :K)
        @test haskey(high_symmetry_points(build_lattice(Honeycomb, 4, 4)), :K)
        # Unmapped topology falls back to the singleton Gamma.
        hsp_kg = high_symmetry_points(build_lattice(Kagome, 4, 4))
        @test collect(keys(hsp_kg)) == [:Gamma]
    end

    @testset "high_symmetry_points: Square X is at b1/2" begin
        lat = build_lattice(Square, 4, 4)
        hsp = high_symmetry_points(lat)
        ml = reciprocal_lattice(lat)
        B = reciprocal_basis(ml)
        @test hsp[:X] ≈ B[:, 1] / 2 atol = 1e-10
        @test hsp[:M] ≈ (B[:, 1] + B[:, 2]) / 2 atol = 1e-10
        @test hsp[:Gamma] ≈ zeros(2) atol = 1e-10
    end

    @testset "plot_brillouin_zone returns a Plots.Plot for every Bravais topology" begin
        for Topo in (Square, Triangular, Honeycomb)
            lat = build_lattice(Topo, 4, 4)
            p = plot_brillouin_zone(lat)
            @test p isa Plots.Plot
            # `BZ` polygon series is registered.
            labels = [String(s[:label]) for s in p.series_list]
            @test "BZ" in labels
        end
    end

    @testset "plot_brillouin_zone with show_mesh overlays the mesh series" begin
        lat = build_lattice(Square, 6, 6)
        p = plot_brillouin_zone(lat; show_mesh=true)
        labels = [String(s[:label]) for s in p.series_list]
        @test "mesh" in labels
    end

    @testset "plot_brillouin_zone with explicit ml overrides the default" begin
        lat = build_lattice(Square, 4, 4)
        ml = gamma_centered(reciprocal_basis(reciprocal_lattice(lat)), (8, 8))
        p = plot_brillouin_zone(lat; show_mesh=true, ml=ml)
        labels = [String(s[:label]) for s in p.series_list]
        @test "mesh" in labels
    end

    @testset "plot_brillouin_zone with show_high_symmetry adds high-sym points" begin
        lat = build_lattice(Triangular, 4, 4)
        p = plot_brillouin_zone(lat; show_high_symmetry=true)
        labels = [String(s[:label]) for s in p.series_list]
        @test "high-sym" in labels
    end
end
