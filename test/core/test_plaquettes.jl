@testset "plaquette enumeration on Lattice" begin
    @testset "Square: 1 plaquette per unit cell under PBC" begin
        for (Lx, Ly) in ((3, 3), (4, 4), (5, 3))
            lat = build_lattice(Square, Lx, Ly)
            @test num_elements(lat, PlaquetteCenter()) == Lx * Ly
            ps = collect(plaquettes(lat))
            @test length(ps) == Lx * Ly
            @test all(p.type === :square for p in ps)
            @test all(length(p.vertices) == 4 for p in ps)
        end
    end

    @testset "Square: centroid is the cell + (0.5, 0.5)" begin
        lat = build_lattice(Square, 4, 4)
        ps = collect(plaquettes(lat))
        # The plaquette anchored at cell (1, 1) has centroid (0.5, 0.5).
        p1 = first(p for p in ps if p.center ≈ SVector(0.5, 0.5))
        @test p1.type === :square
        @test length(p1.vertices) == 4
    end

    @testset "Triangular: 2 plaquettes per unit cell, unit edges" begin
        lat = build_lattice(Triangular, 6, 6)
        @test num_elements(lat, PlaquetteCenter()) == 2 * 6 * 6
        ps = collect(plaquettes(lat))
        @test count(p -> p.type === :up_triangle, ps) == 6 * 6
        @test count(p -> p.type === :down_triangle, ps) == 6 * 6

        # Every interior up-triangle has unit edges.
        center = SVector(3.0, 3.0)
        up = first(
            p for p in ps if p.type === :up_triangle && norm(p.center - center) < 1.5
        )
        @test length(up.vertices) == 3
        for i in 1:3
            a = position(lat, up.vertices[i])
            b = position(lat, up.vertices[mod1(i + 1, 3)])
            @test norm(b - a) ≈ 1.0 atol = 1e-10
        end
    end

    @testset "Honeycomb: 1 hexagon per cell, 6 vertices, unit edges" begin
        lat = build_lattice(Honeycomb, 6, 6)
        @test num_elements(lat, PlaquetteCenter()) == 6 * 6
        ps = collect(plaquettes(lat))
        @test all(p.type === :hexagon for p in ps)
        @test all(length(p.vertices) == 6 for p in ps)

        # Pick an interior hexagon and check every edge is unit length.
        target = SVector(3.0, 3.0)
        k = argmin(norm(p.center - target) for p in ps)
        p = ps[k]
        n = length(p.vertices)
        for i in 1:n
            a = position(lat, p.vertices[i])
            b = position(lat, p.vertices[mod1(i + 1, n)])
            @test norm(b - a) ≈ 1.0 atol = 1e-10
        end
    end

    @testset "Honeycomb: every site lives on exactly 3 hexagons (PBC)" begin
        lat = build_lattice(Honeycomb, 4, 4)
        ps = collect(plaquettes(lat))
        counts = zeros(Int, num_sites(lat))
        for p in ps
            for v in p.vertices
                counts[v] += 1
            end
        end
        @test all(counts .== 3)
    end

    @testset "Kagome: 2 triangles per cell (hexagon not yet declared)" begin
        lat = build_lattice(Kagome, 6, 6)
        @test num_elements(lat, PlaquetteCenter()) == 2 * 6 * 6
        ps = collect(plaquettes(lat))
        @test count(p -> p.type === :up_triangle, ps) == 6 * 6
        @test count(p -> p.type === :down_triangle, ps) == 6 * 6
        @test all(length(p.vertices) == 3 for p in ps)
    end

    @testset "OBC drops plaquettes whose corners leave the sample" begin
        # 4×4 honeycomb OBC: the hexagon rule uses corner offsets
        # (0,0), (1,-1), (1,0), (0,1), (0,0). Valid anchor cells have
        # cx ∈ [1, Lx-1] (dx=+1 reachable) and cy ∈ [2, Ly-1]
        # (dy=-1 and dy=+1 both reachable), giving (Lx-1)*(Ly-2)
        # surviving plaquettes.
        for (Lx, Ly) in ((4, 4), (5, 4), (4, 5))
            lat = build_lattice(Honeycomb, Lx, Ly; boundary=OpenAxis())
            ps = collect(plaquettes(lat))
            @test length(ps) == (Lx - 1) * (Ly - 2)
            @test num_elements(lat, PlaquetteCenter()) == (Lx - 1) * (Ly - 2)
        end

        # Comparison to PBC: PBC version has Lx*Ly plaquettes.
        @test length(collect(plaquettes(build_lattice(Honeycomb, 4, 4)))) == 16
    end

    @testset "topologies without plaquette rules return 0" begin
        for Topo in (Lieb, ShastrySutherland, Dice, UnionJack)
            lat = build_lattice(Topo, 3, 3)
            @test num_elements(lat, PlaquetteCenter()) == 0
            @test isempty(collect(plaquettes(lat)))
        end
    end

    @testset "PlaquetteCenter delegates to plaquettes via element_position" begin
        lat = build_lattice(Square, 4, 4)
        for i in 1:num_elements(lat, PlaquetteCenter())
            @test element_position(lat, PlaquetteCenter(), i) ==
                collect(plaquettes(lat))[i].center
        end
    end

    @testset "incident: bond ↔ plaquette round trip on Square" begin
        lat = build_lattice(Square, 3, 3)
        bs = collect(bonds(lat))
        # Pick the first plaquette; its 4 boundary bonds must each
        # list that plaquette in their own plaquette-incidence.
        bond_ids = incident(lat, PlaquetteCenter(), BondCenter(), 1)
        @test length(bond_ids) == 4
        for b in bond_ids
            @test 1 in incident(lat, BondCenter(), PlaquetteCenter(), b)
        end
    end
end
