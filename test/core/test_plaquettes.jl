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

    @testset "Kagome: 3 plaquette kinds (up/down triangle + hexagon)" begin
        lat = build_lattice(Kagome, 6, 6)
        @test num_elements(lat, PlaquetteCenter()) == 3 * 6 * 6
        ps = collect(plaquettes(lat))
        @test count(p -> p.type === :up_triangle, ps) == 6 * 6
        @test count(p -> p.type === :down_triangle, ps) == 6 * 6
        @test count(p -> p.type === :hexagon, ps) == 6 * 6
        @test count(length(p.vertices) == 3 for p in ps) == 2 * 6 * 6
        @test count(length(p.vertices) == 6 for p in ps) == 6 * 6
    end

    @testset "Kagome hexagon: 6 unit edges, interior sample" begin
        lat = build_lattice(Kagome, 8, 8; boundary=OpenAxis())
        ps = collect(plaquettes(lat))
        # Kagome sites are at half the Bravais spacing, so unit-cell
        # NN distance is 0.5 (the hexagon edge length).
        target = SVector(5.0, 3.0)
        hex = first(p for p in ps if p.type === :hexagon && norm(p.center - target) < 1.0)
        for i in 1:6
            a = position(lat, hex.vertices[i])
            b = position(lat, hex.vertices[mod1(i + 1, 6)])
            @test norm(b - a) ≈ 0.5 atol = 1e-10
        end
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

    @testset "Lieb: one 8-vertex square per cell" begin
        lat = build_lattice(Lieb, 4, 4)
        @test num_elements(lat, PlaquetteCenter()) == 4 * 4
        ps = collect(plaquettes(lat))
        @test all(p.type === :square for p in ps)
        @test all(length(p.vertices) == 8 for p in ps)

        # Interior plaquette anchored at cell (1,1) has centroid
        # (1.0, 1.0) and every edge has unit length.
        p = first(p for p in ps if p.center ≈ SVector(1.0, 1.0))
        for i in 1:8
            a = position(lat, p.vertices[i])
            b = position(lat, p.vertices[mod1(i + 1, 8)])
            @test norm(b - a) ≈ 1.0 atol = 1e-10
        end
    end

    @testset "ShastrySutherland: 4 small squares per cell, half dimer" begin
        for (Lx, Ly) in ((3, 3), (4, 4))
            lat = build_lattice(ShastrySutherland, Lx, Ly)
            @test num_elements(lat, PlaquetteCenter()) == 4 * Lx * Ly
            ps = collect(plaquettes(lat))
            @test count(p -> p.type === :dimer_square, ps) == 2 * Lx * Ly
            @test count(p -> p.type === :square, ps) == 2 * Lx * Ly
            @test all(length(p.vertices) == 4 for p in ps)
        end
    end

    @testset "ShastrySutherland: every dimer_square has a J' dimer inside" begin
        lat = build_lattice(ShastrySutherland, 3, 3)
        ps = collect(plaquettes(lat))
        # Dimer bonds are those with type :type_2 in Lattice2D's bond
        # tagging (Connection.type = 2 for the diagonals).
        dimer_edges = Set{Tuple{Int,Int}}()
        for b in bonds(lat)
            if b.type === :type_2
                push!(dimer_edges, (min(b.i, b.j), max(b.i, b.j)))
            end
        end

        for p in ps
            p.type === :dimer_square || continue
            # CCW-ordered 4 corners ⇒ opposite corners are (1,3) and (2,4).
            found = false
            for (i, j) in ((1, 3), (2, 4))
                edge = (
                    min(p.vertices[i], p.vertices[j]), max(p.vertices[i], p.vertices[j])
                )
                if edge in dimer_edges
                    found = true
                    break
                end
            end
            @test found
        end
    end

    @testset "UnionJack: 4 triangle kinds per cell, 3 vertices each" begin
        lat = build_lattice(UnionJack, 4, 4)
        @test num_elements(lat, PlaquetteCenter()) == 4 * 4 * 4
        ps = collect(plaquettes(lat))
        for ty in (:triangle_south, :triangle_east, :triangle_north, :triangle_west)
            @test count(p -> p.type === ty, ps) == 16
        end
        @test all(length(p.vertices) == 3 for p in ps)
    end

    @testset "UnionJack: every triangle has one body + two corners" begin
        lat = build_lattice(UnionJack, 4, 4)
        ps = collect(plaquettes(lat))
        for p in ps
            sub_ids = [sublattice(lat, v) for v in p.vertices]
            @test count(==(1), sub_ids) == 2
            @test count(==(2), sub_ids) == 1
        end
    end

    @testset "Dice still has no plaquette rules declared" begin
        # Dice (T3) rhombus plaquettes are left for a follow-up PR.
        lat = build_lattice(Dice, 3, 3)
        @test num_elements(lat, PlaquetteCenter()) == 0
        @test isempty(collect(plaquettes(lat)))
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
