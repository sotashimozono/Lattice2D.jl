@testset "Dice T3 plaquettes + bond_type / num_bonds / num_plaquettes" begin
    @testset "Dice: 3 rhombi per cell (east, northeast, northwest)" begin
        lat = build_lattice(Dice, 4, 4)
        @test num_plaquettes(lat) == 3 * 16
        ps = collect(plaquettes(lat))
        @test count(p -> p.type === :rhombus_east, ps) == 16
        @test count(p -> p.type === :rhombus_northeast, ps) == 16
        @test count(p -> p.type === :rhombus_northwest, ps) == 16
        @test all(length(p.vertices) == 4 for p in ps)
    end

    @testset "Dice rhombus: all 4 edges are the NN distance 1/√3" begin
        lat = build_lattice(Dice, 6, 6; boundary=OpenAxis())
        ps = collect(plaquettes(lat))
        # Pick an interior rhombus of each orientation.
        target = SVector(3.0, 2.5)
        for ty in (:rhombus_east, :rhombus_northeast, :rhombus_northwest)
            p = first(p for p in ps if p.type === ty)
            for i in 1:4
                a = position(lat, p.vertices[i])
                b = position(lat, p.vertices[mod1(i + 1, 4)])
                @test norm(b - a) ≈ 1 / sqrt(3) atol = 1e-10
            end
        end
    end

    @testset "Dice rhombus: 2 hubs + 2 rims, one rim1 and one rim2" begin
        lat = build_lattice(Dice, 4, 4)
        ps = collect(plaquettes(lat))
        for p in ps
            sub_ids = [sublattice(lat, v) for v in p.vertices]
            # Hubs (sub 1) × 2, rim1 (sub 2) × 1, rim2 (sub 3) × 1.
            @test count(==(1), sub_ids) == 2
            @test count(==(2), sub_ids) == 1
            @test count(==(3), sub_ids) == 1
        end
    end

    @testset "num_bonds matches cached bond count" begin
        for (Topo, per_site_half_degree) in (
            (Square, 2),              # 4 NN / 2 = 2 per site
            (Triangular, 3),          # 6 / 2 = 3
            (Honeycomb, 3 / 2),        # 3 / 2 = 1.5 but per cell it's 3 per 2 sites = 3/cell
            (Kagome, 2),               # 4 / 2 = 2 (each kagome site has 4 NN)
        )
            lat = build_lattice(Topo, 4, 4)
            # Cross-check: num_bonds equals length(bonds(lat)) (cached
            # Vector so this is O(1)) and equals
            # num_elements(lat, BondCenter()).
            @test num_bonds(lat) == length(bonds(lat))
            @test num_bonds(lat) == num_elements(lat, BondCenter())
        end
    end

    @testset "num_plaquettes matches cached plaquette count" begin
        for Topo in (Square, Triangular, Honeycomb, Kagome, Lieb, ShastrySutherland, Dice, UnionJack)
            lat = build_lattice(Topo, 4, 4)
            @test num_plaquettes(lat) == length(plaquettes(lat))
            @test num_plaquettes(lat) == num_elements(lat, PlaquetteCenter())
        end
    end

    @testset "bond_type: retrieves the tag inherited from Connection.type" begin
        # Square has Connection.type = 1 (right) and 2 (up), so every
        # bond is either :type_1 or :type_2.
        sq = build_lattice(Square, 4, 4)
        for b in bonds(sq)
            @test bond_type(sq, b.i, b.j) === b.type
            # Symmetric in i, j.
            @test bond_type(sq, b.j, b.i) === b.type
        end
    end

    @testset "bond_type on ShastrySutherland: dimer vs NN" begin
        lat = build_lattice(ShastrySutherland, 3, 3)
        bs = collect(bonds(lat))
        nn_count = count(b -> bond_type(lat, b.i, b.j) === :type_1, bs)
        dimer_count = count(b -> bond_type(lat, b.i, b.j) === :type_2, bs)
        @test nn_count > 0
        @test dimer_count > 0
        @test nn_count + dimer_count == length(bs)
    end

    @testset "bond_type on non-bond throws ArgumentError" begin
        lat = build_lattice(Square, 3, 3)
        # Two non-adjacent sites on a 3×3 PBC square: 1 and 5 (across
        # the diagonal) — they aren't connected by a NN bond.
        @test_throws ArgumentError bond_type(lat, 1, 5)
    end
end
