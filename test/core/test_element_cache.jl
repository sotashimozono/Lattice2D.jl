@testset "Lattice cache + O(local) element_api" begin
    @testset "cached bonds / plaquettes are stable across repeated calls" begin
        lat = build_lattice(Square, 5, 5)
        # Fresh cache on first access.
        b1 = bonds(lat)
        b2 = bonds(lat)
        @test b1 === b2                                  # same Vector object
        @test b1 == collect(
            b for i in 1:num_sites(lat) for b in neighbor_bonds(lat, i) if b.j > b.i
        )

        p1 = plaquettes(lat)
        p2 = plaquettes(lat)
        @test p1 === p2
    end

    @testset "num_elements via cache agrees with declared plaquette counts" begin
        for (Topo, per_cell) in (
            (Square, 1),
            (Triangular, 2),
            (Honeycomb, 1),
            (Kagome, 3),
            (Lieb, 1),
            (ShastrySutherland, 4),
            (UnionJack, 4),
        )
            lat = build_lattice(Topo, 4, 4)
            @test num_elements(lat, PlaquetteCenter()) == per_cell * 16
        end

        # Dice still has no rules → empty cache bucket → 0.
        @test num_elements(build_lattice(Dice, 3, 3), PlaquetteCenter()) == 0
    end

    @testset "element_position(BondCenter) O(1) matches bond_center" begin
        lat = build_lattice(Triangular, 4, 4)
        bs = collect(bonds(lat))
        for i in 1:length(bs)
            @test element_position(lat, BondCenter(), i) == bond_center(lat, bs[i])
        end
    end

    @testset "element_position(PlaquetteCenter) O(1) matches plaquette.center" begin
        lat = build_lattice(Honeycomb, 3, 3)
        ps = collect(plaquettes(lat))
        for i in 1:length(ps)
            @test element_position(lat, PlaquetteCenter(), i) == ps[i].center
        end
    end

    @testset "element_neighbors(BondCenter) = bonds sharing a vertex" begin
        lat = build_lattice(Square, 4, 4)
        bs = collect(bonds(lat))
        for i in 1:length(bs)
            nbrs = element_neighbors(lat, BondCenter(), i)
            b_self = bs[i]
            for k in nbrs
                other = bs[k]
                @test other.i == b_self.i ||
                    other.i == b_self.j ||
                    other.j == b_self.i ||
                    other.j == b_self.j
                @test k != i
            end
            # Symmetry: j ∈ neighbors(i) ⇒ i ∈ neighbors(j).
            for k in nbrs
                @test i in element_neighbors(lat, BondCenter(), k)
            end
        end
    end

    @testset "element_neighbors(PlaquetteCenter) = dual graph" begin
        # On a PBC Square sample each unit square shares an edge with
        # exactly 4 neighbours (left, right, up, down plaquette).
        lat = build_lattice(Square, 4, 4)
        @test num_elements(lat, PlaquetteCenter()) == 16
        for i in 1:16
            nbrs = element_neighbors(lat, PlaquetteCenter(), i)
            @test length(nbrs) == 4
            # Symmetry.
            for k in nbrs
                @test i in element_neighbors(lat, PlaquetteCenter(), k)
            end
        end
    end

    @testset "incident: bond ↔ plaquette round-trip on every topology" begin
        for Topo in (Square, Triangular, Honeycomb, Kagome, Lieb, UnionJack)
            lat = build_lattice(Topo, 4, 4)
            np = num_elements(lat, PlaquetteCenter())
            np == 0 && continue
            for p_idx in 1:np
                bond_ids = incident(lat, PlaquetteCenter(), BondCenter(), p_idx)
                @test !isempty(bond_ids)
                # Every bond on p's boundary lists p as one of its
                # bounding plaquettes.
                for b in bond_ids
                    @test p_idx in incident(lat, BondCenter(), PlaquetteCenter(), b)
                end
            end
        end
    end

    @testset "incident: vertex ↔ bond round-trip" begin
        lat = build_lattice(Kagome, 3, 3)
        N = num_sites(lat)
        for i in 1:N
            bond_ids = incident(lat, VertexCenter(), BondCenter(), i)
            for b in bond_ids
                @test i in incident(lat, BondCenter(), VertexCenter(), b)
            end
        end
    end

    @testset "incident: vertex ↔ plaquette round-trip" begin
        lat = build_lattice(Honeycomb, 4, 4)
        N = num_sites(lat)
        for i in 1:N
            plaq_ids = incident(lat, VertexCenter(), PlaquetteCenter(), i)
            for p in plaq_ids
                @test i in incident(lat, PlaquetteCenter(), VertexCenter(), p)
            end
        end
    end

    @testset "incident: same-centring falls through to element_neighbors" begin
        lat = build_lattice(Square, 4, 4)
        for i in 1:num_elements(lat, BondCenter())
            @test incident(lat, BondCenter(), BondCenter(), i) ==
                element_neighbors(lat, BondCenter(), i)
        end
        for i in 1:num_sites(lat)
            @test incident(lat, VertexCenter(), VertexCenter(), i) ==
                element_neighbors(lat, VertexCenter(), i)
        end
    end

    @testset "cache survives repeated calls without re-materialising" begin
        # Side-effect check: after the first call, subsequent calls
        # must return *the same* Vector object (===), not a fresh
        # materialisation.
        lat = build_lattice(Kagome, 5, 5)
        b1 = bonds(lat)
        p1 = plaquettes(lat)
        @test b1 === bonds(lat)
        @test p1 === plaquettes(lat)
    end

    @testset "distinct Lattice instances hold independent caches" begin
        lat_a = build_lattice(Square, 4, 4)
        lat_b = build_lattice(Square, 5, 5)
        @test bonds(lat_a) !== bonds(lat_b)
        @test num_elements(lat_a, BondCenter()) != num_elements(lat_b, BondCenter())
    end
end
