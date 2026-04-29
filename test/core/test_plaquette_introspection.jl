@testset "plaquette introspection helpers" begin
    @testset "get_plaquette_rules ↔ get_unit_cell.plaquettes" begin
        for Topo in (
            Square, Triangular, Honeycomb, Kagome, Lieb, ShastrySutherland, Dice, UnionJack
        )
            rules = get_plaquette_rules(Topo)
            uc_rules = get_unit_cell(Topo).plaquettes
            @test rules isa Vector{PlaquetteRule}
            @test length(rules) == length(uc_rules)
            for (a, b) in zip(rules, uc_rules)
                @test a.corners == b.corners
                @test a.type == b.type
            end
        end
    end

    @testset "kwarg UnitCell constructor matches positional" begin
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        conns = Connection[Connection(1, 1, 1, 0, 1), Connection(1, 1, 0, 1, 1)]
        plaqs = PlaquetteRule[PlaquetteRule(
            [(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)], :square
        )]

        uc_pos = UnitCell{2,Float64}([a1, a2], [[0.0, 0.0]], conns, plaqs)
        uc_kw = UnitCell{2,Float64}(;
            basis=[a1, a2],
            sublattice_positions=[[0.0, 0.0]],
            connections=conns,
            plaquettes=plaqs,
        )
        @test uc_pos.basis == uc_kw.basis
        @test uc_pos.sublattice_positions == uc_kw.sublattice_positions
        @test length(uc_pos.connections) == length(uc_kw.connections)
        @test length(uc_pos.plaquettes) == length(uc_kw.plaquettes)
        @test uc_pos.plaquettes[1].type == uc_kw.plaquettes[1].type
        @test uc_pos.plaquettes[1].corners == uc_kw.plaquettes[1].corners

        # plaquettes default
        uc_kw_default = UnitCell{2,Float64}(;
            basis=[a1, a2], sublattice_positions=[[0.0, 0.0]], connections=conns
        )
        @test uc_kw_default.plaquettes == PlaquetteRule[]
    end

    @testset "num_bonds / num_elements consistency" begin
        for Topo in (Square, Honeycomb, Kagome)
            lat = build_lattice(Topo, 4, 4)
            @test num_bonds(lat) == num_elements(lat, BondCenter())
            @test num_plaquettes(lat) == num_elements(lat, PlaquetteCenter())
            @test length(get_plaquette_rules(Topo)) ≥ 1
        end
    end
end
