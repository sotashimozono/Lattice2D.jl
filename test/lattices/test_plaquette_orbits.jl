using Lattice2D
using StaticArrays
using Test

const _PO_TOPOS = (
    Square, Triangular, Honeycomb, Kagome, Lieb, ShastrySutherland, Dice, UnionJack
)

@testset "plaquette_orbits reproduce the unit cell's plaquette rules" begin
    for Topo in _PO_TOPOS
        rules = get_unit_cell(Topo).plaquettes
        fin = build_lattice(Topo, 4, 4)
        inf = InfiniteLattice(Topo)

        for lat in (fin, inf)
            po = collect(plaquette_orbits(lat))
            @test length(po) == length(rules)
            @test [(p.corners, p.type) for p in po] == [(r.corners, r.type) for r in rules]
            # Routed through the centring trait too (compared by value —
            # PlaquetteRule has no value `==`).
            eo = collect(element_orbits(lat, PlaquetteCenter()))
            @test [(p.corners, p.type) for p in eo] == [(p.corners, p.type) for p in po]
        end
    end
end

@testset "plaquette-centre position is the centroid of its corners" begin
    for Topo in _PO_TOPOS
        inf = InfiniteLattice(Topo)
        for pr in plaquette_orbits(inf)
            # Independent reconstruction of the centroid from the corner
            # (sublattice, dx, dy) specs via cell_position.
            corners = [cell_position(inf, (dx, dy), sub) for (sub, dx, dy) in pr.corners]
            expected = sum(corners) / length(corners)
            @test element_orbit_position(inf, PlaquetteCenter(), pr) ≈ expected
        end
    end
    # Concrete value: the square plaquette is centred at (0.5, 0.5).
    sq = only(plaquette_orbits(InfiniteLattice(Square)))
    @test element_orbit_position(InfiniteLattice(Square), PlaquetteCenter(), sq) ≈
        SVector(0.5, 0.5)
end

@testset "finite and infinite agree on the plaquette centring" begin
    for Topo in _PO_TOPOS
        fin = build_lattice(Topo, 4, 4)
        inf = InfiniteLattice(Topo)
        pf = collect(plaquette_orbits(fin))
        pi = collect(plaquette_orbits(inf))
        @test [(p.corners, p.type) for p in pf] == [(p.corners, p.type) for p in pi]
        for (a, b) in zip(pf, pi)
            @test element_orbit_position(fin, PlaquetteCenter(), a) ≈
                element_orbit_position(inf, PlaquetteCenter(), b)
        end
    end
end

@testset "all three centrings are reachable on the infinite lattice" begin
    inf = InfiniteLattice(Kagome)
    @test length(collect(element_orbits(inf, VertexCenter()))) == 3      # sublattices
    @test length(collect(element_orbits(inf, BondCenter()))) ==
        length(get_unit_cell(Kagome).connections)
    @test length(collect(element_orbits(inf, PlaquetteCenter()))) ==
        length(get_unit_cell(Kagome).plaquettes)
end
