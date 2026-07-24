using Lattice2D
using Lattice2D: CellSite, CellBond
using StaticArrays
using Test

const _TC_TOPOS = (
    Square, Triangular, Honeycomb, Kagome, Lieb, ShastrySutherland, Dice, UnionJack
)

# (displacement vector, bond-type symbol) pairs from a finite interior
# site of sublattice `b`, via the *independently implemented*
# `neighbor_bonds` path (which derives its `:type_N` tag through the
# separate `_conn_type_symbols` cache, not `_cell_bonds_from_unitcell`).
# A 7×7 OBC sample guarantees a full-coordination interior site for
# every sublattice (motif offsets span at most ±1 cell), and OBC means
# the stored displacement is unwrapped.
function _finite_neighbor_pairs(Topo, b)
    lat = build_lattice(Topo, 7, 7; boundary=OpenAxis())
    best_i, best_c = 0, -1
    for i in 1:num_sites(lat)
        sublattice(lat, i) == b || continue
        c = length(neighbor_bonds(lat, i))
        c > best_c && ((best_i, best_c) = (i, c))
    end
    pairs = Set(
        (round.(nb.vector; digits=6), nb.type) for nb in neighbor_bonds(lat, best_i)
    )
    return pairs, best_c
end

# Same (vector, type) set, reconstructed purely from the motif via the
# lazy accessors on the infinite lattice. Comparing types here — not
# just displacements — catches a mangled `:type_N` mapping that the
# geometry alone would miss.
function _motif_neighbor_pairs(Topo, b)
    inf = InfiniteLattice(Topo)
    p_b = basis_position(inf, b)
    ib = incident_cell_bonds(inf, CellSite((0, 0), b))
    return Set(
        (
            round.(cell_position(inf, CellSite(Tuple(cb.offset), cb.dst)) .- p_b; digits=6),
            cb.type,
        ) for cb in ib
    )
end

@testset "finite Lattice motif interface — all topologies" begin
    for Topo in _TC_TOPOS
        lat = build_lattice(Topo, 4, 4)
        uc = get_unit_cell(Topo)
        nsub = num_sublattices(lat)

        @test translation_vectors(lat) == basis_vectors(lat)
        @test num_basis_sites(lat) == nsub
        @test collect(site_orbits(lat)) == collect(1:nsub)
        @test length(collect(cell_bonds(lat))) == length(uc.connections)

        # basis_position matches the topology's sublattice offsets.
        for b in 1:nsub
            @test basis_position(lat, b) ≈ SVector{2,Float64}(uc.sublattice_positions[b])
        end
        @test_throws ArgumentError basis_position(lat, nsub + 1)
    end
end

@testset "motif neighbours (vector + type) match the finite neighbor_bonds path" begin
    # The strong independent check: the coordination shell reconstructed
    # from the unit-cell motif — displacement AND bond type — must equal
    # the one the established finite neighbour code produces for an
    # interior site.
    for Topo in _TC_TOPOS
        for b in 1:num_sublattices(build_lattice(Topo, 2, 2))
            finite_pairs, finite_c = _finite_neighbor_pairs(Topo, b)
            motif_pairs = _motif_neighbor_pairs(Topo, b)
            @test length(motif_pairs) == finite_c
            @test motif_pairs == finite_pairs
        end
    end
end

@testset "cell_bonds reproduces the raw UnitCell Connections" begin
    # Independent of `_cell_bonds_from_unitcell`: rebuild the expected
    # motif directly from each topology's raw `Connection` fields and
    # diff. Catches a src/dst swap, a wrong offset, or a mangled type
    # index that the geometry-only neighbour check could absorb.
    for Topo in _TC_TOPOS
        conns = get_unit_cell(Topo).connections
        expected = Set(
            (c.src_sub, c.dst_sub, (c.dx, c.dy), Symbol("type_", c.type)) for c in conns
        )
        got = Set(
            (cb.src, cb.dst, Tuple(cb.offset), cb.type) for
            cb in cell_bonds(InfiniteLattice(Topo))
        )
        @test got == expected
    end
end

@testset "InfiniteLattice — traits and undefined linear index" begin
    for Topo in _TC_TOPOS
        inf = InfiniteLattice(Topo)
        @test size_trait(inf) == InfiniteSize()
        @test !is_finite(inf)
        @test periodicity(inf) == Periodic()

        @test_throws DomainError num_sites(inf)
        @test_throws DomainError position(inf, 1)
        @test_throws DomainError neighbors(inf, 1)
        # No linear index ⇒ sublattice(i) must throw, not return 1.
        @test_throws DomainError sublattice(inf, 1)
        # Out-of-range basis fails loudly.
        @test_throws ArgumentError basis_position(inf, num_basis_sites(inf) + 1)

        fin = build_lattice(Topo, 3, 3)
        @test collect(site_orbits(inf)) == collect(site_orbits(fin))
        @test translation_vectors(inf) == basis_vectors(fin)

        # `is_bipartite` is topology-intrinsic: the infinite lattice must
        # agree with a finite (PBC, even) sample — an independent path
        # (OBC 6×6 BFS vs PBC 4×4 BFS).
        @test is_bipartite(inf) == is_bipartite(build_lattice(Topo, 4, 4))
    end
end

@testset "InfiniteLattice — materialize bridge and layout propagation" begin
    inf = InfiniteLattice(Honeycomb)
    fin = materialize(inf; dims=(4, 4))
    @test fin isa Lattice
    @test num_sites(fin) == 4 * 4 * 2      # honeycomb has 2 sublattices
    @test is_finite(fin)
    @test periodicity(fin) == Periodic()

    # Custom layout flows through infinite → materialize → finite.
    inf2 = InfiniteLattice(Square; layout=UniformLayout(XYSite()))
    @test site_type(inf2, 1) == XYSite()
    fin2 = materialize(inf2; dims=(2, 2))
    @test site_type(fin2, 1) == XYSite()
end

@testset "materialize — non-square dims, every topology" begin
    # Non-square dims catch an Lx/Ly transposition in materialize; the
    # loop covers the topologies the Honeycomb/Square cases above miss.
    for Topo in _TC_TOPOS
        nsub = num_basis_sites(InfiniteLattice(Topo))
        fin = materialize(InfiniteLattice(Topo); dims=(3, 5))
        @test fin isa Lattice
        @test num_sites(fin) == 3 * 5 * nsub
        @test periodicity(fin) == Periodic()
    end
end

@testset "lazy accessors also work on the finite Lattice" begin
    # The orbit accessors are generic over AbstractLattice, so they must
    # dispatch on a finite `Lattice` too (not only InfiniteLattice).
    lat = build_lattice(Square, 4, 4)
    nb = neighbors_at(lat, CellSite((2, 2), 1))
    @test length(nb) == 4
    @test Set(n.cell .- SVector(2, 2) for n in nb) ==
        Set(SVector.([(1, 0), (-1, 0), (0, 1), (0, -1)]))
    @test length(collect(incident_cell_bonds(lat, CellSite((2, 2), 1)))) == 4
end

@testset "InfiniteLattice — lazy access is translation invariant" begin
    inf = InfiniteLattice(Kagome)
    base = Set(
        n.cell - c for n in neighbors_at(inf, CellSite((0, 0), 1)) for c in (SVector(0, 0),)
    )
    for cell in ((5, 5), (-3, 8), (100, -100))
        shifted = Set(
            n.cell - SVector(cell...) for n in neighbors_at(inf, CellSite(cell, 1))
        )
        @test shifted == base
        # basis index of a neighbour is a valid sublattice id.
        @test all(
            1 <= n.basis <= num_basis_sites(inf) for
            n in neighbors_at(inf, CellSite(cell, 1))
        )
    end
end
