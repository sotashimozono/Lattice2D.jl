using Lattice2D
using Lattice2D: CellSite, CellBond
using StaticArrays
using Test

const _TC_TOPOS = (
    Square, Triangular, Honeycomb, Kagome, Lieb, ShastrySutherland, Dice, UnionJack
)

# Real-space displacement vectors from a finite interior site of
# sublattice `b`, via the *independently implemented* `neighbor_bonds`
# path. A 7×7 OBC sample guarantees a full-coordination interior site
# for every sublattice (motif offsets span at most ±1 cell), and OBC
# means the stored displacement is unwrapped.
function _finite_neighbor_vecs(Topo, b)
    lat = build_lattice(Topo, 7, 7; boundary=OpenAxis())
    best_i, best_c = 0, -1
    for i in 1:num_sites(lat)
        sublattice(lat, i) == b || continue
        c = length(neighbor_bonds(lat, i))
        c > best_c && ((best_i, best_c) = (i, c))
    end
    vecs = Set(round.(nb.vector; digits=6) for nb in neighbor_bonds(lat, best_i))
    return vecs, best_c
end

# Same displacement set, but reconstructed purely from the motif via the
# lazy accessors on the infinite lattice.
function _motif_neighbor_vecs(Topo, b)
    inf = InfiniteLattice(Topo)
    p_b = basis_position(inf, b)
    ib = incident_cell_bonds(inf, CellSite((0, 0), b))
    vecs = Set(
        round.(cell_position(inf, CellSite(Tuple(cb.offset), cb.dst)) .- p_b; digits=6) for
        cb in ib
    )
    return vecs
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

@testset "motif neighbours match the finite neighbor_bonds path" begin
    # The strong independent check: the coordination shell reconstructed
    # from the unit-cell motif must equal the one the established finite
    # neighbour code produces for an interior site.
    for Topo in _TC_TOPOS
        for b in 1:num_sublattices(build_lattice(Topo, 2, 2))
            finite_vecs, finite_c = _finite_neighbor_vecs(Topo, b)
            motif_vecs = _motif_neighbor_vecs(Topo, b)
            @test length(motif_vecs) == finite_c
            @test motif_vecs == finite_vecs
        end
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

        fin = build_lattice(Topo, 3, 3)
        @test collect(site_orbits(inf)) == collect(site_orbits(fin))
        @test translation_vectors(inf) == basis_vectors(fin)
        # Same motif as the finite lattice (src, dst, offset, type).
        key(cb) = (cb.src, cb.dst, cb.offset, cb.type)
        @test Set(key.(cell_bonds(inf))) == Set(key.(cell_bonds(fin)))
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
