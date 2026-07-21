using Lattice2D
using LatticeCore
using Test

# The point of `blocking` is that it needs no site list. The check that it is *correct* is that it
# agrees with `cell_partition`, which does: on a finite lattice both must classify every bond the
# same way as internal to a coarse cell or crossing between two.

siteij(b) = (getfield(b, 1), getfield(b, 2))

"internal/crossing counts per coarse cell, computed the finite way"
function via_cell_partition(lat, f)
    cells = cell_partition(lat, 1)
    cellof = Dict{Int,Int}()
    for (c, g) in enumerate(cells), i in g
        cellof[i] = c
    end
    nint = 0
    ncross = 0
    for b in bonds(lat)
        i, j = siteij(b)
        cellof[i] == cellof[j] ? (nint += 1) : (ncross += 1)
    end
    return length(cells), nint, ncross
end

"the same counts predicted from the blocked UnitCell alone"
function via_blocking(Topo, Lx, Ly, f)
    bu = blocking(get_unit_cell(Topo), f)
    ncell = (Lx ÷ f) * (Ly ÷ f)
    nint = count(c -> c.dx == 0 && c.dy == 0, bu.connections) * ncell
    ncross = count(c -> !(c.dx == 0 && c.dy == 0), bu.connections) * ncell
    return ncell, nint, ncross
end

@testset "blocking needs no lattice" begin
    # dispatches on a TYPE — nothing is built, no size exists
    bu = blocking(Square, 2)
    @test bu isa UnitCell
    @test length(bu.sublattice_positions) == 4
    @test bu.basis == [[2.0, 0.0], [0.0, 2.0]]
    @test blocking(get_unit_cell(Square), 2).connections == bu.connections
    @test_throws ArgumentError blocking(Square, 1)
    @test_throws ArgumentError blocking(Square, 0)
end

@testset "square: the worked example" begin
    bu = blocking(Square, 2)
    # the dx=+1 connection, from each of the four block offsets
    dxc = filter(c -> c.type == 1, bu.connections)
    got = Set((c.src_sub, c.dst_sub, c.dx, c.dy) for c in dxc)
    @test got == Set([(1, 2, 0, 0), (2, 1, 1, 0), (3, 4, 0, 0), (4, 3, 1, 0)])
    # half the bonds of each type become internal
    @test count(c -> c.dx == 0 && c.dy == 0, bu.connections) == 4
    @test length(bu.connections) == 8
end

@testset "sublattice count and basis scale correctly" begin
    for T in AVAILABLE_LATTICES, f in (2, 3)
        uc = get_unit_cell(T)
        nsub = length(uc.sublattice_positions)
        bu = blocking(uc, f)
        @test length(bu.sublattice_positions) == nsub * f^2
        @test bu.basis == [f .* a for a in uc.basis]
        # bond count per cell is preserved: f² cells' worth of the old connections
        @test length(bu.connections) == length(uc.connections) * f^2
        # every new sublattice index is in range and every one is used exactly once as a source
        srcs = sort([c.src_sub for c in bu.connections])
        @test all(1 .<= srcs .<= nsub * f^2)
    end
end

@testset "blocking agrees with cell_partition — the layering is consistent" begin
    for T in AVAILABLE_LATTICES, (Lx, Ly) in ((4, 4), (8, 4), (4, 8))
        lat = build_lattice(T, Lx, Ly)
        a = via_cell_partition(lat, 2)
        b = via_blocking(T, Lx, Ly, 2)
        @test a == b
    end
end

@testset "idempotence: blocking twice by 2 matches blocking once by 4" begin
    for T in AVAILABLE_LATTICES
        uc = get_unit_cell(T)
        twice = blocking(blocking(uc, 2), 2)
        once = blocking(uc, 4)
        @test length(twice.sublattice_positions) == length(once.sublattice_positions)
        @test twice.basis == once.basis
        @test length(twice.connections) == length(once.connections)
        # and both classify the same number of bonds as internal
        int2 = count(c -> c.dx == 0 && c.dy == 0, twice.connections)
        int4 = count(c -> c.dx == 0 && c.dy == 0, once.connections)
        @test int2 == int4
    end
end

@testset "plaquette rules are dropped, not silently carried" begin
    # a plaquette of the fine lattice is not a plaquette of the coarse one
    uc = get_unit_cell(Square)
    @test !isempty(uc.plaquettes)
    @test isempty(blocking(uc, 2).plaquettes)
end
