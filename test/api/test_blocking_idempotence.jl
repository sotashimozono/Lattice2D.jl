using Lattice2D
using LatticeCore
using Test

# The first idempotence test compared only COUNTS (number of sublattices, number of connections,
# number of internal bonds) and the basis. Two unit cells can match on all of those and describe
# completely different lattices, so that test could not fail for the reason it claimed to check.
#
# The actual claim is that `blocking(blocking(uc, 2), 2)` and `blocking(uc, 4)` describe THE SAME
# LATTICE — identical up to a relabelling of sublattices. So: match the sublattices by their real
# positions, then require the connection sets to agree under that relabelling.
#
# The unit cells are also perturbed away from their library values. A wrong index convention in
# `block_sublattice` can still produce the right counts on a symmetric cell; it cannot survive a
# position-matched comparison on a cell whose sublattices sit at generic points.

const TOL = 1e-9
key(p) = Tuple(round.(p; digits=9))

"""
Relabelling from `a`'s sublattice indices to `b`'s, matched by position.
Returns `nothing` if the position multisets differ at all.
"""
function position_map(a::UnitCell, b::UnitCell)
    length(a.sublattice_positions) == length(b.sublattice_positions) || return nothing
    lookup = Dict{Any,Int}()
    for (i, p) in enumerate(b.sublattice_positions)
        haskey(lookup, key(p)) && return nothing        # degenerate positions: matching is ambiguous
        lookup[key(p)] = i
    end
    m = Vector{Int}(undef, length(a.sublattice_positions))
    for (i, p) in enumerate(a.sublattice_positions)
        j = get(lookup, key(p), 0)
        j == 0 && return nothing
        m[i] = j
    end
    return m
end

function conn_set(uc::UnitCell, relabel=nothing)
    return Set(
        (
            relabel === nothing ? c.src_sub : relabel[c.src_sub],
            relabel === nothing ? c.dst_sub : relabel[c.dst_sub],
            c.dx,
            c.dy,
            c.type,
        ) for c in uc.connections
    )
end

"do two unit cells describe the same lattice?"
function same_lattice(a::UnitCell, b::UnitCell)
    all(isapprox.(reduce(vcat, a.basis), reduce(vcat, b.basis); atol=TOL)) ||
        return false, "basis"
    m = position_map(a, b)
    m === nothing && return false, "sublattice positions"
    conn_set(a, m) == conn_set(b) || return false, "connections"
    return true, ""
end

"""
Perturb a unit cell: shift every sublattice by a distinct generic vector, so that no accidental
symmetry can make a wrong index convention look right.
"""
function perturb(uc::UnitCell{2,T}, seed) where {T}
    pos = [
        p .+ [0.017 * (seed + 3k), 0.023 * (seed + 5k)] for
        (k, p) in enumerate(uc.sublattice_positions)
    ]
    return UnitCell{2,T}(uc.basis, pos, uc.connections, uc.plaquettes)
end

@testset "blocking(·,2)∘blocking(·,2) describes the same lattice as blocking(·,4)" begin
    for T in AVAILABLE_LATTICES
        uc = get_unit_cell(T)
        ok, why = same_lattice(blocking(blocking(uc, 2), 2), blocking(uc, 4))
        @test ok || "$T: $why" == ""          # surface which part differed
    end
end

@testset "… and still with generic sublattice positions" begin
    for T in AVAILABLE_LATTICES, seed in 1:3
        uc = perturb(get_unit_cell(T), seed)
        ok, why = same_lattice(blocking(blocking(uc, 2), 2), blocking(uc, 4))
        @test ok || "$T seed=$seed: $why" == ""
    end
end

@testset "blocking(·,2)∘blocking(·,3) matches blocking(·,6), both orders" begin
    # a non-square composite, and asymmetric so the two orders are a real check
    for T in AVAILABLE_LATTICES
        uc = perturb(get_unit_cell(T), 1)
        a, _ = same_lattice(blocking(blocking(uc, 2), 3), blocking(uc, 6))
        b, _ = same_lattice(blocking(blocking(uc, 3), 2), blocking(uc, 6))
        @test a
        @test b
    end
end

@testset "the position map is a genuine bijection, not a coincidence of counts" begin
    # guard the guard: `position_map` must actually be injective and total
    for T in AVAILABLE_LATTICES
        uc = perturb(get_unit_cell(T), 2)
        m = position_map(blocking(blocking(uc, 2), 2), blocking(uc, 4))
        @test m !== nothing
        @test sort(m) == collect(1:length(m))
    end
end

@testset "the check can actually fail" begin
    # A test that never fails proves nothing. Corrupt one connection and confirm `same_lattice`
    # reports it — otherwise the comparisons above are vacuous.
    uc = get_unit_cell(Square)
    good = blocking(uc, 4)
    bad = UnitCell{2,Float64}(
        good.basis,
        good.sublattice_positions,
        [
            Connection(c.src_sub, c.dst_sub, c.dx + (i == 1), c.dy, c.type) for
            (i, c) in enumerate(good.connections)
        ],
        good.plaquettes,
    )
    @test !first(same_lattice(blocking(blocking(uc, 2), 2), bad))
    # and a shifted sublattice must be caught too
    shifted = UnitCell{2,Float64}(
        good.basis,
        [k == 1 ? p .+ [0.5, 0.0] : p for (k, p) in enumerate(good.sublattice_positions)],
        good.connections,
        good.plaquettes,
    )
    @test !first(same_lattice(blocking(blocking(uc, 2), 2), shifted))
end

@testset "blocking(uc,4) agrees with cell_partition at k=2" begin
    # PRECONDITION, learned the hard way: the finite/infinite correspondence only holds when periodic
    # wrapping does not identify distinct coarse cells. On a 4×4 lattice coarsened by 4 there is a
    # SINGLE coarse cell, so every bond that leaves a cell in the infinite lattice wraps back into the
    # same one and counts as internal — `blocking` and `cell_partition` then legitimately disagree.
    # The comparison is meaningful only when every coarse dimension is at least 2.
    for T in AVAILABLE_LATTICES, (Lx, Ly) in ((8, 8), (16, 8))
        cdims = (Lx ÷ 4, Ly ÷ 4)
        @assert all(cdims .>= 2) "test setup: coarse dims $cdims would let the lattice wrap onto itself"
        lat = build_lattice(T, Lx, Ly)
        cells = cell_partition(lat, 2)
        cellof = Dict{Int,Int}()
        for (c, g) in enumerate(cells), i in g
            cellof[i] = c
        end
        nint = count(b -> cellof[getfield(b, 1)] == cellof[getfield(b, 2)], bonds(lat))
        ncross = length(bonds(lat)) - nint

        bu = blocking(get_unit_cell(T), 4)
        ncell = prod(cdims)
        @test length(cells) == ncell
        @test count(c -> c.dx == 0 && c.dy == 0, bu.connections) * ncell == nint
        @test count(c -> !(c.dx == 0 && c.dy == 0), bu.connections) * ncell == ncross
    end
end

@testset "and the correspondence really does break when the lattice is too small" begin
    # Pin the precondition itself: on a lattice that coarsens to one cell, cell_partition must call
    # every bond internal while blocking still reports crossing connections. If this ever stops
    # failing, the precondition above has silently become unnecessary and the comment is stale.
    lat = build_lattice(Square, 4, 4)
    cells = cell_partition(lat, 2)
    @test length(cells) == 1
    bu = blocking(get_unit_cell(Square), 4)
    @test count(c -> !(c.dx == 0 && c.dy == 0), bu.connections) > 0   # infinite view: bonds leave
    cellof = Dict(i => 1 for i in 1:num_sites(lat))
    @test count(b -> cellof[getfield(b, 1)] != cellof[getfield(b, 2)], bonds(lat)) == 0  # finite: none
end
