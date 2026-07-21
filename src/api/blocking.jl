# Blocking on the LOCAL description: the scale change that does not need a site list.
#
# `cell_partition` works on site indices and therefore only on a finite lattice. A `UnitCell` plus its
# `Connection`s already describes the lattice completely without any size, so blocking can be defined
# there instead — and then it applies to an infinite lattice as directly as to a finite one.
#
# The arithmetic is one integer division. Blocking by `f` in D dimensions: the new cell holds `f^D`
# old cells, so its sublattices are labelled by (old sublattice, block offset o ∈ {0,…,f−1}^D), and
# the new basis is `f·aᵢ`. A `Connection` with offset `d` leaving block offset `o` lands on
# `target = o + d`, which splits as
#     new block offset = mod(target, f)      (which sublattice of the new cell)
#     new cell offset  = fld(target, f)      (which new cell)
# `fld` rather than `div`, so that negative offsets floor the right way.

"""
    block_sublattice(nsub, f, sub, o) -> Int

Index of the new sublattice carrying old sublattice `sub` at block offset `o` (a `(dx, dy)` in
`0:f-1`). Ordering is old-sublattice-fastest, then `ox`, then `oy`, so `sub = 1, o = (0,0)` is 1.
"""
function block_sublattice(nsub::Int, f::Int, sub::Int, o::NTuple{2,Int})
    return sub + nsub * (o[1] + f * o[2])
end

"""
    blocking(uc::UnitCell{2,T}, f::Integer = 2) -> UnitCell{2,T}

The unit cell of the lattice coarse-grained by a factor `f` per axis: `f²` old cells become one new
cell, whose basis is `f` times the old one and whose sublattices are the old sublattices at each of
the `f²` block offsets.

Every `Connection` is re-expressed in the new cell. A bond whose two ends land in the same new cell
becomes an internal connection (offset `(0,0)`); one that crosses gets the corresponding new offset.
No site indices are involved, so this applies to an infinite lattice exactly as it does to a finite
one — unlike [`cell_partition`](@ref), which enumerates `1:num_sites`.

```jldoctest
julia> using Lattice2D

julia> uc = get_unit_cell(Square);

julia> b = blocking(uc, 2);

julia> length(b.sublattice_positions)      # 2×2 block of a one-site cell
4

julia> b.basis                              # the basis is doubled
2-element Vector{Vector{Float64}}:
 [2.0, 0.0]
 [0.0, 2.0]

julia> count(c -> c.dx == 0 && c.dy == 0, b.connections)   # bonds now internal to the new cell
4
```
"""
function blocking(uc::UnitCell{2,T}, f::Integer=2) where {T}
    f > 1 || throw(ArgumentError("blocking needs f > 1; got $f"))
    fi = Int(f)
    nsub = length(uc.sublattice_positions)

    basis = [fi .* a for a in uc.basis]
    # sublattice positions: the old offsets, plus the block displacement in real space
    pos = Vector{Vector{T}}(undef, nsub * fi^2)
    for oy in 0:(fi - 1), ox in 0:(fi - 1), s in 1:nsub
        disp = ox .* uc.basis[1] .+ oy .* uc.basis[2]
        pos[block_sublattice(nsub, fi, s, (ox, oy))] = uc.sublattice_positions[s] .+ disp
    end

    conns = Connection[]
    for c in uc.connections, oy in 0:(fi - 1), ox in 0:(fi - 1)
        tx, ty = ox + c.dx, oy + c.dy
        src = block_sublattice(nsub, fi, c.src_sub, (ox, oy))
        dst = block_sublattice(nsub, fi, c.dst_sub, (mod(tx, fi), mod(ty, fi)))
        push!(conns, Connection(src, dst, fld(tx, fi), fld(ty, fi), c.type))
    end

    # plaquette rules are stated in old-cell coordinates and are not transported: a plaquette of the
    # fine lattice is not a plaquette of the coarse one. Left empty rather than silently reused.
    return UnitCell{2,T}(basis, pos, conns, PlaquetteRule[])
end

"""
    blocking(::Type{Topo}, f::Integer = 2) -> UnitCell

Convenience form: `blocking(get_unit_cell(Topo), f)`. Takes a topology *type*, so no lattice — and in
particular no lattice size — has to exist first.
"""
blocking(Topo::Type{<:AbstractTopology{2}}, f::Integer=2) = blocking(get_unit_cell(Topo), f)
