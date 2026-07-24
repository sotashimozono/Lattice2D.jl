# Translation-cell (unit-cell motif) layer for Lattice2D.
#
# Implements LatticeCore's translation-cell interface
# (`translation_vectors`, `basis_position`, `cell_bonds`) on the finite
# `Lattice{Topo,T}` by reading the topology's `get_unit_cell`, and adds
# `InfiniteLattice{Topo,T}` — the thermodynamic-limit counterpart. With
# those four methods the whole 2D catalogue gains LatticeCore's lazy
# orbit accessors (`site_orbits`, `bond_orbits`, `cell_position`,
# `neighbors_at`, `incident_cell_bonds`) for free, so a tensor network
# can be placed one tensor per site/bond orbit and contracted along the
# motif without materialising the lattice.
#
# The motif is exactly the topology's `UnitCell`: `sublattice_positions`
# are the basis-site offsets and each `Connection(src, dst, dx, dy,
# type)` is one undirected motif bond. Finiteness and boundary
# conditions only change how that motif is tiled, not the motif itself —
# so `Lattice` (finite) and `InfiniteLattice` share the same geometry
# source and a downstream builder written against the orbit accessors
# runs unchanged across both.

# ---- Shared: UnitCell.Connection -> LatticeCore.CellBond -------------
#
# `build_lattice` tags emitted bonds `:type_N` from the integer
# `Connection.type`; we mirror that convention here so a motif bond and
# its instantiated `Bond` carry the same type symbol.
#
# Cached `Topo -> Tuple{CellBond{2}...}` and reused by every
# `cell_bonds` / `incident_cell_bonds` / `neighbors_at` call, mirroring
# the `_SUB_OFFSETS_CACHE` / `_CONN_TYPE_SYMBOLS` discipline in
# `lattice.jl` — these accessors are meant for hot tensor-network build
# loops. An immutable tuple is returned so the shared value cannot be
# mutated by a caller.
const _CELL_BONDS_CACHE = IdDict{DataType,Any}()

function _cell_bonds_from_unitcell(::Type{Topo}) where {Topo<:AbstractTopology{2}}
    cached = get(_CELL_BONDS_CACHE, Topo, nothing)
    cached === nothing || return cached::Tuple{Vararg{CellBond{2}}}
    conns = get_unit_cell(Topo).connections
    val = Tuple(
        CellBond(c.src_sub, c.dst_sub, (c.dx, c.dy), Symbol("type_", c.type)) for c in conns
    )
    _CELL_BONDS_CACHE[Topo] = val
    return val
end

# ---- Finite Lattice: motif methods ----------------------------------
#
# `num_basis_sites` is not overridden: LatticeCore's default delegates
# to `num_sublattices`, which `Lattice` already implements.

LatticeCore.translation_vectors(lat::Lattice) = basis_vectors(lat)

function LatticeCore.basis_position(lat::Lattice{Topo,T}, b::Int) where {Topo,T}
    offs = _sub_offsets(Lattice{Topo,T})
    1 <= b <= length(offs) ||
        throw(ArgumentError("basis site $b out of range 1:$(length(offs)) for $(Topo)"))
    return offs[b]
end

LatticeCore.cell_bonds(::Lattice{Topo}) where {Topo} = _cell_bonds_from_unitcell(Topo)

# ---- InfiniteLattice: the thermodynamic-limit object ----------------

"""
    InfiniteLattice{Topo, T, L}(; layout = UniformLayout(IsingSite()))

A truly infinite 2D lattice of topology `Topo` — the
thermodynamic-limit counterpart of the finite [`Lattice`](@ref) under
periodic boundary conditions. It carries no size (`size_trait` is
`InfiniteSize`, [`num_sites`](@ref) throws) and is described by, and
accessed through, its unit-cell motif:

- `site_orbits` / `bond_orbits` — the finite fundamental domain (one
  entry per sublattice / per `Connection` of `get_unit_cell(Topo)`);
- `cell_position`, `neighbors_at`, `incident_cell_bonds` — on-demand
  access to any [`CellSite`](@ref), computed without materialising the
  lattice.

If a finite sample is needed, [`materialize`](@ref) tiles the motif
into a periodic [`Lattice`](@ref):

```julia
inf = InfiniteLattice(Honeycomb)
fin = materialize(inf; dims = (8, 8))   # 8×8 PBC honeycomb Lattice
```

but the lazy accessors above never require it.

Positions use `Float64`, matching the finite [`Lattice`](@ref) produced
by [`build_lattice`](@ref).

# Example
```julia
inf = InfiniteLattice(Kagome)
site_orbits(inf)                         # 1:3  (three sublattices)
length(collect(bond_orbits(inf)))        # motif bond count
neighbors_at(inf, CellSite((0, 0), 1))   # neighbours of basis-1 site
```
"""
struct InfiniteLattice{Topo<:AbstractTopology{2},T<:AbstractFloat,L<:AbstractSiteLayout} <:
       AbstractLattice{2,T}
    layout::L
end

function InfiniteLattice(
    ::Type{Topo}; layout::AbstractSiteLayout=UniformLayout(IsingSite())
) where {Topo<:AbstractTopology{2}}
    return InfiniteLattice{Topo,Float64,typeof(layout)}(layout)
end

# ---- Motif interface (same geometry source as the finite Lattice) ---

function LatticeCore.translation_vectors(::InfiniteLattice{Topo,T}) where {Topo,T}
    a1, a2 = get_unit_cell(Topo).basis
    return SMatrix{2,2,T}(T(a1[1]), T(a1[2]), T(a2[1]), T(a2[2]))
end

function LatticeCore.num_sublattices(::InfiniteLattice{Topo}) where {Topo}
    return length(get_unit_cell(Topo).sublattice_positions)
end

function LatticeCore.basis_position(::InfiniteLattice{Topo,T}, b::Int) where {Topo,T}
    subs = get_unit_cell(Topo).sublattice_positions
    1 <= b <= length(subs) ||
        throw(ArgumentError("basis site $b out of range 1:$(length(subs)) for $(Topo)"))
    p = subs[b]
    return SVector{2,T}(T(p[1]), T(p[2]))
end

function LatticeCore.cell_bonds(::InfiniteLattice{Topo}) where {Topo}
    return _cell_bonds_from_unitcell(Topo)
end

# ---- Size / traits --------------------------------------------------

LatticeCore.size_trait(::InfiniteLattice) = InfiniteSize()
LatticeCore.periodicity(::InfiniteLattice) = Periodic()
LatticeCore.reciprocal_support(::InfiniteLattice) = HasReciprocal()
LatticeCore.site_layout(l::InfiniteLattice) = l.layout

function LatticeCore.topology(::InfiniteLattice{Topo}) where {Topo}
    return TopologyTrait{topology_name(Topo())}()
end

# Conceptually the thermodynamic limit of full periodic boundaries.
LatticeCore.boundary(::InfiniteLattice) = LatticeBoundary((PeriodicAxis(), PeriodicAxis()))

# `is_bipartite` is a topology-intrinsic property, so the infinite
# lattice must agree with the finite one rather than fall through to the
# generic `false` default. Reuse the finite BFS 2-colouring on an open
# 6×6 patch: every motif offset spans at most ±1 cell, so the patch
# contains the shortest cycle, and OBC avoids the boundary-wrap odd
# cycles that a periodic sample could introduce — matching the infinite
# graph.
function LatticeCore.is_bipartite(::InfiniteLattice{Topo}) where {Topo}
    return is_bipartite(build_lattice(Topo, 6, 6; boundary=OpenAxis()))
end

# There is no linear site index, so `sublattice(inf, i)` is undefined.
# Throw the same way as `position` / `neighbors` / `num_sites` rather
# than return the misleading generic default of `1`; a site's sublattice
# is the `.basis` field of its `CellSite`.
function LatticeCore.sublattice(::InfiniteLattice, ::Int)
    return throw(
        DomainError(
            InfiniteSize(),
            "InfiniteLattice has no linear site index; a site's sublattice is the " *
            "`.basis` field of its `CellSite`",
        ),
    )
end

# ---- Linear-index API is undefined for an infinite lattice ----------

function LatticeCore.num_sites(::InfiniteLattice)
    return throw(
        DomainError(
            InfiniteSize(),
            "InfiniteLattice has no finite site count; use `site_orbits` for the " *
            "unit-cell basis or `materialize(lat; dims)` for a finite sample",
        ),
    )
end

function LatticeCore.position(::InfiniteLattice, ::Int)
    return throw(
        DomainError(
            InfiniteSize(),
            "InfiniteLattice has no linear site index; address sites by `CellSite` " *
            "and use `cell_position(lat, site)`",
        ),
    )
end

function LatticeCore.neighbors(::InfiniteLattice, ::Int)
    return throw(
        DomainError(
            InfiniteSize(),
            "InfiniteLattice has no linear site index; use `neighbors_at(lat, ::CellSite)`",
        ),
    )
end

# ---- Bridge: materialise into a finite periodic sample --------------

"""
    materialize(lat::InfiniteLattice{Topo}; dims::NTuple{2, Int}) → Lattice{Topo}

Tile the motif into a `dims[1] × dims[2]` periodic [`Lattice`](@ref) of
the same topology, carrying the infinite lattice's site layout. This is
the optional infinite → finite bridge; the lazy translation-cell
accessors do not require it.
"""
function LatticeCore.materialize(
    lat::InfiniteLattice{Topo}; dims::NTuple{2,Int}
) where {Topo}
    return build_lattice(Topo, dims[1], dims[2]; boundary=PeriodicAxis(), layout=lat.layout)
end
