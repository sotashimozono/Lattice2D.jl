"""
    AbstractTopology{D}

Abstract supertype for 2D lattice topologies (Square, Triangular,
Honeycomb, ...). Each concrete subtype is a singleton that acts as
a dispatch key for [`get_unit_cell`](@ref), and, through
`TopologyTrait` in LatticeCore, as the `topology(lat)` value of the
resulting [`PeriodicLattice2D`](@ref).
"""
abstract type AbstractTopology{D} end

"""
    Connection(src_sub, dst_sub, dx, dy, type)

A unit-cell-level **connection rule** — the description of an edge
between sublattices inside a single unit cell, or between a
sublattice in one cell and a sublattice in a neighbouring cell.

This is distinct from `LatticeCore.Bond`:

- A `Connection` is a **template** on the unit cell (`src_sub`,
  `dst_sub` are sublattice ids, `dx`, `dy` are *relative* cell
  offsets). It is topology data, known statically from
  [`get_unit_cell`](@ref).
- A `LatticeCore.Bond` is an **instantiated** edge on a concrete
  `Lx × Ly` sample, with absolute site indices and a wrapped
  displacement vector. It is the per-sample output of `build_lattice`.

Fields
- `src_sub::Int` — 1-based sublattice id of the source site
- `dst_sub::Int` — 1-based sublattice id of the destination site
- `dx::Int` — x-axis cell offset (0 = same unit cell)
- `dy::Int` — y-axis cell offset
- `type::Int` — bond type tag. Currently stored as an `Int` on the
  `Connection` side for backward compatibility; `build_lattice`
  converts this to a `Symbol` (`:type_N`) when it emits
  `LatticeCore.Bond`.
"""
struct Connection
    src_sub::Int
    dst_sub::Int
    dx::Int
    dy::Int
    type::Int
end

"""
    UnitCell{D, T}

Static topology description for a `D`-dimensional Bravais lattice
with an arbitrary basis. Contains

- `basis::Vector{Vector{T}}` — the `D` primitive vectors
- `sublattice_positions::Vector{Vector{T}}` — one offset per
  geometric sublattice inside the unit cell
- `connections::Vector{Connection}` — the full list of
  intra- and inter-cell connection rules for this topology
- `plaquettes::Vector{LatticeCore.PlaquetteRule}` — declarative list
  of plaquette rules (one per plaquette **kind**) anchored at the
  reference cell. Empty for topologies that haven't been wired up to
  the plaquette API yet.

Produced by [`get_unit_cell`](@ref) on an `AbstractTopology`
singleton.
"""
struct UnitCell{D,T}
    basis::Vector{Vector{T}}
    sublattice_positions::Vector{Vector{T}}
    connections::Vector{Connection}
    plaquettes::Vector{PlaquetteRule}
end

# Back-compat 3-arg constructor — legacy unit cells without any
# declared plaquettes still work, they just contribute 0 plaquettes
# to their lattice.
function UnitCell{D,T}(
    basis::Vector{Vector{T}},
    sublattice_positions::Vector{Vector{T}},
    connections::Vector{Connection},
) where {D,T}
    return UnitCell{D,T}(basis, sublattice_positions, connections, PlaquetteRule[])
end

"""
    get_unit_cell(::Type{T}) where {T <: AbstractTopology}

Return the [`UnitCell`](@ref) describing topology `T`. Concrete
topology types (`Square`, `Triangular`, ...) specialise this
method. The fallback throws.
"""
function get_unit_cell(::Type{T}) where {T<:AbstractTopology}
    return error("UnitCell not defined for $T")
end
