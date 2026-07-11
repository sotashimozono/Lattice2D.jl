"""
    Lattice2D dual-lattice construction
    ===================================

The (planar) **dual** of a lattice places a vertex at the centre of every
face of the original and joins two dual vertices whenever their faces share
an edge. For the regular tilings this is a fixed map between topologies:

| original      | dual          |
|:--------------|:--------------|
| `Square`      | `Square`  (self-dual) |
| `Triangular`  | `Honeycomb`   |
| `Honeycomb`   | `Triangular`  |
| `Kagome`      | `Dice`        |
| `Dice`        | `Kagome`      |

The map is an involution — the dual of the dual is the original topology
(`dual_topology(dual_topology(T)) === T`).

For every supported pair the number of faces of an `Lx × Ly` sample equals
the number of dual sites of an `Lx × Ly` sample of the dual topology, so
`dual_lattice` keeps the `(Lx, Ly)` extents and the boundary condition of
the original. See issue #69.
"""

# ---- topology-level duality -----------------------------------------

"""
    dual_topology(::Type{T}) where {T <: AbstractTopology} -> Type

Return the topology type dual to `T`. Defined for `Square` (self-dual),
`Triangular`/`Honeycomb`, and `Kagome`/`Dice`. Throws `ArgumentError` for
topologies whose dual is not itself a supported tiling (`Lieb`,
`ShastrySutherland`, `UnionJack`).

The map is an involution: `dual_topology(dual_topology(T)) === T`.
"""
dual_topology(::Type{Square}) = Square
dual_topology(::Type{Triangular}) = Honeycomb
dual_topology(::Type{Honeycomb}) = Triangular
dual_topology(::Type{Kagome}) = Dice
dual_topology(::Type{Dice}) = Kagome
function dual_topology(::Type{T}) where {T<:AbstractTopology}
    return throw(
        ArgumentError(
            "dual_topology is not defined for topology $(nameof(T)); " *
            "supported: Square, Triangular, Honeycomb, Kagome, Dice",
        ),
    )
end

"""
    dual_topology(t::AbstractTopology) -> AbstractTopology

Instance form: returns the dual topology singleton, e.g.
`dual_topology(Triangular()) === Honeycomb()`.
"""
dual_topology(::T) where {T<:AbstractTopology} = dual_topology(T)()

# Recover the topology *type* of a Lattice from its type parameter.
_topology_type(::Lattice{Topo}) where {Topo} = Topo

# ---- lattice-level dual ---------------------------------------------

"""
    dual_lattice(lat::Lattice; kwargs...) -> Lattice

Construct the dual lattice of `lat` (see [`dual_topology`](@ref) for the
topology map). The dual keeps the `(Lx, Ly)` extents and the boundary
condition of `lat`; the indexing strategy is preserved unless overridden.
Any extra `kwargs` are forwarded to [`build_lattice`](@ref).

Because the topology map is an involution, `dual_lattice(dual_lattice(lat))`
has the same topology, extents and boundary as `lat`.

# Example
```julia
tri  = triangular(6, 6)
hexy = dual_lattice(tri)          # honeycomb, 6x6, same boundary
topology(dual_lattice(hexy))      # back to the triangular topology
```

Throws `ArgumentError` if the dual of `lat`'s topology is not a supported
tiling (`Lieb`, `ShastrySutherland`, `UnionJack`).

See also [`dual_topology`](@ref).
"""
function dual_lattice(
    lat::Lattice; boundary=LatticeCore.boundary(lat), indexing=lat.indexing, kwargs...
)
    dual_topo = dual_topology(_topology_type(lat))
    return build_lattice(
        dual_topo, lat.Lx, lat.Ly; boundary=boundary, indexing=indexing, kwargs...
    )
end
