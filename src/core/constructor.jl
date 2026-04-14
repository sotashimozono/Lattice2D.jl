"""
    _resolve_boundary(axis::LatticeCore.AbstractAxisBC)
    _resolve_boundary(boundary::LatticeCore.LatticeBoundary)

Normalise the user-facing `boundary` argument of [`build_lattice`](@ref)
to a `LatticeBoundary{2}`. A single axis BC is broadcast to both axes
with a `NoModifier`; an explicit `LatticeBoundary` is returned as-is.
"""
_resolve_boundary(axis::AbstractAxisBC) = LatticeBoundary((axis, axis), NoModifier())
_resolve_boundary(boundary::LatticeBoundary) = boundary

"""
    build_lattice(Topology, Lx, Ly;
                  boundary = PeriodicAxis(),
                  indexing = RowMajor(),
                  layout   = UniformLayout(IsingSite())) → Lattice

Construct a finite 2D lattice of topology `Topology` on an
`Lx × Ly` sample. The `boundary` argument accepts either a single
`LatticeCore.AbstractAxisBC` (broadcast to both axes) or an
explicit `LatticeCore.LatticeBoundary` for mixed-axis setups such
as cylinders.

All geometric / connectivity data (positions, neighbours, bonds) are
computed lazily by the accessors defined on [`Lattice`](@ref);
this constructor only validates and wraps its arguments, so building
is O(1) regardless of `Lx × Ly`.

# Examples

```julia
# 4 × 4 periodic square lattice, row-major linearisation, Ising sites
lat = build_lattice(Square, 4, 4)

# Open honeycomb 6 × 6
open_hc = build_lattice(Honeycomb, 6, 6; boundary = OpenAxis())

# Cylinder: square lattice with PBC in x, OBC in y
cyl = build_lattice(Square, 4, 4;
    boundary = LatticeBoundary((PeriodicAxis(), OpenAxis())))

# Kagome with XY sites
xy_kg = build_lattice(Kagome, 4, 4; layout = UniformLayout(XYSite()))
```
"""
function build_lattice(
    Topology::Type{<:AbstractTopology{2}},
    Lx::Int,
    Ly::Int;
    boundary=PeriodicAxis(),
    indexing::AbstractIndexing=RowMajor(),
    layout::AbstractSiteLayout=UniformLayout(IsingSite()),
)
    bc = _resolve_boundary(boundary)
    T = Float64
    return Lattice{Topology,T,typeof(bc),typeof(indexing),typeof(layout)}(
        Lx, Ly, bc, indexing, layout,
    )
end

"""
    Lattice2D.Lattice(Topology, Lx, Ly; kwargs...)

Convenience alias for [`build_lattice`](@ref). Kept so that
`Lattice2D.Lattice(Square, 4, 4)` reads naturally at the call site.
"""
function Lattice(Topology::Type{<:AbstractTopology{2}}, Lx::Int, Ly::Int; kwargs...)
    return build_lattice(Topology, Lx, Ly; kwargs...)
end
