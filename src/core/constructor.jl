"""
    calc_reciprocal_vectors(basis::Vector{<:Vector}) → Vector{Vector{Float64}}

Return the reciprocal-lattice basis associated with a real-space
`basis` expressed as a vector of vectors. Uses
``B = 2\\pi (A^{-\\top})`` where the columns of ``A`` are the
primitive real-space vectors. Provided for parity with the old
public API; new code should prefer `LatticeCore.reciprocal_lattice(lat)`.
"""
function calc_reciprocal_vectors(basis::AbstractVector{<:AbstractVector})
    A = hcat(basis...)
    B = 2π * inv(A')
    return [B[:, i] for i in 1:size(B, 2)]
end

"""
    check_bipartite_bfs(N::Int, neighbors::Vector{Vector{Int}}) → Bool

BFS two-colouring of an adjacency list. Returns `true` iff the
graph is bipartite.
"""
function check_bipartite_bfs(N::Int, neighbors::Vector{Vector{Int}})
    colors = zeros(Int, N)
    for i in 1:N
        colors[i] == 0 || continue
        colors[i] = 1
        queue = [i]
        while !isempty(queue)
            u = popfirst!(queue)
            for v in neighbors[u]
                if colors[v] == 0
                    colors[v] = -colors[u]
                    push!(queue, v)
                elseif colors[v] == colors[u]
                    return false
                end
            end
        end
    end
    return true
end

# ---- Private helpers ------------------------------------------------
#
# Translate a (cell_x, cell_y, sublattice) triple into a 1-based site
# index using a LatticeCore indexing strategy. All three shipped
# strategies — RowMajor, ColMajor, Snake — reduce to the 2D formulas
# already implemented in LatticeCore, so we just wrap them.

@inline function _coord_to_index(
    indexing::AbstractIndexing, x::Int, y::Int, s::Int, Lx::Int, Ly::Int, nsub::Int
)
    return site_index(indexing, (Lx, Ly), nsub, LatticeCoord{2}((x, y), s))
end

"""
    _resolve_boundary(axis::LatticeCore.AbstractAxisBC)
    _resolve_boundary(boundary::LatticeCore.LatticeBoundary)

Normalise the user-facing `boundary` argument of
[`build_lattice`](@ref) to a `LatticeBoundary{2}`. A single axis BC
is broadcast to both axes with a `NoModifier`; an explicit
`LatticeBoundary` is returned as-is.
"""
_resolve_boundary(axis::AbstractAxisBC) = LatticeBoundary((axis, axis), NoModifier())
_resolve_boundary(boundary::LatticeBoundary) = boundary

# ---- Main constructor -----------------------------------------------

"""
    build_lattice(Topology, Lx, Ly;
                  boundary = PeriodicAxis(),
                  indexing = RowMajor(),
                  layout   = UniformLayout(IsingSite())) → PeriodicLattice2D

Construct a finite 2D lattice of topology `Topology` on an
`Lx × Ly` sample. The `boundary` argument accepts either a single
`LatticeCore.AbstractAxisBC` (broadcast to both axes) or an
explicit `LatticeCore.LatticeBoundary` for mixed-axis setups such
as cylinders.

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
    return _build(Topology, Lx, Ly, bc, indexing, layout)
end

function _build(
    Topology::Type{<:AbstractTopology{2}},
    Lx::Int,
    Ly::Int,
    bc::LatticeBoundary,
    indexing::AbstractIndexing,
    layout::AbstractSiteLayout,
)
    T = Float64
    uc = get_unit_cell(Topology)
    nsub = length(uc.sublattice_positions)
    N = Lx * Ly * nsub

    # --- Geometry: positions + sublattice ids -------------------------
    positions = Vector{SVector{2,T}}(undef, N)
    sublattice_ids = Vector{Int}(undef, N)
    site_map = Matrix{Int}(undef, Lx, Ly)

    a1 = SVector{2,T}(T(uc.basis[1][1]), T(uc.basis[1][2]))
    a2 = SVector{2,T}(T(uc.basis[2][1]), T(uc.basis[2][2]))
    basis_matrix = SMatrix{2,2,T}(a1[1], a1[2], a2[1], a2[2])

    sub_offsets = [SVector{2,T}(T(p[1]), T(p[2])) for p in uc.sublattice_positions]

    for y in 1:Ly, x in 1:Lx
        base_id = _coord_to_index(indexing, x, y, 1, Lx, Ly, nsub)
        site_map[x, y] = base_id
        cell_origin = (x - 1) * a1 + (y - 1) * a2
        for s in 1:nsub
            i = _coord_to_index(indexing, x, y, s, Lx, Ly, nsub)
            positions[i] = cell_origin + sub_offsets[s]
            sublattice_ids[i] = s
        end
    end

    # --- Connectivity: expand each Connection into per-site bonds -----
    bx, by = bc.axes
    bonds = Bond{2,T}[]
    sizehint!(bonds, length(uc.connections) * Lx * Ly)
    nn_table = [Int[] for _ in 1:N]

    for y in 1:Ly, x in 1:Lx
        for conn in uc.connections
            tx = x + conn.dx
            ty = y + conn.dy

            # Apply per-axis BC. If either axis rejects the step,
            # the candidate bond is dropped.
            nx, ok_x = apply_axis_bc(bx, tx, Lx)
            ok_x || continue
            ny, ok_y = apply_axis_bc(by, ty, Ly)
            ok_y || continue

            src = _coord_to_index(indexing, x, y, conn.src_sub, Lx, Ly, nsub)
            dst = _coord_to_index(indexing, nx, ny, conn.dst_sub, Lx, Ly, nsub)

            # The wrapped displacement vector: compute it from the
            # *ideal* (dx, dy) cell offset + intra-cell sublattice
            # displacement, not from the wrapped positions. This gives
            # the true local bond direction regardless of wrapping.
            d_vec =
                conn.dx * a1 +
                conn.dy * a2 +
                (sub_offsets[conn.dst_sub] - sub_offsets[conn.src_sub])

            push!(nn_table[src], dst)
            push!(nn_table[dst], src)
            push!(bonds, Bond{2,T}(src, dst, d_vec, Symbol("type_", conn.type)))
        end
    end

    is_bipartite_flag = check_bipartite_bfs(N, nn_table)

    return PeriodicLattice2D{Topology,T,typeof(bc),typeof(indexing),typeof(layout)}(
        Lx,
        Ly,
        N,
        positions,
        nn_table,
        bonds,
        basis_matrix,
        sublattice_ids,
        is_bipartite_flag,
        site_map,
        Topology(),
        bc,
        indexing,
        layout,
    )
end

"""
    Lattice2D.Lattice(Topology, Lx, Ly; kwargs...)

Alias for [`build_lattice`](@ref). Kept for readability; returns a
[`PeriodicLattice2D`](@ref).
"""
function Lattice(Topology::Type{<:AbstractTopology{2}}, Lx::Int, Ly::Int; kwargs...)
    return build_lattice(Topology, Lx, Ly; kwargs...)
end
