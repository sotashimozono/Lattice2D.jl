"""
    Lattice{Topo, T, B, I, L}

Finite 2D lattice obtained by tiling a unit cell (described by
[`get_unit_cell`](@ref) on a topology singleton) over an `Lx × Ly`
sample, with a `LatticeCore.LatticeBoundary` defining per-axis
boundary conditions. Subtype of `LatticeCore.AbstractLattice{2, T}`.

The struct stores only the parameters that define the lattice; all
geometric and connectivity data (`position`, `sublattice`, `neighbors`,
`bonds`, ...) are derived on demand from the `topology` type parameter
and the `LatticeCore` interface (`lattice_coord`, `apply_axis_bc`,
`to_real`). This mirrors the reference `SimpleSquareLattice` in
LatticeCore and keeps the type lightweight.

Type parameters

- `Topo<:AbstractTopology{2}` — topology singleton
- `T<:AbstractFloat` — numeric type for positions
- `B<:LatticeCore.LatticeBoundary` — composite boundary condition
- `I<:LatticeCore.AbstractIndexing` — linearisation strategy
- `L<:LatticeCore.AbstractSiteLayout` — site-type layout

Prefer constructing via [`build_lattice`](@ref).
"""
struct Lattice{
    Topo<:AbstractTopology{2},
    T<:AbstractFloat,
    B<:LatticeBoundary,
    I<:AbstractIndexing,
    L<:AbstractSiteLayout,
} <: AbstractLattice{2,T}
    Lx::Int
    Ly::Int
    boundary::B
    indexing::I
    layout::L
end

# ---- Internal helpers -----------------------------------------------

@inline _dims(lat::Lattice) = (lat.Lx, lat.Ly)

@inline function _basis_sv(::Type{Lattice{Topo,T,B,I,L}}) where {Topo,T,B,I,L}
    uc = get_unit_cell(Topo)
    a1 = SVector{2,T}(T(uc.basis[1][1]), T(uc.basis[1][2]))
    a2 = SVector{2,T}(T(uc.basis[2][1]), T(uc.basis[2][2]))
    return a1, a2
end

@inline function _sub_offsets(::Type{Lattice{Topo,T,B,I,L}}) where {Topo,T,B,I,L}
    uc = get_unit_cell(Topo)
    return [SVector{2,T}(T(p[1]), T(p[2])) for p in uc.sublattice_positions]
end

# ---- LatticeCore required interface ---------------------------------

function LatticeCore.num_sublattices(::Lattice{Topo}) where {Topo}
    return length(get_unit_cell(Topo).sublattice_positions)
end

LatticeCore.num_sites(lat::Lattice) = lat.Lx * lat.Ly * num_sublattices(lat)

function LatticeCore.position(lat::Lattice, i::Int)
    coord = lattice_coord(lat.indexing, _dims(lat), num_sublattices(lat), i)
    return to_real(lat, coord).x
end

function LatticeCore.sublattice(lat::Lattice, i::Int)
    coord = lattice_coord(lat.indexing, _dims(lat), num_sublattices(lat), i)
    return coord.sublattice
end

LatticeCore.boundary(lat::Lattice) = lat.boundary

LatticeCore.site_layout(lat::Lattice) = lat.layout

LatticeCore.size_trait(lat::Lattice) = FiniteSize((lat.Lx, lat.Ly))

# ---- Neighbour walk -------------------------------------------------
#
# Shared driver used by `neighbors` and `neighbor_bonds`: iterate the
# outgoing + incoming unit-cell connections at the cell containing
# `i`, apply per-axis BCs, and yield `(j, d_vec, type)` for every
# valid step. This is the one place where the unit-cell topology meets
# the sample boundary.

function _connection_steps(lat::Lattice{Topo,T}, i::Int) where {Topo,T}
    Lx, Ly = _dims(lat)
    nsub = num_sublattices(lat)
    coord = lattice_coord(lat.indexing, (Lx, Ly), nsub, i)
    cx, cy = coord.cell
    s = coord.sublattice

    uc = get_unit_cell(Topo)
    bx, by = lat.boundary.axes
    a1, a2 = _basis_sv(typeof(lat))
    subs = _sub_offsets(typeof(lat))

    steps = Tuple{Int,SVector{2,T},Symbol}[]
    for conn in uc.connections
        # Outgoing: this cell's `conn.src_sub` → neighbour cell's `conn.dst_sub`.
        if conn.src_sub == s
            tx, ty = cx + conn.dx, cy + conn.dy
            nx, ok_x = apply_axis_bc(bx, tx, Lx)
            if ok_x
                ny, ok_y = apply_axis_bc(by, ty, Ly)
                if ok_y
                    j = site_index(
                        lat.indexing,
                        (Lx, Ly),
                        nsub,
                        LatticeCoord{2}((nx, ny), conn.dst_sub),
                    )
                    d_vec =
                        conn.dx * a1 +
                        conn.dy * a2 +
                        (subs[conn.dst_sub] - subs[conn.src_sub])
                    push!(steps, (j, d_vec, Symbol("type_", conn.type)))
                end
            end
        end
        # Incoming: neighbour cell's `conn.src_sub` → this cell's `conn.dst_sub`.
        if conn.dst_sub == s
            tx, ty = cx - conn.dx, cy - conn.dy
            nx, ok_x = apply_axis_bc(bx, tx, Lx)
            if ok_x
                ny, ok_y = apply_axis_bc(by, ty, Ly)
                if ok_y
                    j = site_index(
                        lat.indexing,
                        (Lx, Ly),
                        nsub,
                        LatticeCoord{2}((nx, ny), conn.src_sub),
                    )
                    d_vec =
                        -(conn.dx * a1 + conn.dy * a2) +
                        (subs[conn.src_sub] - subs[conn.dst_sub])
                    push!(steps, (j, d_vec, Symbol("type_", conn.type)))
                end
            end
        end
    end
    return steps
end

"""
    neighbors(lat::Lattice, i::Int) → Vector{Int}
    neighbors(lat::Lattice, i::Int; shell::Int) → Vector{Int}

Neighbour indices of site `i`.

- Without the `shell` keyword, returns **all declared unit-cell
  connections**. For most topologies this is the geometric NN; for
  topologies whose unit cell mixes bond types at different distances
  (e.g. [`ShastrySutherland`](@ref), which declares both square NN
  and dimer J′ bonds), the returned set is the union.
- With `shell=k` (k ≥ 1), returns the `k`-th **geometric** neighbour
  shell, ranked by Euclidean distance. For `ShastrySutherland` this
  separates the square NN (`shell=1`) from the dimer partner
  (`shell=2`).

The geometric search is bounded by the sample size — for PBC the
minimum-image distance is used, for OBC only in-sample candidates
are considered.
"""
function LatticeCore.neighbors(lat::Lattice, i::Int; shell::Union{Nothing,Int}=nothing)
    if shell === nothing
        out = Int[]
        seen = Set{Int}()
        for (j, _, _) in _connection_steps(lat, i)
            if j != i && !(j in seen)
                push!(out, j)
                push!(seen, j)
            end
        end
        return out
    else
        shell ≥ 1 || throw(ArgumentError("shell must be ≥ 1, got $shell"))
        return _neighbors_by_shell(lat, i, shell)
    end
end

# Compute geometric k-th shell neighbours by scanning a full
# `(−Lx, Lx) × (−Ly, Ly)` window of unit-cell displacements, filtering
# candidates through the axis BCs and ranking unique distances.
function _neighbors_by_shell(lat::Lattice{Topo,T}, i::Int, k::Int) where {Topo,T}
    Lx, Ly = _dims(lat)
    nsub = num_sublattices(lat)
    coord_i = lattice_coord(lat.indexing, (Lx, Ly), nsub, i)
    cx, cy = coord_i.cell
    si = coord_i.sublattice

    bx, by = lat.boundary.axes
    a1, a2 = _basis_sv(typeof(lat))
    subs = _sub_offsets(typeof(lat))

    # `dist_map[j]` is the minimum distance at which we've reached
    # site j from site i across all unwrapped (dx, dy, sub') candidates.
    dist_map = Dict{Int,T}()

    for dy in -(Ly - 1):(Ly - 1), dx in -(Lx - 1):(Lx - 1), s in 1:nsub
        tx = cx + dx
        ty = cy + dy
        nx, ok_x = apply_axis_bc(bx, tx, Lx)
        ok_x || continue
        ny, ok_y = apply_axis_bc(by, ty, Ly)
        ok_y || continue
        j = site_index(lat.indexing, (Lx, Ly), nsub, LatticeCoord{2}((nx, ny), s))
        j == i && continue
        d_vec = dx * a1 + dy * a2 + (subs[s] - subs[si])
        d = norm(d_vec)
        cur = get(dist_map, j, typemax(T))
        if d < cur
            dist_map[j] = d
        end
    end

    isempty(dist_map) && return Int[]

    # Group unique distances with a tolerance (float noise from basis
    # matrix multiplication can jitter the last few ulps).
    sorted_d = sort(collect(values(dist_map)))
    shells = T[sorted_d[1]]
    for d in sorted_d
        ref = shells[end]
        if d - ref > 10 * eps(T) * max(one(T), ref)
            push!(shells, d)
        end
    end

    k > length(shells) && return Int[]
    target = shells[k]
    tol = 10 * eps(T) * max(one(T), target)
    result = Int[]
    for (j, d) in dist_map
        if abs(d - target) ≤ tol
            push!(result, j)
        end
    end
    return sort!(result)
end

function LatticeCore.neighbor_bonds(lat::Lattice{Topo,T}, i::Int) where {Topo,T}
    out = Bond{2,T}[]
    seen = Set{Int}()
    for (j, d_vec, type) in _connection_steps(lat, i)
        if j != i && !(j in seen)
            push!(out, Bond{2,T}(i, j, d_vec, type))
            push!(seen, j)
        end
    end
    return out
end

function LatticeCore.bonds(lat::Lattice)
    return (b for i in 1:num_sites(lat) for b in neighbor_bonds(lat, i) if b.j > b.i)
end

# ---- Graph / topology traits ----------------------------------------

function LatticeCore.is_bipartite(lat::Lattice)
    N = num_sites(lat)
    colors = zeros(Int, N)
    for i in 1:N
        colors[i] == 0 || continue
        colors[i] = 1
        queue = [i]
        while !isempty(queue)
            u = popfirst!(queue)
            for v in neighbors(lat, u)
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

function LatticeCore.periodicity(lat::Lattice)
    return all(!(ax isa OpenAxis) for ax in lat.boundary.axes) ? Periodic() : Aperiodic()
end

function LatticeCore.reciprocal_support(lat::Lattice)
    return if all(!(ax isa OpenAxis) for ax in lat.boundary.axes)
        HasReciprocal()
    else
        NoReciprocal()
    end
end

function LatticeCore.topology(::Lattice{Topo}) where {Topo}
    return TopologyTrait{topology_name(Topo())}()
end

topology_name(::Square) = :square
topology_name(::Triangular) = :triangular
topology_name(::Honeycomb) = :honeycomb
topology_name(::Kagome) = :kagome
topology_name(::Lieb) = :lieb
topology_name(::ShastrySutherland) = :shastry_sutherland
topology_name(::Dice) = :dice
topology_name(::UnionJack) = :union_jack

# ---- Basis / reciprocal ---------------------------------------------

function LatticeCore.basis_vectors(lat::Lattice{Topo,T}) where {Topo,T}
    a1, a2 = _basis_sv(typeof(lat))
    return SMatrix{2,2,T}(a1[1], a1[2], a2[1], a2[2])
end

function LatticeCore.reciprocal_lattice(lat::Lattice{Topo,T}) where {Topo,T}
    reciprocal_support(lat) isa HasReciprocal || throw(
        ArgumentError("Lattice has no reciprocal lattice unless every axis is periodic")
    )
    A = basis_vectors(lat)
    B = SMatrix{2,2,T}(T(2π) * inv(transpose(A)))
    return monkhorst_pack(B, (lat.Lx, lat.Ly))
end

# ---- Coordinate conversions -----------------------------------------

function LatticeCore.to_real(lat::Lattice{Topo,T}, coord::LatticeCoord{2}) where {Topo,T}
    cx, cy = coord.cell
    s = coord.sublattice
    a1, a2 = _basis_sv(typeof(lat))
    subs = _sub_offsets(typeof(lat))
    return RealSpace{2,T}((cx - 1) * a1 + (cy - 1) * a2 + subs[s])
end
