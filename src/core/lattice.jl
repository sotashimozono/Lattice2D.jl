"""
    Lattice{Topo, T}

Finite 2D lattice obtained by tiling a unit cell (described by
[`get_unit_cell`](@ref) on a topology singleton) over an `Lx × Ly`
sample, with a `LatticeCore.LatticeBoundary` defining per-axis
boundary conditions. Subtype of `LatticeCore.AbstractLattice{2, T}`.

The struct stores only the parameters that define the lattice; all
pure-function-of-input data (`position`, `sublattice`, `neighbors`,
`basis_vectors`, ...) are derived on demand through the LatticeCore
interface (`lattice_coord`, `apply_axis_bc`, `to_real`). This mirrors
the reference `SimpleSquareLattice` in LatticeCore and keeps the type
lightweight.

Boundary-condition-dependent aggregates that appear repeatedly
(`bonds`, `plaquettes`, the bond and plaquette reverse-lookup tables)
are memoised through a lazy `cache` field — populated on first access
by [`_get_cache`](@ref) and then reused by the O(local) element-center
and incidence overrides.

Type parameters

- `Topo<:AbstractTopology{2}` — topology singleton
- `T<:AbstractFloat` — numeric type for positions

The boundary, indexing, and layout are stored as **fields** with
abstract eltypes (`LatticeBoundary`, `AbstractIndexing`,
`AbstractSiteLayout`). Hot paths that depend on the concrete types of
these fields (`_connection_steps`, `_neighbors_by_shell`,
`_materialise_plaquettes`) use a function-barrier pattern: they
extract the fields and immediately forward to a specialised kernel,
so the JIT specialises on the runtime types and the inner loops stay
type-stable. This trades 5 type parameters for 2 to cut dispatch /
compile-time cost without measurable hot-path regression. See issue
#48 for context.

Prefer constructing via [`build_lattice`](@ref).
"""
struct Lattice{Topo<:AbstractTopology{2},T<:AbstractFloat} <: AbstractLattice{2,T}
    Lx::Int
    Ly::Int
    boundary::LatticeBoundary
    indexing::AbstractIndexing
    layout::AbstractSiteLayout
    # Lazy cache for bonds / plaquettes / reverse-lookup tables.
    # Stored as Ref{Any} because `LatticeCache{T}` isn't defined at
    # this include point; the getter `_get_cache(lat)` re-asserts the
    # concrete type. See `core/cache.jl`.
    cache::Base.RefValue{Any}
end

# ---- Internal helpers -----------------------------------------------

@inline _dims(lat::Lattice) = (lat.Lx, lat.Ly)

@inline function _basis_sv(::Type{<:Lattice{Topo,T}}) where {Topo,T}
    uc = get_unit_cell(Topo)
    a1 = SVector{2,T}(T(uc.basis[1][1]), T(uc.basis[1][2]))
    a2 = SVector{2,T}(T(uc.basis[2][1]), T(uc.basis[2][2]))
    return a1, a2
end

@inline function _sub_offsets(::Type{<:Lattice{Topo,T}}) where {Topo,T}
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
# `i`, apply per-axis BCs, and yield a `Step{T}` for every valid step.
# This is the one place where the unit-cell topology meets the sample
# boundary.

"""
    Step{T}

Internal, concretely-typed record returned by [`_connection_steps`](@ref):
the result of resolving a unit-cell `Connection` against the sample
boundary. Fields:

- `j::Int` — neighbour site index
- `d_vec::SVector{2,T}` — wrapped displacement vector
- `type::Symbol` — bond type tag (`:type_N`)

Replaces the previous `Tuple{Int,SVector{2,T},Symbol}` representation
to make the eltype concrete and remove the per-step dynamic
`Symbol("type_", conn.type)` from the hot path.
"""
struct Step{T<:AbstractFloat}
    j::Int
    d_vec::SVector{2,T}
    type::Symbol
end

# Cached `connection-type-id → Symbol` lookup for a topology. Built
# once per `Topo` and reused on every `_connection_steps` call so the
# hot path never calls `Symbol("type_", n)` again.
const _CONN_TYPE_SYMBOLS = IdDict{DataType,Vector{Symbol}}()

@inline function _conn_type_symbols(::Type{Topo}) where {Topo<:AbstractTopology}
    cached = get(_CONN_TYPE_SYMBOLS, Topo, nothing)
    cached === nothing || return cached
    uc = get_unit_cell(Topo)
    maxtype = 0
    for conn in uc.connections
        conn.type > maxtype && (maxtype = conn.type)
    end
    syms = Vector{Symbol}(undef, maxtype)
    @inbounds for k in 1:maxtype
        syms[k] = Symbol("type_", k)
    end
    _CONN_TYPE_SYMBOLS[Topo] = syms
    return syms
end

function _connection_steps(lat::Lattice{Topo,T}, i::Int) where {Topo,T}
    # Function barrier: forward to a kernel that captures the runtime
    # types of `boundary.axes` / `indexing` so the inner loops type-
    # specialise even though those fields are abstract on `Lattice`
    # (issue #48 — Lattice has 2 type params; B/I/L are field-only).
    bx, by = lat.boundary.axes
    return _connection_steps_kernel(
        Topo, T, lat.Lx, lat.Ly, num_sublattices(lat), lat.indexing, bx, by, i
    )
end

function _connection_steps_kernel(
    ::Type{Topo}, ::Type{T}, Lx::Int, Ly::Int, nsub::Int, indexing, bx, by, i::Int
) where {Topo,T}
    coord = lattice_coord(indexing, (Lx, Ly), nsub, i)
    cx, cy = coord.cell
    s = coord.sublattice

    uc = get_unit_cell(Topo)
    a1, a2 = _basis_sv(Lattice{Topo,T})
    subs = _sub_offsets(Lattice{Topo,T})
    type_syms = _conn_type_symbols(Topo)

    steps = Step{T}[]
    for conn in uc.connections
        type_sym = @inbounds type_syms[conn.type]
        # Outgoing: this cell's `conn.src_sub` → neighbour cell's `conn.dst_sub`.
        if conn.src_sub == s
            tx, ty = cx + conn.dx, cy + conn.dy
            nx, ok_x = apply_axis_bc(bx, tx, Lx)
            if ok_x
                ny, ok_y = apply_axis_bc(by, ty, Ly)
                if ok_y
                    j = site_index(
                        indexing, (Lx, Ly), nsub, LatticeCoord{2}((nx, ny), conn.dst_sub)
                    )
                    d_vec =
                        conn.dx * a1 +
                        conn.dy * a2 +
                        (subs[conn.dst_sub] - subs[conn.src_sub])
                    push!(steps, Step{T}(j, d_vec, type_sym))
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
                        indexing, (Lx, Ly), nsub, LatticeCoord{2}((nx, ny), conn.src_sub)
                    )
                    d_vec =
                        -(conn.dx * a1 + conn.dy * a2) +
                        (subs[conn.src_sub] - subs[conn.dst_sub])
                    push!(steps, Step{T}(j, d_vec, type_sym))
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
        for step in _connection_steps(lat, i)
            j = step.j
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
    bx, by = lat.boundary.axes
    return _neighbors_by_shell_kernel(
        Topo, T, lat.Lx, lat.Ly, num_sublattices(lat), lat.indexing, bx, by, i, k
    )
end

function _neighbors_by_shell_kernel(
    ::Type{Topo}, ::Type{T}, Lx::Int, Ly::Int, nsub::Int, indexing, bx, by, i::Int, k::Int
) where {Topo,T}
    coord_i = lattice_coord(indexing, (Lx, Ly), nsub, i)
    cx, cy = coord_i.cell
    si = coord_i.sublattice

    a1, a2 = _basis_sv(Lattice{Topo,T})
    subs = _sub_offsets(Lattice{Topo,T})

    # `dist_map[j]` is the minimum distance at which we've reached
    # site j from site i across all unwrapped (dx, dy, sub') candidates.
    dist_map = Dict{Int,T}()

    for dy in (-(Ly - 1)):(Ly - 1), dx in (-(Lx - 1)):(Lx - 1), s in 1:nsub
        tx = cx + dx
        ty = cy + dy
        nx, ok_x = apply_axis_bc(bx, tx, Lx)
        ok_x || continue
        ny, ok_y = apply_axis_bc(by, ty, Ly)
        ok_y || continue
        j = site_index(indexing, (Lx, Ly), nsub, LatticeCoord{2}((nx, ny), s))
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
    for step in _connection_steps(lat, i)
        j = step.j
        if j != i && !(j in seen)
            push!(out, Bond{2,T}(i, j, step.d_vec, step.type))
            push!(seen, j)
        end
    end
    return out
end

# `bonds(lat)` and `plaquettes(lat)` for `Lattice` live in
# `core/element_api.jl` where they pull from the lazy
# [`LatticeCache`](@ref) built by [`_get_cache`](@ref). The
# underlying materialiser for plaquettes stays here because the
# cache builder needs it in its include-order position.

function _materialise_plaquettes(lat::Lattice{Topo,T}) where {Topo,T}
    bx, by = lat.boundary.axes
    return _materialise_plaquettes_kernel(
        Topo, T, lat.Lx, lat.Ly, num_sublattices(lat), lat.indexing, bx, by
    )
end

function _materialise_plaquettes_kernel(
    ::Type{Topo}, ::Type{T}, Lx::Int, Ly::Int, nsub::Int, indexing, bx, by
) where {Topo,T}
    uc = get_unit_cell(Topo)
    a1, a2 = _basis_sv(Lattice{Topo,T})
    subs = _sub_offsets(Lattice{Topo,T})

    out = Plaquette{2,T}[]
    isempty(uc.plaquettes) && return out

    for cy in 1:Ly, cx in 1:Lx
        for rule in uc.plaquettes
            vertices = Int[]
            positions_acc = SVector{2,T}[]
            valid = true
            for (s, dx, dy) in rule.corners
                tx, ty = cx + dx, cy + dy
                nx, ok_x = apply_axis_bc(bx, tx, Lx)
                if !ok_x
                    valid = false
                    break
                end
                ny, ok_y = apply_axis_bc(by, ty, Ly)
                if !ok_y
                    valid = false
                    break
                end
                j = site_index(indexing, (Lx, Ly), nsub, LatticeCoord{2}((nx, ny), s))
                push!(vertices, j)
                # Centroid uses the **unwrapped** corner position
                # (relative to the anchor cell) so the centre of a
                # PBC-wrapped plaquette sits just outside the sample
                # edge instead of inside its interior — same rule
                # `bond_center` uses.
                anchor = (cx - 1) * a1 + (cy - 1) * a2
                corner_pos = anchor + dx * a1 + dy * a2 + subs[s]
                push!(positions_acc, corner_pos)
            end
            if valid
                center = sum(positions_acc) / length(positions_acc)
                push!(out, Plaquette{2,T}(vertices, center, rule.type))
            end
        end
    end
    return out
end

# `num_elements(::Lattice, ::PlaquetteCenter)` and the rest of the
# element-center / incidence API for `Lattice` live in
# `core/element_api.jl`, powered by the cached plaquette vector.

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
