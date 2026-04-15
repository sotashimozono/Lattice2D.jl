"""
    LatticeCache{T}

Lazy, per-lattice cache of bonds, plaquettes, and the reverse-lookup
tables needed to make the element-center / incidence API O(local)
instead of O(num_bonds + num_plaquettes) per call.

The cache is constructed once on first access through
[`_get_cache`](@ref) and lives inside a `Ref` field on the parent
[`Lattice`](@ref), so the struct itself stays immutable and cheap to
copy. For sizes typical of MC work (10⁴–10⁶ sites), the cache
occupies a few hundred kilobytes and is negligible relative to the
state vectors MC actually keeps hot.

Fields

- `bonds::Vector{Bond{2, T}}` — materialised bond list, ordered by
  `(cell_y, cell_x, connection_index)` just like the generic
  `bonds(lat)` iteration.
- `plaquettes::Vector{Plaquette{2, T}}` — materialised plaquette
  list, ordered by `(cell_y, cell_x, rule_index)`.
- `bond_index_of::Dict{Tuple{Int, Int}, Int}` — `(min_i, max_j) →
  integer bond index`. Keyed by the canonical (smaller, larger)
  endpoint tuple so that `bond_index_of[(min(i, j), max(i, j))]`
  always hits regardless of which endpoint walked first.
- `plaquettes_by_vertex::Vector{Vector{Int}}` — per-site list of
  plaquette indices whose boundary contains that site.
- `plaquettes_by_bond::Dict{Tuple{Int, Int}, Vector{Int}}` — per-bond
  list of plaquette indices that have that bond on their boundary.
  Same canonical `(min, max)` keying as `bond_index_of`.
"""
struct LatticeCache{T<:AbstractFloat}
    bonds::Vector{Bond{2,T}}
    plaquettes::Vector{Plaquette{2,T}}
    bond_index_of::Dict{Tuple{Int,Int},Int}
    plaquettes_by_vertex::Vector{Vector{Int}}
    plaquettes_by_bond::Dict{Tuple{Int,Int},Vector{Int}}
end

@inline _edge_key(i::Int, j::Int) = (min(i, j), max(i, j))

"""
    _build_lattice_cache(lat::Lattice{Topo, T}) → LatticeCache{T}

Populate the full cache in a single pass over the lattice. Called
once on first cache access through [`_get_cache`](@ref).
"""
function _build_lattice_cache(lat::Lattice{Topo,T}) where {Topo,T}
    # ---- Bonds via the existing generator --------------------------
    bond_vec = collect(
        b for i in 1:num_sites(lat) for b in neighbor_bonds(lat, i) if b.j > b.i
    )
    bond_index_of = Dict{Tuple{Int,Int},Int}()
    sizehint!(bond_index_of, length(bond_vec))
    for (k, b) in enumerate(bond_vec)
        bond_index_of[_edge_key(b.i, b.j)] = k
    end

    # ---- Plaquettes via the existing materialiser ------------------
    plaq_vec = _materialise_plaquettes(lat)

    N = num_sites(lat)
    plaquettes_by_vertex = [Int[] for _ in 1:N]
    plaquettes_by_bond = Dict{Tuple{Int,Int},Vector{Int}}()

    for (p_idx, p) in enumerate(plaq_vec)
        vs = p.vertices
        n = length(vs)
        for v in vs
            push!(plaquettes_by_vertex[v], p_idx)
        end
        for k in 1:n
            key = _edge_key(vs[k], vs[mod1(k + 1, n)])
            lst = get!(() -> Int[], plaquettes_by_bond, key)
            push!(lst, p_idx)
        end
    end

    return LatticeCache{T}(
        bond_vec, plaq_vec, bond_index_of, plaquettes_by_vertex, plaquettes_by_bond
    )
end

"""
    _get_cache(lat::Lattice{Topo, T}) → LatticeCache{T}

Return the lattice's per-sample cache, populating it on first access.
"""
function _get_cache(lat::Lattice{Topo,T}) where {Topo,T}
    ref = lat.cache
    if ref[] === nothing
        ref[] = _build_lattice_cache(lat)
    end
    return ref[]::LatticeCache{T}
end
