"""
    element_api.jl

Specialised element-center and incidence methods for
[`Lattice2D.Lattice`](@ref) that pull from the lazy
[`LatticeCache`](@ref) instead of materialising everything on every
call. Include order-wise this file has to come after `lattice.jl`
(for the `Lattice` type and `_materialise_plaquettes`) and
`cache.jl` (for `LatticeCache` and `_get_cache`).

The overrides here replace LatticeCore's correctness-only O(N)
defaults — the lattice still implements the same contract, callers
just stop paying O(num_bonds) per element-position lookup.

## Naming convention: element-API vs. Lattice2D helpers

Two coexisting query styles are exposed for bonds and plaquettes:

* **Generic, LatticeCore-owned element API** —
  `num_elements(lat, BondCenter())`, `element_position(lat, BondCenter(), i)`,
  `elements(lat, BondCenter())`, `incident(lat, BondCenter(), VertexCenter(), i)`,
  and so on. These are uniform across centring kinds (`VertexCenter`,
  `BondCenter`, `PlaquetteCenter`, …) and are what generic algorithms
  written against the `AbstractLattice` interface should call.

* **Lattice2D-local convenience helpers** — `num_bonds(lat)`,
  `num_plaquettes(lat)`, `bond_type(lat, i, j)`. These are thin,
  human-friendly aliases over the element-API methods above; they are
  exported by `Lattice2D` (not `LatticeCore`) and exist so that
  hand-written code at the Lattice2D layer can read naturally
  (`num_bonds(lat)` reads better than
  `num_elements(lat, BondCenter())` at a call site that only ever
  cares about bonds).

The two never disagree by construction (`num_bonds` literally calls
`num_elements(lat, BondCenter())`). A future refactor may move the
helpers under a sub-namespace to make the distinction more explicit;
that change is breaking and tracked separately.
"""

# ---- Bond and plaquette enumeration through the cache --------------

LatticeCore.bonds(lat::Lattice) = _get_cache(lat).bonds

LatticeCore.plaquettes(lat::Lattice) = _get_cache(lat).plaquettes

# ---- Element counts ------------------------------------------------

LatticeCore.num_elements(lat::Lattice, ::VertexCenter) = num_sites(lat)
LatticeCore.num_elements(lat::Lattice, ::BondCenter) = length(_get_cache(lat).bonds)
function LatticeCore.num_elements(lat::Lattice, ::PlaquetteCenter)
    return length(_get_cache(lat).plaquettes)
end

"""
    num_bonds(lat::Lattice) → Int

Convenience alias for `num_elements(lat, BondCenter())` / the length
of the cached bond Vector. O(1) after the first cache access.
"""
num_bonds(lat::Lattice) = num_elements(lat, BondCenter())

"""
    num_plaquettes(lat::Lattice) → Int

Convenience alias for `num_elements(lat, PlaquetteCenter())`.
Returns 0 for topologies whose unit cell declares no
[`PlaquetteRule`](@ref) (currently this is just `Dice` — wait, now
it is not; every shipped topology has a plaquette list as of Iter 6).
"""
num_plaquettes(lat::Lattice) = num_elements(lat, PlaquetteCenter())

"""
    bond_type(lat::Lattice, i::Int, j::Int) → Symbol

Return the bond-type tag of the edge connecting sites `i` and `j` on
`lat`. Uses the cached bond-index reverse lookup, so each call is
O(1) after the first cache access. Throws `ArgumentError` if there is
no declared bond between `i` and `j`.

The tag is inherited from the `Connection.type` of the underlying
unit cell, so it distinguishes anisotropic edges such as the
ShastrySutherland J′ dimer (`:type_2`) from its square-lattice NN
bonds (`:type_1`).
"""
function bond_type(lat::Lattice, i::Int, j::Int)
    cache = _get_cache(lat)
    key = _edge_key(i, j)
    idx = get(cache.bond_index_of, key, 0)
    idx == 0 &&
        throw(ArgumentError("no bond between sites $i and $j on $(typeof(lat).name.name)"))
    return cache.bonds[idx].type
end

# ---- Element iteration --------------------------------------------

LatticeCore.elements(lat::Lattice, ::VertexCenter) = 1:num_sites(lat)
LatticeCore.elements(lat::Lattice, ::BondCenter) = _get_cache(lat).bonds
LatticeCore.elements(lat::Lattice, ::PlaquetteCenter) = _get_cache(lat).plaquettes

# ---- Element positions: O(1) via cached Vector --------------------

# VertexCenter override: explicit so the call site does not have to
# round-trip through LatticeCore's generic AbstractLattice fallback
# (which already calls `position(lat, i)`). The behaviour is the same,
# this method just makes the element-API contract self-documenting on
# `Lattice` and matches the issue #42 request for an explicit
# specialisation.
LatticeCore.element_position(lat::Lattice, ::VertexCenter, i::Int) = position(lat, i)

function LatticeCore.element_position(lat::Lattice, ::BondCenter, i::Int)
    b = _get_cache(lat).bonds[i]
    return bond_center(lat, b)
end

function LatticeCore.element_position(lat::Lattice, ::PlaquetteCenter, i::Int)
    return _get_cache(lat).plaquettes[i].center
end

# ---- Element-position iterators: lazy, no full materialisation ----
#
# LatticeCore's generic `element_positions(lat, e)` is a lazy
# generator over `1:num_elements(lat, e)` that calls
# `element_position(lat, e, i)` on each step. For `Lattice` we can do
# better: `BondCenter` and `PlaquetteCenter` already have a cached
# Vector, so we hand back a generator that indexes the cache directly
# (no per-step `_nth` walk, no extra dispatch through the generic
# fallback). For `VertexCenter` we forward to `positions(lat)` so the
# call shares whatever iteration strategy the lattice picks for sites.

LatticeCore.element_positions(lat::Lattice, ::VertexCenter) = positions(lat)

function LatticeCore.element_positions(lat::Lattice, ::BondCenter)
    bs = _get_cache(lat).bonds
    return (bond_center(lat, b) for b in bs)
end

function LatticeCore.element_positions(lat::Lattice, ::PlaquetteCenter)
    ps = _get_cache(lat).plaquettes
    return (p.center for p in ps)
end

# ---- Same-centring adjacency: O(degree) / O(boundary) -------------

function LatticeCore.element_neighbors(lat::Lattice, ::BondCenter, i::Int)
    cache = _get_cache(lat)
    b = cache.bonds[i]
    seen = Set{Int}()
    # Bonds incident to either endpoint, excluding b itself.
    for endpoint in (b.i, b.j)
        for nb in neighbor_bonds(lat, endpoint)
            key = _edge_key(nb.i, nb.j)
            idx = get(cache.bond_index_of, key, 0)
            idx == 0 && continue
            idx == i && continue
            push!(seen, idx)
        end
    end
    return sort!(collect(seen))
end

function LatticeCore.element_neighbors(lat::Lattice, ::PlaquetteCenter, i::Int)
    cache = _get_cache(lat)
    p = cache.plaquettes[i]
    seen = Set{Int}()
    vs = p.vertices
    n = length(vs)
    for k in 1:n
        key = _edge_key(vs[k], vs[mod1(k + 1, n)])
        for q in get(cache.plaquettes_by_bond, key, Int[])
            q == i && continue
            push!(seen, q)
        end
    end
    return sort!(collect(seen))
end

# ---- Cross-centring incidence: O(local) ---------------------------

# Vertex → Bond: bonds touching site i. `neighbor_bonds` gives us
# the incident Bond objects directly; translate each to its canonical
# integer index via the cache's reverse-lookup dict.
function LatticeCore.incident(lat::Lattice, ::VertexCenter, ::BondCenter, i::Int)
    cache = _get_cache(lat)
    out = Int[]
    for b in neighbor_bonds(lat, i)
        key = _edge_key(b.i, b.j)
        idx = get(cache.bond_index_of, key, 0)
        idx == 0 || push!(out, idx)
    end
    return out
end

# Bond → Vertex: endpoints are in the cached Bond object.
function LatticeCore.incident(lat::Lattice, ::BondCenter, ::VertexCenter, i::Int)
    b = _get_cache(lat).bonds[i]
    return [b.i, b.j]
end

# Vertex → Plaquette: cached per-vertex list.
function LatticeCore.incident(lat::Lattice, ::VertexCenter, ::PlaquetteCenter, i::Int)
    return copy(_get_cache(lat).plaquettes_by_vertex[i])
end

# Plaquette → Vertex: cached boundary walk.
function LatticeCore.incident(lat::Lattice, ::PlaquetteCenter, ::VertexCenter, i::Int)
    return copy(_get_cache(lat).plaquettes[i].vertices)
end

# Bond → Plaquette: cached per-bond list.
function LatticeCore.incident(lat::Lattice, ::BondCenter, ::PlaquetteCenter, i::Int)
    cache = _get_cache(lat)
    b = cache.bonds[i]
    key = _edge_key(b.i, b.j)
    return copy(get(cache.plaquettes_by_bond, key, Int[]))
end

# Plaquette → Bond: walk the plaquette's boundary edges, map each to
# a bond index via the cache's edge → bond dict.
function LatticeCore.incident(lat::Lattice, ::PlaquetteCenter, ::BondCenter, i::Int)
    cache = _get_cache(lat)
    p = cache.plaquettes[i]
    vs = p.vertices
    n = length(vs)
    out = Int[]
    for k in 1:n
        key = _edge_key(vs[k], vs[mod1(k + 1, n)])
        idx = get(cache.bond_index_of, key, 0)
        idx == 0 || push!(out, idx)
    end
    return out
end
