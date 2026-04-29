"""
    _edge_key(i::Int, j::Int) → Tuple{Int, Int}

Canonical, orientation-independent key for an undirected edge between
sites `i` and `j`, defined as `(min(i, j), max(i, j))`.

This helper is used by the per-lattice [`LatticeCache`](@ref) and the
element / incidence API in `core/element_api.jl` to look up bonds and
plaquette-by-bond entries without caring which endpoint walked first.

It is shared as a single internal utility (rather than re-defined ad
hoc in each file) so that any future change to the canonical-form
convention (e.g. promoting it to a struct) only needs to happen here.

This is **internal** to `Lattice2D` — it is not exported and should
not be relied on by downstream packages. If a similar helper is needed
in `LatticeCore` itself it should be promoted there in a separate PR.
"""
@inline _edge_key(i::Int, j::Int) = (min(i, j), max(i, j))
