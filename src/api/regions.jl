"""
    Lattice2D edge / bulk region API
    ================================

Helpers that split the sites (and bonds) of a [`Lattice`](@ref) into an
**edge** region and a **bulk** region. They exist to support
symmetry-protected-topological (SPT) analyses — e.g. separating the edge
modes of a Haldane / Kitaev chain from the gapped bulk — where an explicit
edge/bulk partition of the sample is needed (see issue #38).

The partition is defined purely from the declared connectivity, so it works
for any boundary condition and any topology without geometric thresholds:

* A **boundary site** is one that is *under-coordinated* — its coordination
  number is strictly below the bulk coordination number of its own
  sublattice. On a fully periodic sample every site is fully coordinated, so
  the boundary set (and hence the edge) is empty.
* Comparing per **sublattice** (rather than to a single global maximum) keeps
  the partition correct on lattices whose bulk coordination is not uniform
  across sublattices (e.g. Lieb, dice), where a low-coordination *bulk* site
  must not be mistaken for an edge site.
* `depth` peels successive layers inward: `depth = 1` is exactly the
  under-coordinated boundary ring; each increment adds every site within one
  more graph step of that ring.

All three helpers delegate to existing primitives ([`neighbors`](@ref),
[`bonds`](@ref), [`coordination_number`](@ref), [`sublattice`](@ref)).
"""

# ---- boundary detection ---------------------------------------------

"""
    _boundary_sites(lat::Lattice) -> Vector{Int}

Sorted site indices that are under-coordinated relative to the bulk
coordination number of their own sublattice. Internal helper backing
[`edge_sites`](@ref) at `depth = 1`.
"""
function _boundary_sites(lat::Lattice)
    N = num_sites(lat)
    N == 0 && return Int[]
    deg = coordination_number(lat)                      # per-site degree
    # Bulk coordination = maximum degree seen within each sublattice.
    bulk_deg = Dict{Int,Int}()
    @inbounds for i in 1:N
        s = sublattice(lat, i)
        d = deg[i]
        bulk_deg[s] = haskey(bulk_deg, s) ? max(bulk_deg[s], d) : d
    end
    out = Int[]
    @inbounds for i in 1:N
        if deg[i] < bulk_deg[sublattice(lat, i)]
            push!(out, i)
        end
    end
    return out
end

# ---- edge / bulk sites ----------------------------------------------

"""
    edge_sites(lat::Lattice; depth::Int = 1) -> Vector{Int}

Return the sorted site indices that make up the **edge** region of `lat`.

A site is on the edge if it lies within graph distance `depth - 1` of an
*under-coordinated* boundary site (a site whose coordination number is below
the bulk coordination of its sublattice). With `depth = 1` this is exactly
the boundary ring; each increment of `depth` peels one more layer of sites
inward.

On a fully periodic sample there are no under-coordinated sites, so the edge
is empty for every `depth`.

`depth` must be `≥ 1`.

# Example
```julia
chain = square(10, 1; boundary=OpenAxis())   # a 1D open chain
edge_sites(chain)                             # -> [1, 10]  (the two ends)
edge_sites(chain; depth=2)                    # -> [1, 2, 9, 10]
```

See also [`bulk_sites`](@ref), [`edge_bonds`](@ref).
"""
function edge_sites(lat::Lattice; depth::Int=1)
    depth ≥ 1 || throw(ArgumentError("depth must be ≥ 1, got $depth"))
    N = num_sites(lat)
    N == 0 && return Int[]
    boundary_ring = _boundary_sites(lat)
    isempty(boundary_ring) && return Int[]

    # BFS outward from the boundary ring, recording graph distance, and keep
    # every site reached within `depth - 1` steps.
    dist = fill(-1, N)
    queue = Int[]
    for i in boundary_ring
        dist[i] = 0
        push!(queue, i)
    end
    head = 1
    while head ≤ length(queue)
        u = queue[head]
        head += 1
        dist[u] == depth - 1 && continue        # do not expand past the last layer
        for v in neighbors(lat, u)
            if dist[v] == -1
                dist[v] = dist[u] + 1
                push!(queue, v)
            end
        end
    end
    return findall(d -> 0 ≤ d ≤ depth - 1, dist)
end

"""
    bulk_sites(lat::Lattice; depth::Int = 1) -> Vector{Int}

Return the sorted site indices in the **bulk** region of `lat`, i.e. the
complement of [`edge_sites`](@ref)`(lat; depth)`. On a fully periodic sample
this is every site.

`depth` must be `≥ 1`.

See also [`edge_sites`](@ref).
"""
function bulk_sites(lat::Lattice; depth::Int=1)
    edge = Set(edge_sites(lat; depth=depth))
    return [i for i in 1:num_sites(lat) if i ∉ edge]
end

# ---- edge bonds -----------------------------------------------------

"""
    edge_bonds(lat::Lattice; depth::Int = 1) -> Vector{Bond}

Return the bonds of `lat` that are incident to at least one edge site (see
[`edge_sites`](@ref)). These are the bonds touching the boundary region —
useful for restricting an environment / boundary Hamiltonian to the edge.

Bonds are returned in the order produced by [`bonds`](@ref). On a fully
periodic sample the edge is empty, so this returns an empty vector.

`depth` must be `≥ 1`.

See also [`edge_sites`](@ref), [`bulk_sites`](@ref).
"""
function edge_bonds(lat::Lattice; depth::Int=1)
    edge = Set(edge_sites(lat; depth=depth))
    isempty(edge) && return Bond[]
    return filter(b -> b.i ∈ edge || b.j ∈ edge, bonds(lat))
end
