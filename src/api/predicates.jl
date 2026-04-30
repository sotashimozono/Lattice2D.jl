"""
    Lattice2D structural predicates
    ===============================

This file groups one-shot helpers that report structural / graph-theoretic
properties of a [`Lattice`](@ref) without touching its construction. They
are pure functions of the lattice value and its declared connectivity, and
they delegate to existing primitives ([`neighbors`](@ref), [`bonds`](@ref),
[`num_sites`](@ref)) so any future change to those primitives is reflected
here automatically.

The helpers are intentionally cheap to write call sites for — common
sanity checks ("is this bipartite?", "what is the mean coordination of
this sample?", "give me the first three neighbour shells") collapse to a
single call.
"""

# ---- is_bipartite ---------------------------------------------------
#
# Implementation moved here from `core/lattice.jl` to keep all
# structural-predicate methods in one file. The signature and behaviour
# are unchanged: BFS 2-colouring on the declared adjacency, returning
# `false` on the first odd cycle.

"""
    is_bipartite(lat::Lattice) -> Bool

Return `true` iff the declared adjacency graph of `lat` is bipartite,
i.e. its sites admit a 2-colouring such that every bond connects
different colours.

Implemented as a BFS 2-colouring over [`neighbors`](@ref), starting a
fresh BFS from each unvisited site so disconnected components are
handled. `O(num_sites + num_bonds)`.
"""
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

# ---- coordination number --------------------------------------------

"""
    coordination_number(lat::Lattice, i::Int) -> Int

Degree of site `i` in the declared adjacency graph of `lat` (number of
distinct neighbours, no double-counting). For PBC samples on uniform
lattices this matches the textbook coordination number (Square = 4,
Triangular = 6, Honeycomb = 3, Kagome = 4, ...). For OBC samples it
returns the local degree, which is smaller at the boundary.
"""
coordination_number(lat::Lattice, i::Int) = length(neighbors(lat, i))

"""
    coordination_number(lat::Lattice) -> Vector{Int}

Per-site degrees, in site-index order. Length `num_sites(lat)`.
"""
function coordination_number(lat::Lattice)
    N = num_sites(lat)
    out = Vector{Int}(undef, N)
    @inbounds for i in 1:N
        out[i] = coordination_number(lat, i)
    end
    return out
end

"""
    mean_coordination(lat::Lattice) -> Float64

Average coordination number across all sites of `lat`. For periodic
uniform lattices this equals the bulk coordination number; for samples
with open boundaries it is strictly smaller.
"""
function mean_coordination(lat::Lattice)
    N = num_sites(lat)
    N == 0 && return 0.0
    s = 0
    @inbounds for i in 1:N
        s += coordination_number(lat, i)
    end
    return s / N
end

# ---- bond distances -------------------------------------------------

"""
    bond_distances(lat::Lattice) -> Vector{Float64}

Euclidean length of every bond returned by [`bonds`](@ref), in the
same order. Lengths are not deduplicated, so the multiset of distances
encodes both the bond-type spectrum and the bond multiplicity.
"""
function bond_distances(lat::Lattice)
    bs = bonds(lat)
    out = Vector{Float64}(undef, length(bs))
    @inbounds for k in eachindex(bs)
        out[k] = Float64(norm(bs[k].vector))
    end
    return out
end

# ---- shells ---------------------------------------------------------

"""
    shells(lat::Lattice, i::Int; n_shells::Int = 3) -> Vector{Vector{Int}}

Return the first `n_shells` geometric neighbour shells of site `i` as
a list of site-index vectors. Each entry corresponds to one Euclidean
distance class, in increasing order; the `k`-th entry is exactly
`neighbors(lat, i; shell = k)`.

Trailing empty shells (e.g. on a small OBC sample where the requested
shell count exceeds the geometric reach) are returned as empty
vectors, so the result always has length `n_shells`.

`n_shells` must be `≥ 1`.
"""
function shells(lat::Lattice, i::Int; n_shells::Int=3)
    n_shells ≥ 1 || throw(ArgumentError("n_shells must be ≥ 1, got $n_shells"))
    out = Vector{Vector{Int}}(undef, n_shells)
    @inbounds for k in 1:n_shells
        out[k] = neighbors(lat, i; shell=k)
    end
    return out
end
