"""
    dilution.jl

Site- and bond-dilution wrappers for any `AbstractLattice{D,T}`.

The wrapper type [`DilutedLattice`](@ref) takes an existing lattice
("base") plus

* `active_sites::BitVector` — `true` for sites kept on the diluted
  lattice (length = `num_sites(base)`)
* `killed_bonds::Set{Tuple{Int, Int}}` — canonical `(min, max)` pairs
  of *base* site indices whose bond is removed even when both
  endpoints are active

and exposes the standard `AbstractLattice` interface in terms of a
**remapped** index space `1..N_active`. Removed sites are not skipped
or masked — they disappear from the API, and the remaining sites are
renumbered contiguously. This keeps every consumer that walks
`1:num_sites(lat)` correct without it having to know about the
underlying disorder.

Two constructors realise the standard random-disorder ensembles:

* [`dilute_sites`](@ref) — independently delete each site with
  probability `p`
* [`dilute_bonds`](@ref) — keep all sites; independently delete each
  bond with probability `p`

Both accept an `rng` keyword for reproducibility (default
`Random.default_rng()`).

The wrapper is itself an `AbstractLattice{D,T}`, so it composes:
diluting a `DilutedLattice` again is well-defined and yields a new
`DilutedLattice` over the same `base` (the inner masks compose).
"""

# Note: this file deliberately does not `using Random`. `Random` is in
# Julia's stdlib but is **not** declared in Lattice2D's `[deps]` —
# loading it from `src/` would fail at precompile. The disorder
# generators below therefore accept an opaque `rng` keyword and
# forward it to `rand(rng, ...)`, which works for any AbstractRNG
# the caller chooses to plumb in (e.g. `Random.MersenneTwister(seed)`
# from the test or downstream package). The default `rng = nothing`
# falls back to the argument-less `rand()` (Julia's own thread-local
# RNG), which keeps the API ergonomic without forcing `Random` to be
# a hard dependency. See issue #36 / discussion in the dilution PR.

# ---- Type ----------------------------------------------------------

"""
    DilutedLattice{D, T, L<:AbstractLattice{D, T}} <: AbstractLattice{D, T}

Wrapper lattice with disordered sites / bonds. See module docstring of
`disorder/dilution.jl`.

Fields
------

- `base::L` — the underlying lattice
- `active_sites::BitVector` — length `num_sites(base)`; `true` where
  the site is kept
- `killed_bonds::Set{Tuple{Int, Int}}` — canonical `(min_i, max_j)`
  base-index pairs of bonds that are explicitly removed (in addition
  to bonds whose endpoints are inactive)
- `new_to_old::Vector{Int}` — `new_to_old[k]` is the base-lattice
  index of the `k`-th active site (`length == count(active_sites)`)
- `old_to_new::Vector{Int}` — `old_to_new[i]` is the new index of
  base site `i`, or `0` if `i` is inactive

The `new_to_old` / `old_to_new` arrays are derived from
`active_sites` and cached at construction; the constructors below are
the only sanctioned entry points so the invariants stay tight.
"""
struct DilutedLattice{D,T,L<:AbstractLattice{D,T}} <: AbstractLattice{D,T}
    base::L
    active_sites::BitVector
    killed_bonds::Set{Tuple{Int,Int}}
    new_to_old::Vector{Int}
    old_to_new::Vector{Int}
end

"""
    DilutedLattice(base, active_sites, killed_bonds = Set{Tuple{Int, Int}}())

Construct a [`DilutedLattice`](@ref) from a base lattice plus a
site-activity mask and an optional explicit bond-removal set. The
remap arrays are derived from `active_sites`. `killed_bonds` keys are
canonicalised to `(min, max)` form.

Throws `ArgumentError` if `active_sites` does not match
`num_sites(base)`, or if any entry of `killed_bonds` references a
site outside `1:num_sites(base)`.
"""
function DilutedLattice(
    base::AbstractLattice{D,T},
    active_sites::BitVector,
    killed_bonds::Set{Tuple{Int,Int}}=Set{Tuple{Int,Int}}(),
) where {D,T}
    N = num_sites(base)
    length(active_sites) == N || throw(
        ArgumentError(
            "active_sites length $(length(active_sites)) does not match num_sites(base) = $N",
        ),
    )

    # Canonicalise killed_bonds and validate endpoint range.
    canon = Set{Tuple{Int,Int}}()
    for (a, b) in killed_bonds
        (1 ≤ a ≤ N && 1 ≤ b ≤ N) ||
            throw(ArgumentError("killed_bonds entry ($a, $b) out of range 1..$N"))
        a == b && throw(ArgumentError("killed_bonds entry ($a, $b) is a self-loop"))
        push!(canon, (min(a, b), max(a, b)))
    end

    new_to_old = Int[]
    old_to_new = zeros(Int, N)
    for i in 1:N
        if active_sites[i]
            push!(new_to_old, i)
            old_to_new[i] = length(new_to_old)
        end
    end

    L = typeof(base)
    return DilutedLattice{D,T,L}(base, active_sites, canon, new_to_old, old_to_new)
end

# ---- Random constructors -------------------------------------------

"""
    dilute_sites(base::AbstractLattice, p::Real; rng=nothing) → DilutedLattice

Delete each site of `base` independently with probability `p`. The
returned [`DilutedLattice`](@ref) renumbers the surviving sites
contiguously as `1..N_active`.

`p` must lie in `[0, 1]`. `p = 0` reproduces the base lattice
(modulo the wrapper); `p = 1` removes every site (the wrapper is
still well-formed but most queries become empty).

`rng` is forwarded as the first positional argument to `rand`, so any
`AbstractRNG` (e.g. `Random.MersenneTwister(seed)`) may be passed for
reproducibility. The default `nothing` falls back to the argument-less
`rand()` (Julia's thread-local RNG).
"""
function dilute_sites(base::AbstractLattice, p::Real; rng=nothing)
    0 ≤ p ≤ 1 || throw(ArgumentError("dilution probability p must lie in [0, 1], got $p"))
    N = num_sites(base)
    active = BitVector(undef, N)
    if rng === nothing
        @inbounds for i in 1:N
            active[i] = rand() ≥ p
        end
    else
        @inbounds for i in 1:N
            active[i] = rand(rng) ≥ p
        end
    end
    return DilutedLattice(base, active)
end

"""
    dilute_bonds(base::AbstractLattice, p::Real; rng=nothing) → DilutedLattice

Keep all sites; delete each bond of `base` independently with
probability `p`. The returned [`DilutedLattice`](@ref) shares the
site indexing of `base` (`num_sites` and `position` are unchanged)
but its bond / neighbour iterators omit the killed edges.

`p` must lie in `[0, 1]`. `p = 0` keeps every bond; `p = 1` deletes
every bond.

`rng` is forwarded as the first positional argument to `rand`, so any
`AbstractRNG` may be passed for reproducibility. The default
`nothing` falls back to the argument-less `rand()`.
"""
function dilute_bonds(base::AbstractLattice, p::Real; rng=nothing)
    0 ≤ p ≤ 1 || throw(ArgumentError("dilution probability p must lie in [0, 1], got $p"))
    N = num_sites(base)
    active = trues(N)
    killed = Set{Tuple{Int,Int}}()
    if rng === nothing
        for b in bonds(base)
            if rand() < p
                push!(killed, (min(b.i, b.j), max(b.i, b.j)))
            end
        end
    else
        for b in bonds(base)
            if rand(rng) < p
                push!(killed, (min(b.i, b.j), max(b.i, b.j)))
            end
        end
    end
    return DilutedLattice(base, active, killed)
end

# ---- Internal helpers ----------------------------------------------

@inline _bond_alive(lat::DilutedLattice, oi::Int, oj::Int) = begin
    lat.active_sites[oi] &&
        lat.active_sites[oj] &&
        !((min(oi, oj), max(oi, oj)) in lat.killed_bonds)
end

# ---- AbstractLattice interface -------------------------------------

LatticeCore.num_sites(lat::DilutedLattice) = length(lat.new_to_old)

function LatticeCore.position(lat::DilutedLattice, i::Int)
    1 ≤ i ≤ num_sites(lat) || throw(BoundsError(lat, i))
    return position(lat.base, lat.new_to_old[i])
end

LatticeCore.boundary(lat::DilutedLattice) = boundary(lat.base)

LatticeCore.size_trait(lat::DilutedLattice) = size_trait(lat.base)

# Sublattice info follows the base lattice's geometric layout (the
# wrapper does not invent new sublattices, it just removes sites).
LatticeCore.num_sublattices(lat::DilutedLattice) = num_sublattices(lat.base)
function LatticeCore.sublattice(lat::DilutedLattice, i::Int)
    return sublattice(lat.base, lat.new_to_old[i])
end

# Topology / periodicity / reciprocal traits: dilution does not
# preserve translation symmetry, but it does preserve the *embedding*
# topology / boundary conditions of the underlying mesh. Downstream
# code that needs Bloch-style reciprocal queries should use the base
# lattice; we conservatively report `:disordered` topology and
# `Aperiodic` periodicity / `NoReciprocal` reciprocal support so no
# generic algorithm assumes Bravais structure on a diluted sample.
LatticeCore.topology(::DilutedLattice) = TopologyTrait{:disordered}()
LatticeCore.periodicity(::DilutedLattice) = Aperiodic()
LatticeCore.reciprocal_support(::DilutedLattice) = NoReciprocal()

# Bipartite check: re-run the BFS on the diluted graph (the base
# lattice's `is_bipartite` answers the wrong question once edges are
# missing or vertices removed).
function LatticeCore.is_bipartite(lat::DilutedLattice)
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

# Site layout: defer to the base layout but reindex through
# `new_to_old`. We deliberately do **not** override `site_layout` to
# avoid mutating layout state; instead we override `site_type`
# directly so the wrapper's site index `i` resolves to the correct
# base site type.
LatticeCore.site_type(lat::DilutedLattice, i::Int) = site_type(lat.base, lat.new_to_old[i])

# `site_layout` is required by some `AbstractLattice` consumers; we
# return the base layout. Because layouts are indexed by site, the
# returned object should not be queried with diluted indices — use
# `site_type(lat, i)` for that. The contract here is: `site_layout`
# tells you *what kind* of layout the lattice has, not how to index
# it.
LatticeCore.site_layout(lat::DilutedLattice) = site_layout(lat.base)

# ---- Neighbours ----------------------------------------------------

function LatticeCore.neighbors(lat::DilutedLattice, i::Int)
    1 ≤ i ≤ num_sites(lat) || throw(BoundsError(lat, i))
    oi = lat.new_to_old[i]
    out = Int[]
    for oj in neighbors(lat.base, oi)
        _bond_alive(lat, oi, oj) || continue
        nj = lat.old_to_new[oj]
        nj == 0 && continue
        push!(out, nj)
    end
    return out
end

function LatticeCore.neighbor_bonds(lat::DilutedLattice{D,T}, i::Int) where {D,T}
    1 ≤ i ≤ num_sites(lat) || throw(BoundsError(lat, i))
    oi = lat.new_to_old[i]
    out = Bond{D,T}[]
    for b in neighbor_bonds(lat.base, oi)
        # `neighbor_bonds(base, oi)` always yields bonds with `b.i == oi`.
        oj = b.j
        _bond_alive(lat, oi, oj) || continue
        nj = lat.old_to_new[oj]
        nj == 0 && continue
        push!(out, Bond{D,T}(i, nj, b.vector, b.type))
    end
    return out
end

# ---- Bonds ---------------------------------------------------------

function LatticeCore.bonds(lat::DilutedLattice{D,T}) where {D,T}
    out = Bond{D,T}[]
    for b in bonds(lat.base)
        _bond_alive(lat, b.i, b.j) || continue
        ni = lat.old_to_new[b.i]
        nj = lat.old_to_new[b.j]
        (ni == 0 || nj == 0) && continue
        push!(out, Bond{D,T}(ni, nj, b.vector, b.type))
    end
    return out
end

# ---- Element API ---------------------------------------------------
#
# Override `num_elements`, `elements`, `element_position` for vertex
# and bond centres so the LatticeCore generic fallbacks don't try to
# walk through the *base* `bonds` iterator and mis-count.

LatticeCore.num_elements(lat::DilutedLattice, ::VertexCenter) = num_sites(lat)
LatticeCore.num_elements(lat::DilutedLattice, ::BondCenter) = length(bonds(lat))

LatticeCore.elements(lat::DilutedLattice, ::VertexCenter) = 1:num_sites(lat)
LatticeCore.elements(lat::DilutedLattice, ::BondCenter) = bonds(lat)

LatticeCore.element_position(lat::DilutedLattice, ::VertexCenter, i::Int) = position(lat, i)

function LatticeCore.element_position(lat::DilutedLattice, ::BondCenter, i::Int)
    bs = bonds(lat)
    1 ≤ i ≤ length(bs) || throw(BoundsError(bs, i))
    return bond_center(lat, bs[i])
end

# ---- Iteration / display -------------------------------------------

Base.length(lat::DilutedLattice) = num_sites(lat)

function Base.iterate(lat::DilutedLattice, state=1)
    return state > num_sites(lat) ? nothing : (state, state + 1)
end

Base.eltype(::Type{<:DilutedLattice}) = Int

function Base.show(io::IO, lat::DilutedLattice)
    Nbase = num_sites(lat.base)
    Nact = num_sites(lat)
    nkill = length(lat.killed_bonds)
    return print(
        io,
        "DilutedLattice over $(typeof(lat.base).name.name): " *
        "$(Nact) / $(Nbase) sites active, $(nkill) explicit bond removals",
    )
end
