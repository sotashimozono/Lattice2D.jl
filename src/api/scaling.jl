# Scale changes for Bravais-family lattices: the LatticeCore scaling interface.
#
# A `Lattice` is fully described by `(Topo, Lx, Ly, boundary, indexing, layout)`, so a scale step is
# just a change of `(Lx, Ly)` with everything else carried over. The one thing that cannot always be
# carried is the layout: a layout that stores data per *site* is sized to the old lattice.

# One step doubles the per-axis cell counts. Halving is available whenever the counts allow it, which
# `rescale` checks rather than rounding.
LatticeCore.scaling_rule(::Lattice) = LinearScaling(2)

# A layout may be transported to a lattice of a different size only if it carries no per-site data.
# `UniformLayout` stores one site type; `SublatticeLayout` stores one per sublattice — both are
# properties of the topology, not of the size. Anything else (e.g. an explicit per-site table) is
# sized to the lattice it was built for and must not be silently reused.
_size_independent_layout(::UniformLayout) = true
_size_independent_layout(::SublatticeLayout) = true
_size_independent_layout(::AbstractSiteLayout) = false

function _rebuild(lat::Lattice{Topo,T}, Lx::Int, Ly::Int) where {Topo,T}
    _size_independent_layout(lat.layout) || throw(
        ArgumentError(
            "cannot rescale a lattice whose layout carries per-site data " *
            "($(typeof(lat.layout))): it is sized to the current $(lat.Lx)×$(lat.Ly) lattice. " *
            "Rebuild the layout for the target size and construct the lattice directly.",
        ),
    )
    return Lattice{Topo,T}(
        Lx, Ly, lat.boundary, lat.indexing, lat.layout, Ref{Any}(nothing)
    )
end

"""
    rescale(lat::Lattice, k::Integer = 1) -> Lattice

The same lattice `k` scale steps larger: `(Lx, Ly)` are multiplied by `factor^k`, where `factor`
comes from [`scaling_rule`](@ref). Topology, boundary condition, indexing and layout are carried
over unchanged.

`k = 0` returns `lat`. Negative `k` steps *down*, and requires both `Lx` and `Ly` to be divisible by
`factor^|k|` — a lattice that cannot be halved exactly raises rather than rounding, since a silently
rounded size would corrupt any sequence built on it.

```jldoctest
julia> using Lattice2D, LatticeCore

julia> lat = build_lattice(Square, 2, 3);

julia> l2 = rescale(lat, 2);

julia> (l2.Lx, l2.Ly)
(8, 12)

julia> num_sites.(size_sequence(lat, 2))
3-element Vector{Int64}:
  6
 24
 96
```
"""
function LatticeCore.rescale(lat::Lattice, k::Integer=1)
    k == 0 && return lat
    f = LatticeCore.scaling_rule(lat).factor
    if k > 0
        m = f^k
        return _rebuild(lat, lat.Lx * m, lat.Ly * m)
    end
    m = f^(-k)
    (lat.Lx % m == 0 && lat.Ly % m == 0) || throw(
        ArgumentError(
            "cannot step down $(-k) level(s): $(lat.Lx)×$(lat.Ly) is not divisible by $m",
        ),
    )
    return _rebuild(lat, lat.Lx ÷ m, lat.Ly ÷ m)
end

"""
    cell_partition(lat::Lattice, k::Integer = 1) -> Vector{Vector{Int}}

Group the sites of `lat` by which unit cell of the lattice `k` steps *coarser* — that is, of
`rescale(lat, -k)` — they fall into. Entry `c` lists the site indices of `lat` inside cell `c` of the
coarser lattice, indexed with `lat`'s own [`AbstractIndexing`](@ref); together the groups partition
`1:num_sites(lat)`.

Each group holds `factor^k × factor^k` unit cells of `lat`, hence
`factor^(2k) * num_sublattices(lat)` sites. Requires `Lx` and `Ly` to be divisible by `factor^k`,
for the same reason `rescale` does.

```jldoctest
julia> using Lattice2D, LatticeCore

julia> groups = cell_partition(build_lattice(Square, 4, 4), 1);

julia> length(groups), length(first(groups))
(4, 4)

julia> sort(reduce(vcat, groups)) == 1:16
true
```
"""
function LatticeCore.cell_partition(lat::Lattice, k::Integer=1)
    k >= 1 ||
        throw(ArgumentError("cell_partition needs k ≥ 1 (a coarser lattice); got k = $k"))
    f = LatticeCore.scaling_rule(lat).factor
    m = f^k
    (lat.Lx % m == 0 && lat.Ly % m == 0) || throw(
        ArgumentError(
            "cannot coarsen by $k level(s): $(lat.Lx)×$(lat.Ly) is not divisible by $m"
        ),
    )
    nsub = num_sublattices(lat)
    cdims = (lat.Lx ÷ m, lat.Ly ÷ m)
    groups = [Int[] for _ in 1:prod(cdims)]
    for i in 1:num_sites(lat)
        lc = lattice_coord(lat.indexing, (lat.Lx, lat.Ly), nsub, i)
        gx = (lc.cell[1] - 1) ÷ m + 1
        gy = (lc.cell[2] - 1) ÷ m + 1
        g = site_index(lat.indexing, cdims, 1, LatticeCoord((gx, gy), 1))
        push!(groups[g], i)
    end
    return groups
end
