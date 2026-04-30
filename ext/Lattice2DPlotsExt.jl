module Lattice2DPlotsExt

using Lattice2D
using LatticeCore
using LinearAlgebra
using Plots
using StaticArrays

"""
    Lattice2DPlotsExt

Lattice2D-side companion to the upstream `LatticeCorePlotsExt`. The
upstream extension already provides `plot_lattice(lat)` and a
default `plot_bonds!(p, lat)` that draws every bond as one line
segment in a single colour. This file adds **bond-aware**
visualisation helpers that are specific to `Lattice2D.Lattice`:

- `plot_bonds(lat; bond_types, color_by, kwargs...)` — standalone,
  returns a fresh figure with bonds (and optionally sites) drawn.
- `plot_bonds!(p, lat::Lattice2D.Lattice; bond_types, color_by, ...)` —
  overlay form; specialises the LatticeCore generic with extra
  Lattice2D-only keywords without breaking call sites that only
  pass `(p, lat; color, lw)`.

Both routines support two grouping modes:

- `color_by = :type` (default) — colour each bond by its
  `Bond.type` symbol (e.g. `:type_1`, `:type_2` for
  `ShastrySutherland`'s J / J' dimers).
- `color_by = :direction` — colour each bond by the rounded
  orientation of its displacement vector. Useful on lattices where
  every bond shares one `:nearest` tag but the geometry has multiple
  inequivalent directions (e.g. `Honeycomb`, `Kagome`).

The new methods do **not** override the existing two-argument
`LatticeCore.plot_bonds!(p, lat; color, lw)` for non-Lattice2D
lattices — they only specialise on `lat::Lattice2D.Lattice`.

# Example

```julia
using Lattice2D, Plots

lat = build_lattice(ShastrySutherland, 4, 4)

# Standalone bond drawing, J vs J' coloured separately.
plot_bonds(lat; color_by=:type, with_sites=true)

# Overlay on an existing lattice plot.
p = plot_lattice(lat; show_bonds=false)
plot_bonds!(p, lat; color_by=:type, lw=2.0)
```
"""
Lattice2DPlotsExt

# ---- Internal helpers ------------------------------------------------

# Quantise a 2D bond vector onto an undirected orientation key in
# `[0, π)` so antiparallel bonds share a colour. Steps of π/12 (15°)
# are fine enough to separate the inequivalent directions on every
# shipped topology while still tolerating floating-point drift in
# `Bond.vector`.
const _DIRECTION_BIN_STEP = π / 12

@inline function _direction_key(v::SVector{2,T}) where {T}
    n = sqrt(v[1]^2 + v[2]^2)
    n == 0 && return 0
    θ = atan(v[2], v[1])  # [-π, π]
    # Fold antiparallel bonds onto the same key.
    θ < 0 && (θ += π)
    θ >= π && (θ -= π)
    return round(Int, θ / _DIRECTION_BIN_STEP)
end

# Per-bond grouping key, returned as a Symbol so we can hand it to
# `Plots.plot!(...; group=keys)` and get one legend entry per group.
function _bond_group_keys(lat::Lattice2D.Lattice, bs::AbstractVector, color_by::Symbol)
    if color_by === :type
        return Symbol[b.type for b in bs]
    elseif color_by === :direction
        return Symbol[Symbol("dir_", _direction_key(b.vector)) for b in bs]
    else
        throw(ArgumentError("color_by must be :type or :direction; got $(color_by)"))
    end
end

# Filter bonds by the user-supplied `bond_types` selector.
function _filter_bonds(lat::Lattice2D.Lattice, bond_types)
    bs = collect(bonds(lat))
    if bond_types === :all
        return bs
    elseif bond_types isa Symbol
        return [b for b in bs if b.type === bond_types]
    elseif bond_types isa AbstractVector || bond_types isa Tuple
        wanted = Set(Symbol.(collect(bond_types)))
        return [b for b in bs if b.type in wanted]
    else
        throw(
            ArgumentError(
                "bond_types must be :all, a Symbol, or a collection " *
                "of Symbols; got $(bond_types)",
            ),
        )
    end
end

# Build (seg_x, seg_y, group_keys_per_segment) by traversing the
# selected bonds. The per-segment group vector has the same
# `NaN`-separated layout as `seg_x` / `seg_y`: three entries per
# bond (two endpoints + NaN sentinel) sharing the same group key,
# which is what `Plots.plot!(...; group=...)` consumes.
function _grouped_segments(lat::Lattice2D.Lattice, bs::AbstractVector, color_by::Symbol)
    seg_x = Float64[]
    seg_y = Float64[]
    grp_x = Symbol[]
    grp_y = Symbol[]
    keys = _bond_group_keys(lat, bs, color_by)
    for (b, k) in zip(bs, keys)
        src = position(lat, b.i)
        dst = src + b.vector
        push!(seg_x, Float64(src[1]), Float64(dst[1]), NaN)
        push!(seg_y, Float64(src[2]), Float64(dst[2]), NaN)
        push!(grp_x, k, k, k)
        # `grp_y` mirrors `grp_x`; we keep both so the call site can
        # pass the matching key for either coordinate axis if needed.
        push!(grp_y, k, k, k)
    end
    return seg_x, seg_y, grp_x
end

# ---- plot_bonds! (overlay, Lattice2D-aware) -------------------------

"""
    LatticeCore.plot_bonds!(p, lat::Lattice2D.Lattice;
                            bond_types=:all, color_by=:type,
                            lw=1.5, kwargs...) → p

Overlay coloured bonds onto an existing `Plots.Plot` `p`.

# Keyword arguments
- `bond_types` — `:all` (default), a single `Symbol`, or a vector /
  tuple of `Symbol`s selecting which `Bond.type` tags to draw.
- `color_by` — `:type` (default) or `:direction`. Controls how bonds
  are grouped for colouring; see the module docstring.
- `lw` — line width forwarded to `Plots.plot!`.
- additional `kwargs...` are forwarded verbatim to `Plots.plot!`.

This method specialises the upstream `LatticeCore.plot_bonds!` for
`Lattice2D.Lattice` and is fully backward-compatible: callers that
omit `bond_types` / `color_by` get the same set of bonds drawn as
before, just with one legend entry per bond-type group.
"""
function LatticeCore.plot_bonds!(
    p, lat::Lattice2D.Lattice; bond_types=:all, color_by::Symbol=:type, lw=1.5, kwargs...
)
    bs = _filter_bonds(lat, bond_types)
    isempty(bs) && return p
    seg_x, seg_y, grp = _grouped_segments(lat, bs, color_by)
    return Plots.plot!(p, seg_x, seg_y; group=grp, lw=lw, kwargs...)
end

# ---- plot_bonds (standalone) ----------------------------------------

"""
    plot_bonds(lat::Lattice2D.Lattice; bond_types=:all, color_by=:type,
               with_sites=false, site_size=4, lw=1.5,
               aspect_ratio=:equal, kwargs...) → Plots.Plot

Build a fresh figure that draws the bonds of `lat` grouped (and thus
coloured) by `color_by`. Pass `with_sites=true` to overlay the site
markers from the upstream `plot_sites!` recipe on top of the bonds.

# Example: combining with the upstream lattice recipe

```julia
using Lattice2D, Plots

lat = build_lattice(ShastrySutherland, 4, 4)

# Standalone — bonds only, coloured by `Bond.type`.
plot_bonds(lat)

# Overlay on `plot_lattice` (which draws sites too).
p = plot_lattice(lat; show_bonds=false)
plot_bonds!(p, lat; color_by=:type)
```

See also: [`LatticeCore.plot_bonds!`](@ref).
"""
function Lattice2D.plot_bonds(
    lat::Lattice2D.Lattice;
    bond_types=:all,
    color_by::Symbol=:type,
    with_sites::Bool=false,
    site_size=4,
    lw=1.5,
    aspect_ratio=:equal,
    title=nothing,
    kwargs...,
)
    p = Plots.plot(;
        aspect_ratio=aspect_ratio,
        xlabel="x",
        ylabel="y",
        legend=:topright,
        title=title,
        kwargs...,
    )
    LatticeCore.plot_bonds!(p, lat; bond_types=bond_types, color_by=color_by, lw=lw)
    if with_sites
        LatticeCore.plot_sites!(p, lat; marker_size=site_size)
    end
    return p
end

# ---- plot_state -----------------------------------------------------

function _site_xy(lat::Lattice2D.Lattice)
    N = num_sites(lat)
    xs = Vector{Float64}(undef, N)
    ys = Vector{Float64}(undef, N)
    for i in 1:N
        p = position(lat, i)
        xs[i] = Float64(p[1])
        ys[i] = Float64(p[2])
    end
    return xs, ys
end

function _bond_segments(lat::Lattice2D.Lattice)
    seg_x = Float64[]
    seg_y = Float64[]
    for b in bonds(lat)
        src = position(lat, b.i)
        dst_local = src + b.vector
        push!(seg_x, Float64(src[1]), Float64(dst_local[1]), NaN)
        push!(seg_y, Float64(src[2]), Float64(dst_local[2]), NaN)
    end
    return seg_x, seg_y
end

_is_discrete_state(::AbstractVector{Bool}) = true
function _is_discrete_state(state::AbstractVector{<:Integer})
    return length(unique(state)) <= 12
end
_is_discrete_state(::AbstractVector) = false

const _DISCRETE_PALETTE = (
    :steelblue,
    :tomato,
    :seagreen,
    :goldenrod,
    :mediumpurple,
    :sandybrown,
    :hotpink,
    :gray,
    :olive,
    :teal,
    :crimson,
    :darkorange,
)

function Lattice2D.plot_state(
    lat::Lattice2D.Lattice,
    state::AbstractVector;
    colormap=:viridis,
    marker_size::Real=12,
    show_bonds::Bool=true,
    bond_color=:lightgray,
    bond_width::Real=1.0,
    title=nothing,
    aspect_ratio=:equal,
    kwargs...,
)
    N = num_sites(lat)
    length(state) == N || throw(
        DimensionMismatch(
            "plot_state expects length(state) == num_sites(lat); got $(length(state)) vs $N",
        ),
    )

    xs, ys = _site_xy(lat)

    p = Plots.plot(;
        aspect_ratio=aspect_ratio,
        xlabel="x",
        ylabel="y",
        title=title,
        legend=:topright,
        kwargs...,
    )

    if show_bonds
        seg_x, seg_y = _bond_segments(lat)
        Plots.plot!(p, seg_x, seg_y; color=bond_color, lw=bond_width, label="")
    end

    if _is_discrete_state(state)
        labels = sort(unique(state))
        for (k, lbl) in enumerate(labels)
            mask = state .== lbl
            colour = _DISCRETE_PALETTE[mod1(k, length(_DISCRETE_PALETTE))]
            Plots.scatter!(
                p,
                xs[mask],
                ys[mask];
                ms=marker_size,
                mc=colour,
                markerstrokewidth=0,
                label=string(lbl),
            )
        end
    else
        zs = Float64.(state)
        Plots.scatter!(
            p,
            xs,
            ys;
            ms=marker_size,
            marker_z=zs,
            color=colormap,
            markerstrokewidth=0,
            label="",
            colorbar=true,
        )
    end

    return p
end

# ---- Brillouin zone -------------------------------------------------

"""
    Lattice2D.brillouin_zone(lat::Lattice2D.Lattice; shell::Int=2)
        -> Vector{SVector{2,Float64}}

Compute the Brillouin zone of `lat` as the Wigner-Seitz cell of its
reciprocal lattice via half-plane intersection over the
`(m, n) ∈ [-shell, shell]² ∖ {(0,0)}` shell of reciprocal vectors.
Returns vertices ordered counter-clockwise around the origin.
"""
function Lattice2D.brillouin_zone(lat::Lattice2D.Lattice; shell::Int=2)
    LatticeCore.reciprocal_support(lat) isa LatticeCore.HasReciprocal || throw(
        ArgumentError(
            "brillouin_zone requires a fully periodic Lattice (got $(lat.boundary))"
        ),
    )
    ml = LatticeCore.reciprocal_lattice(lat)
    B = LatticeCore.reciprocal_basis(ml)
    return _brillouin_zone_from_basis(SMatrix{2,2,Float64}(B), shell)
end

"""
    Lattice2D.high_symmetry_points(lat::Lattice2D.Lattice)
        -> Dict{Symbol,SVector{2,Float64}}

Topology-keyed dictionary of high-symmetry points (Cartesian k).
Populated for `Square`, `Triangular`, `Honeycomb`; falls back to
`Dict(:Gamma => 0)` for other topologies.
"""
function Lattice2D.high_symmetry_points(lat::Lattice2D.Lattice)
    LatticeCore.reciprocal_support(lat) isa LatticeCore.HasReciprocal || throw(
        ArgumentError(
            "high_symmetry_points requires a fully periodic Lattice (got $(lat.boundary))",
        ),
    )
    ml = LatticeCore.reciprocal_lattice(lat)
    B = SMatrix{2,2,Float64}(LatticeCore.reciprocal_basis(ml))
    return _high_symmetry_points(LatticeCore.topology(lat), B)
end

const _BZ_TOL = 1e-10

function _brillouin_zone_from_basis(B::SMatrix{2,2,Float64}, shell::Int)
    Gs = SVector{2,Float64}[]
    for n in (-shell):shell, m in (-shell):shell
        (m == 0 && n == 0) && continue
        push!(Gs, B * SVector{2,Float64}(m, n))
    end

    verts = SVector{2,Float64}[]
    nG = length(Gs)
    for i in 1:nG
        Gi = Gs[i]
        rhs_i = dot(Gi, Gi) / 2
        for j in (i + 1):nG
            Gj = Gs[j]
            rhs_j = dot(Gj, Gj) / 2
            M = SMatrix{2,2,Float64}(Gi[1], Gj[1], Gi[2], Gj[2])
            d = det(M)
            abs(d) < _BZ_TOL && continue
            k = M \ SVector{2,Float64}(rhs_i, rhs_j)
            inside = true
            for G in Gs
                lhs = 2 * dot(k, G)
                rhs = dot(G, G)
                if lhs > rhs + _BZ_TOL * (1 + abs(rhs))
                    inside = false
                    break
                end
            end
            inside || continue
            already = false
            for v in verts
                if norm(v - k) < _BZ_TOL
                    already = true
                    break
                end
            end
            already || push!(verts, k)
        end
    end

    isempty(verts) && return verts

    cx = sum(v[1] for v in verts) / length(verts)
    cy = sum(v[2] for v in verts) / length(verts)
    sort!(verts; by=v -> atan(v[2] - cy, v[1] - cx))
    return verts
end

function _high_symmetry_points(
    ::LatticeCore.TopologyTrait{:square}, B::SMatrix{2,2,Float64}
)
    b1 = B[:, 1]
    b2 = B[:, 2]
    return Dict{Symbol,SVector{2,Float64}}(
        :Gamma => SVector{2,Float64}(0.0, 0.0),
        :X => SVector{2,Float64}(b1 / 2),
        :M => SVector{2,Float64}((b1 + b2) / 2),
    )
end

function _high_symmetry_points(
    ::LatticeCore.TopologyTrait{:triangular}, B::SMatrix{2,2,Float64}
)
    b1 = B[:, 1]
    b2 = B[:, 2]
    return Dict{Symbol,SVector{2,Float64}}(
        :Gamma => SVector{2,Float64}(0.0, 0.0),
        :M => SVector{2,Float64}(b1 / 2),
        :K => SVector{2,Float64}((b1 + 2b2) / 3),
        :K_prime => SVector{2,Float64}((2b1 + b2) / 3),
    )
end

function _high_symmetry_points(
    ::LatticeCore.TopologyTrait{:honeycomb}, B::SMatrix{2,2,Float64}
)
    return _high_symmetry_points(LatticeCore.TopologyTrait{:triangular}(), B)
end

function _high_symmetry_points(::LatticeCore.TopologyTrait, ::SMatrix{2,2,Float64})
    return Dict{Symbol,SVector{2,Float64}}(:Gamma => SVector{2,Float64}(0.0, 0.0))
end

function Lattice2D.plot_brillouin_zone(
    lat::Lattice2D.Lattice;
    show_mesh::Bool=false,
    ml::Union{Nothing,LatticeCore.AbstractMomentumLattice}=nothing,
    show_high_symmetry::Bool=false,
    shell::Int=2,
    bz_color=:black,
    bz_lw=1.5,
    mesh_color=:steelblue,
    mesh_size=3,
    label_color=:crimson,
    title=nothing,
    kwargs...,
)
    verts = Lattice2D.brillouin_zone(lat; shell=shell)
    isempty(verts) && throw(
        ErrorException(
            "Brillouin zone polygon is empty -- bump `shell` or check the basis"
        ),
    )

    xs = Float64[v[1] for v in verts]
    ys = Float64[v[2] for v in verts]
    push!(xs, xs[1])
    push!(ys, ys[1])

    p = Plots.plot(;
        aspect_ratio=:equal,
        xlabel="kx",
        ylabel="ky",
        legend=:topright,
        title=title === nothing ? "Brillouin zone" : title,
        kwargs...,
    )
    Plots.plot!(p, xs, ys; color=bz_color, lw=bz_lw, label="BZ")

    if show_mesh
        mesh = ml === nothing ? LatticeCore.reciprocal_lattice(lat) : ml
        M = LatticeCore.num_k_points(mesh)
        kxs = [Float64(LatticeCore.k_point(mesh, i)[1]) for i in 1:M]
        kys = [Float64(LatticeCore.k_point(mesh, i)[2]) for i in 1:M]
        Plots.scatter!(
            p, kxs, kys; mc=mesh_color, ms=mesh_size, markerstrokewidth=0, label="mesh"
        )
    end

    if show_high_symmetry
        hsp = Lattice2D.high_symmetry_points(lat)
        hxs = Float64[]
        hys = Float64[]
        labels = String[]
        for (name, k) in hsp
            push!(hxs, Float64(k[1]))
            push!(hys, Float64(k[2]))
            push!(labels, String(name))
        end
        Plots.scatter!(
            p, hxs, hys; mc=label_color, ms=5, markerstrokewidth=0, label="high-sym"
        )
        for (x, y, name) in zip(hxs, hys, labels)
            Plots.annotate!(p, x, y, Plots.text(" $name", 10, label_color, :left, :bottom))
        end
    end

    return p
end

end # module Lattice2DPlotsExt
