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

end # module Lattice2DPlotsExt
