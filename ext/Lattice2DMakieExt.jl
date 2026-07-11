module Lattice2DMakieExt

using Lattice2D
using LatticeCore
using Makie

"""
    Lattice2DMakieExt

Makie backend for `Lattice2D`, loaded automatically once `Makie` is in scope
(e.g. `using CairoMakie` or `using GLMakie`). It implements the three
`makie_*` entry points declared in the core module:

- [`makie_lattice`](@ref) — geometry (sites + bonds), sublattice colouring,
  bond highlighting;
- [`makie_state`](@ref) — per-site scalar field as a colour-mapped scatter,
  optional XY-spin arrows;
- [`makie_structure_factor`](@ref) — static structure factor `S(k)` heatmap.

Each returns a `Makie.Figure`; the active Makie backend (CairoMakie /
GLMakie) decides how it is rendered or saved. The `makie_` prefix keeps these
disjoint from the `Plots` backend's `plot_*` methods, so both can be loaded
at once.
"""
Lattice2DMakieExt

# ---- shared helpers --------------------------------------------------

@inline _xy(lat, i) = (p=LatticeCore.position(lat, i); (Float64(p[1]), Float64(p[2])))

function _site_coords(lat)
    N = LatticeCore.num_sites(lat)
    xs = Vector{Float64}(undef, N)
    ys = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        xs[i], ys[i] = _xy(lat, i)
    end
    return xs, ys
end

# Bond segments drawn from site `i` along the per-bond *wrapped* displacement
# `b.vector`, so periodic samples show no boundary-crossing long lines.
function _bond_segments(lat, indices)
    seg = Point2f[]
    bs = LatticeCore.bonds(lat)
    for k in indices
        b = bs[k]
        xi, yi = _xy(lat, b.i)
        push!(seg, Point2f(xi, yi))
        push!(seg, Point2f(xi + Float64(b.vector[1]), yi + Float64(b.vector[2])))
    end
    return seg
end

# ---- makie_lattice ---------------------------------------------------

function Lattice2D.makie_lattice(
    lat::Lattice2D.Lattice;
    colorby::Symbol=:sublattice,
    highlight_bonds=nothing,
    markersize::Real=12,
    show_sites::Bool=true,
    figure=(;),
    axis=(;),
)
    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), axis...)

    nb = length(LatticeCore.bonds(lat))
    linesegments!(ax, _bond_segments(lat, 1:nb); color=(:gray, 0.6), linewidth=1.0)

    if highlight_bonds !== nothing
        linesegments!(ax, _bond_segments(lat, highlight_bonds); color=:red, linewidth=2.5)
    end

    if show_sites
        xs, ys = _site_coords(lat)
        if colorby === :sublattice && LatticeCore.num_sublattices(lat) > 1
            cols = [LatticeCore.sublattice(lat, i) for i in eachindex(xs)]
            scatter!(ax, xs, ys; color=cols, colormap=:tab10, markersize=markersize)
        else
            scatter!(ax, xs, ys; color=:steelblue, markersize=markersize)
        end
    end
    return fig
end

# ---- makie_state -----------------------------------------------------

function Lattice2D.makie_state(
    lat::Lattice2D.Lattice,
    state::AbstractVector;
    colormap=:RdBu,
    arrows::Bool=false,
    markersize::Real=15,
    figure=(;),
    axis=(;),
)
    N = LatticeCore.num_sites(lat)
    length(state) == N || throw(
        DimensionMismatch("state has length $(length(state)) but lattice has $N sites")
    )
    xs, ys = _site_coords(lat)

    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), axis...)
    sc = scatter!(
        ax, xs, ys; color=Float64.(state), colormap=colormap, markersize=markersize
    )
    Colorbar(fig[1, 2], sc)

    if arrows
        seg = Point2f[]
        @inbounds for i in 1:N
            θ = Float64(state[i])
            dx, dy = 0.4 * cos(θ), 0.4 * sin(θ)
            push!(seg, Point2f(xs[i] - dx, ys[i] - dy))
            push!(seg, Point2f(xs[i] + dx, ys[i] + dy))
        end
        linesegments!(ax, seg; color=:black, linewidth=1.5)
    end
    return fig
end

# ---- makie_structure_factor -----------------------------------------

# S(k) = |Σ_j state_j exp(-i k·r_j)|² / N on a resolution×resolution k-grid.
# Split out from the plot so the numerics can be tested without rendering.
function _structure_factor_grid(lat, state::AbstractVector; k_range, resolution::Int)
    N = LatticeCore.num_sites(lat)
    length(state) == N || throw(
        DimensionMismatch("state has length $(length(state)) but lattice has $N sites")
    )
    resolution ≥ 1 || throw(ArgumentError("resolution must be ≥ 1, got $resolution"))
    xs, ys = _site_coords(lat)
    st = ComplexF64.(state)
    ks = range(Float64(k_range[1]), Float64(k_range[2]); length=resolution)
    S = Matrix{Float64}(undef, resolution, resolution)
    @inbounds for a in 1:resolution
        kx = ks[a]
        for b in 1:resolution
            ky = ks[b]
            acc = zero(ComplexF64)
            for j in 1:N
                acc += st[j] * cis(-(kx * xs[j] + ky * ys[j]))
            end
            S[a, b] = abs2(acc) / N
        end
    end
    return ks, S
end

function Lattice2D.makie_structure_factor(
    lat::Lattice2D.Lattice,
    state::AbstractVector;
    k_range=(-π, π),
    resolution::Int=200,
    colormap=:viridis,
    figure=(;),
    axis=(;),
)
    ks, S = _structure_factor_grid(lat, state; k_range=k_range, resolution=resolution)

    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="kₓ", ylabel="k_y", axis...)
    hm = heatmap!(ax, ks, ks, S; colormap=colormap)
    Colorbar(fig[1, 2], hm)
    return fig
end

end # module Lattice2DMakieExt
