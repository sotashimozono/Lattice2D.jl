module Lattice2DFFTWExt

using Lattice2D
using LatticeCore
using FFTW
using StaticArrays
using LinearAlgebra

"""
    Lattice2DFFTWExt

Optional `Lattice2D` extension that opts the package's `Lattice{Topo,T}`
topologies into the FFT-based fast path of
`structure_factor(lat, state, ml)`.

Loaded automatically once the user issues `using FFTW` (with the
`LatticeCore` FFT extension already wired up by the same `using FFTW`
upstream).

Two orthogonal mechanisms are used here:

- **Single-sublattice topologies** (`Square`, `Triangular`): they
  share the row-major site-index ↔ Bravais-cell layout used by
  `LineLattice` / `SimpleSquareLattice` upstream, so we register the
  upstream opt-in hooks `LatticeCore._has_known_grid` /
  `LatticeCore._reshape_state` and let `LatticeCoreFFTWExt` drive the
  FFT directly.

- **Multi-sublattice topologies** (`Honeycomb`, `Kagome`, `Lieb`,
  `UnionJack`, `Dice`, `ShastrySutherland`): the natural state layout
  is a 3-axis grid `(nsub, Lx, Ly)`, which doesn't match the upstream
  hook contract (a 2D reshape sized to `ml.mesh`). Instead we
  override `LatticeCore._structure_factor_fast(::HasReciprocal,
  ::Lattice{Topo,T}, ...)` ourselves: one column-major `reshape`
  splits the state by sublattice, each per-sublattice plane is
  FFT'd with the same prephase trick used upstream, and the
  per-sublattice contributions are combined coherently in k-space
  with the basis-position phase factor `exp(-i k · d_α)` before
  taking `|·|² / N`.

The FFT path requires:

- both axes periodic (already enforced by `reciprocal_support(lat)`
  returning `HasReciprocal()` only when neither axis is open),
- `RowMajor` indexing (the default), since the closed-form
  `state[i]` ↔ `(cx, cy, α)` decomposition is what the FFT consumes,
- and `ml::PeriodicMomentumLattice` with `ml.mesh == (Lx, Ly)`.

If any precondition fails we fall back to the naive O(N · M) loop
through `LatticeCore._structure_factor_naive`, so loading the
extension never changes results — only performance.
"""
Lattice2DFFTWExt

# ---- Single-sublattice opt-in ---------------------------------------
#
# `Square` and `Triangular` have one site per unit cell, so the
# state-vector layout matches the LatticeCore reference contract
# (`reshape(state, (Lx, Ly))`). Register the upstream opt-in hooks
# directly — the FFT path for these topologies then reuses
# `LatticeCoreFFTWExt._structure_factor_fft` unchanged.

function LatticeCore._has_known_grid(lat::Lattice2D.Lattice{Lattice2D.Square})
    return _is_rowmajor_periodic(lat)
end
function LatticeCore._has_known_grid(lat::Lattice2D.Lattice{Lattice2D.Triangular})
    return _is_rowmajor_periodic(lat)
end

function LatticeCore._reshape_state(::Lattice2D.Lattice{Lattice2D.Square}, state, dims)
    return reshape(state, dims)
end
function LatticeCore._reshape_state(::Lattice2D.Lattice{Lattice2D.Triangular}, state, dims)
    return reshape(state, dims)
end

# ---- Multi-sublattice override --------------------------------------
#
# For these topologies the upstream 2D-reshape hook isn't expressive
# enough; we override the trait-dispatched fast path directly. The
# method below is strictly more specific than
# `LatticeCoreFFTWExt`'s own `(::HasReciprocal, ::AbstractLattice, ...)`
# method so dispatch picks ours for `Lattice{Topo,T}` whenever this
# extension is loaded.
const _MultiSublatticeTopo = Union{
    Lattice2D.Honeycomb,
    Lattice2D.Kagome,
    Lattice2D.Lieb,
    Lattice2D.UnionJack,
    Lattice2D.Dice,
    Lattice2D.ShastrySutherland,
}

function LatticeCore._structure_factor_fast(
    ::LatticeCore.HasReciprocal,
    lat::Lattice2D.Lattice{Topo,T},
    state::AbstractVector,
    ml::LatticeCore.AbstractMomentumLattice,
) where {Topo<:_MultiSublatticeTopo,T}
    if _multi_sub_eligible(lat, ml)
        return _structure_factor_fft_multisub(lat, state, ml)
    else
        return LatticeCore._structure_factor_naive(lat, state, ml)
    end
end

# ---- Eligibility helpers --------------------------------------------

# RowMajor + fully periodic. Cylinder / OBC already gets filtered out
# by `reciprocal_support(lat) == NoReciprocal()` upstream, but we
# guard explicitly so the fast path is robust to future trait
# refinements.
function _is_rowmajor_periodic(lat::Lattice2D.Lattice)
    lat.indexing isa LatticeCore.RowMajor || return false
    return all(ax isa LatticeCore.PeriodicAxis for ax in lat.boundary.axes)
end

function _multi_sub_eligible(
    lat::Lattice2D.Lattice, ml::LatticeCore.AbstractMomentumLattice
)
    ml isa LatticeCore.PeriodicMomentumLattice || return false
    _is_rowmajor_periodic(lat) || return false
    st = LatticeCore.size_trait(lat)
    st isa LatticeCore.FiniteSize || return false
    return st.dims == ml.mesh
end

# ---- FFT core for multi-sublattice Bravais lattices -----------------
#
# For a `Lattice{Topo,T}` with `nsub` sublattices and RowMajor
# indexing,
#
#     site_index(cx, cy, α) = ((cy - 1) * Lx + (cx - 1)) * nsub + α,
#
# so `reshape(state, (nsub, Lx, Ly))` gives a 3-axis array
# `M[α, cx, cy]` whose `α`-th slice is the per-sublattice "image"
# of the state. The position of every site decomposes as
#
#     r_i = (cx - 1) a₁ + (cy - 1) a₂ + d_α,
#
# with `a₁`, `a₂` the Bravais primitive vectors and `d_α` the
# sublattice offset inside the unit cell. Substituting into the
# defining sum,
#
#     Σ_i s_i e^{-i k · r_i}
#       = Σ_α e^{-i k · d_α}
#           · Σ_{cx, cy} M[α, cx, cy] e^{-i k · ((cx-1) a₁ + (cy-1) a₂)}.
#
# The cell-index sum is exactly a 2D DFT of the `α`-th slice once we
# express `k` on the reciprocal mesh `k = B · frac` with
# `B^T A = 2π I`: then `k · a_d = 2π frac_d`, and the prephase
# trick used upstream (multiply by `cis(-2π Σ_d frac0_d (n_d - 1))`)
# absorbs the mesh offset into the input array.
function _structure_factor_fft_multisub(
    lat::Lattice2D.Lattice{Topo,T},
    state::AbstractVector,
    ml::LatticeCore.PeriodicMomentumLattice{2,Tk},
) where {Topo,T,Tk}
    Lx, Ly = lat.Lx, lat.Ly
    dims = (Lx, Ly)
    nsub = LatticeCore.num_sublattices(lat)
    Nsites = Lx * Ly * nsub
    Nk = Lx * Ly

    # `state` is laid out site-by-site with sublattice innermost. The
    # column-major reshape onto `(nsub, Lx, Ly)` therefore satisfies
    # `grid[α, cx, cy] = state[((cy - 1) * Lx + (cx - 1)) * nsub + α]`.
    grid = reshape(state, (nsub, Lx, Ly))

    B = LatticeCore.reciprocal_basis(ml)
    # Fractional offset of the mesh (per axis, scaled by 1/N_d).
    # `idx = 1` corresponds to `n = 0` per axis, so `B \ k_point(ml, 1)`
    # is exactly `offset / N` along each axis.
    frac0 = SVector{2,Float64}(B \ LatticeCore.k_point(ml, 1))

    # Sublattice offsets in real space, as `SVector{2}` for cheap
    # `dot(k, d_α)` later.
    uc = Lattice2D.get_unit_cell(Topo)
    sub_offsets = SVector{2,Float64}[
        SVector{2,Float64}(p[1], p[2]) for p in uc.sublattice_positions
    ]

    # Accumulator in k-space, shaped to match the (Lx, Ly) k-grid.
    acc = zeros(ComplexF64, dims)

    # Workspace for the per-sublattice prephased grid; reused across
    # sublattices to avoid an O(nsub) allocation churn.
    plane = Array{ComplexF64,2}(undef, dims)

    @inbounds for α in 1:nsub
        # Copy slice `grid[α, :, :]` into `plane`, applying the mesh
        # prephase. When the offset is exactly zero (Γ-centred mesh)
        # we skip the multiplication entirely.
        if iszero(frac0[1]) && iszero(frac0[2])
            for cy in 1:Ly, cx in 1:Lx
                plane[cx, cy] = ComplexF64(grid[α, cx, cy])
            end
        else
            for cy in 1:Ly, cx in 1:Lx
                ph = frac0[1] * (cx - 1) + frac0[2] * (cy - 1)
                plane[cx, cy] = ComplexF64(grid[α, cx, cy]) * cis(-2π * ph)
            end
        end

        F = fft(plane)

        # Multiply by the sublattice basis-position phase
        # `e^{-i k · d_α}` and accumulate into `acc`. `k` at index
        # `(idx_x, idx_y)` is `B * frac`, with
        # `frac_d = (idx_d - 1 + offset_d) / N_d`.
        d_α = sub_offsets[α]
        for idx_y in 1:Ly, idx_x in 1:Lx
            frac = SVector{2,Float64}(
                (idx_x - 1) / Lx + frac0[1], (idx_y - 1) / Ly + frac0[2]
            )
            k = B * frac
            phase = cis(-dot(k, d_α))
            acc[idx_x, idx_y] += F[idx_x, idx_y] * phase
        end
    end

    # Normalisation matches the naive `S(k) = |Σ s_i e^{-i k·r_i}|² / N`
    # convention, where `N` is the total number of sites (not k-points).
    out = Vector{Float64}(undef, Nk)
    @inbounds for i in eachindex(acc)
        out[i] = abs2(acc[i]) / Nsites
    end
    return out
end

end # module Lattice2DFFTWExt
