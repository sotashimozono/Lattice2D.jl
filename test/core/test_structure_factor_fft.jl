using FFTW
using LatticeCore
using LinearAlgebra
using Random
using StaticArrays

# Build state values that realise a plane wave at a chosen k-target.
# For a Bravais multi-sublattice lattice with positions r_i,
# `s_i = exp(i k_target · r_i)` makes
#     Σ_i s_i e^{-i k r_i} = Σ_i e^{i (k_target - k) · r_i}
# vanish off `k_target` (mod reciprocal-lattice vectors) and equal
# `Nsites` at `k_target`, so `S(k_target) = Nsites`.
function _plane_wave_state(lat, k_target)
    return [cis(dot(k_target, position(lat, i))) for i in 1:num_sites(lat)]
end

@testset "structure_factor FFT path (Lattice2DFFTWExt)" begin
    Random.seed!(20260429)

    @testset "FFT vs naive on every Bravais topology" begin
        # Random states, MP mesh from `reciprocal_lattice`.
        for Topo in (
            Square, Triangular, Honeycomb, Kagome, Lieb, UnionJack, Dice, ShastrySutherland
        )
            for (Lx, Ly) in ((4, 4), (6, 4))
                lat = build_lattice(Topo, Lx, Ly)
                ml = reciprocal_lattice(lat)
                for trial in 1:2
                    state = randn(num_sites(lat))
                    fast = structure_factor(lat, state, ml)
                    slow = LatticeCore._structure_factor_naive(lat, state, ml)
                    @test fast ≈ slow rtol = 1e-9
                end
            end
        end
    end

    @testset "FFT vs naive on Γ-centred mesh" begin
        # `gamma_centered` exercises the no-prephase branch in the FFT
        # path (frac0 == 0). Cover one single-sublattice topology and
        # one multi-sublattice topology.
        for Topo in (Triangular, Kagome)
            lat = build_lattice(Topo, 6, 6)
            A = basis_vectors(lat)
            B = SMatrix{2,2,Float64}(2π * inv(transpose(A)))
            ml = LatticeCore.gamma_centered(B, (lat.Lx, lat.Ly))
            state = randn(num_sites(lat))
            fast = structure_factor(lat, state, ml)
            slow = LatticeCore._structure_factor_naive(lat, state, ml)
            @test fast ≈ slow rtol = 1e-9
        end
    end

    # ---- Analytical k-space checks ---------------------------------
    #
    # The defining identity `S(k) = (1/N) |Σ_i s_i e^{-i k r_i}|²`
    # gives closed-form values for a few canonical states. We
    # cross-check those against both the naive reference loop and the
    # FFT-accelerated path:
    #
    # - Uniform state (`s_i ≡ 1`): peaks at k = 0 with `S(Γ) = Nsites`,
    #   and `S(k) = 0` at every other Bragg-mesh k inside the BZ
    #   (their cell sum vanishes).
    #
    # - Plane wave at `k_target` (`s_i = exp(i k_target · r_i)`): peaks
    #   at `k = k_target` with `S(k_target) = Nsites`. Useful for
    #   multi-sublattice topologies where the basis-position phase
    #   factor `exp(-i k · d_α)` shows up explicitly.

    @testset "Triangular: uniform state peaks at Γ" begin
        lat = build_lattice(Triangular, 6, 6)
        A = basis_vectors(lat)
        B = SMatrix{2,2,Float64}(2π * inv(transpose(A)))
        ml = LatticeCore.gamma_centered(B, (lat.Lx, lat.Ly))
        N = num_sites(lat)
        state = ones(Float64, N)
        S = structure_factor(lat, state, ml)
        @test S[1] ≈ Float64(N)
        @test maximum(@view S[2:end]) < 1e-9
    end

    @testset "Triangular: stripe (-1)^(cx+cy) peaks at the M-point" begin
        # A stripe of cell-index parity is a plane wave at
        # `k_M = (b₁ + b₂) / 2`, which on the Γ-centred 6×6 mesh
        # lands exactly on a mesh point; S at that k equals Nsites.
        lat = build_lattice(Triangular, 6, 6)
        A = basis_vectors(lat)
        B = SMatrix{2,2,Float64}(2π * inv(transpose(A)))
        ml = LatticeCore.gamma_centered(B, (lat.Lx, lat.Ly))
        N = num_sites(lat)
        state = Vector{Float64}(undef, N)
        for i in 1:N
            coord = LatticeCore.lattice_coord(lat.indexing, (lat.Lx, lat.Ly), 1, i)
            cx, cy = coord.cell
            state[i] = (-1.0)^(cx + cy)
        end
        S = structure_factor(lat, state, ml)
        @test maximum(S) ≈ Float64(N) rtol = 1e-10
        # And the peak must sit at k = b₁/2 + b₂/2 (M-point).
        peak_k = ml.k_points[argmax(S)]
        b1 = SVector{2}(B[:, 1])
        b2 = SVector{2}(B[:, 2])
        @test peak_k ≈ (b1 + b2) / 2 rtol = 1e-10
    end

    @testset "Honeycomb: uniform state peaks at Γ with value Nsites" begin
        lat = build_lattice(Honeycomb, 4, 4)
        A = basis_vectors(lat)
        B = SMatrix{2,2,Float64}(2π * inv(transpose(A)))
        ml = LatticeCore.gamma_centered(B, (lat.Lx, lat.Ly))
        N = num_sites(lat)
        state = ones(Float64, N)
        S = structure_factor(lat, state, ml)
        @test S[1] ≈ Float64(N)
        @test maximum(@view S[2:end]) < 1e-9
    end

    @testset "Honeycomb: plane wave at a chosen mesh k" begin
        # Build s_i = exp(i k_target · r_i) and verify the peak lands
        # exactly at k_target. This exercises the basis-position phase
        # factor `exp(-i k · d_α)` in the multi-sublattice FFT path.
        lat = build_lattice(Honeycomb, 4, 4)
        ml = reciprocal_lattice(lat)
        N = num_sites(lat)
        for target_idx in (1, 4, 7, 13)
            k_target = ml.k_points[target_idx]
            state = _plane_wave_state(lat, k_target)
            S = structure_factor(lat, state, ml)
            @test argmax(S) == target_idx
            @test S[target_idx] ≈ Float64(N) rtol = 1e-10
            # All other k-points must vanish to FP noise.
            S_off = copy(S);
            S_off[target_idx] = 0.0
            @test maximum(S_off) < 1e-9
        end
    end

    @testset "Kagome: uniform state peaks at Γ with value Nsites" begin
        lat = build_lattice(Kagome, 4, 4)
        A = basis_vectors(lat)
        B = SMatrix{2,2,Float64}(2π * inv(transpose(A)))
        ml = LatticeCore.gamma_centered(B, (lat.Lx, lat.Ly))
        N = num_sites(lat)
        state = ones(Float64, N)
        S = structure_factor(lat, state, ml)
        @test S[1] ≈ Float64(N)
        @test maximum(@view S[2:end]) < 1e-9
    end

    @testset "Kagome: plane wave at a chosen mesh k" begin
        lat = build_lattice(Kagome, 4, 4)
        ml = reciprocal_lattice(lat)
        N = num_sites(lat)
        for target_idx in (1, 5, 11)
            k_target = ml.k_points[target_idx]
            state = _plane_wave_state(lat, k_target)
            S = structure_factor(lat, state, ml)
            @test argmax(S) == target_idx
            @test S[target_idx] ≈ Float64(N) rtol = 1e-10
        end
    end

    # ---- Fallback paths -------------------------------------------

    @testset "FFT path falls back on non-RowMajor indexing" begin
        # ColMajor / Snake state layouts don't match the closed-form
        # `(cy-1)*Lx + (cx-1)*nsub + α` site index, so the multi-sub
        # eligibility check should reject them and route through the
        # naive helper. Numerical equality is the contract.
        for indexing in (ColMajor(), Snake())
            lat = build_lattice(Honeycomb, 4, 4; indexing=indexing)
            ml = reciprocal_lattice(lat)
            state = randn(num_sites(lat))
            fast = structure_factor(lat, state, ml)
            slow = LatticeCore._structure_factor_naive(lat, state, ml)
            @test fast ≈ slow rtol = 1e-12
        end

        # Single-sublattice opt-in (`Triangular`) likewise refuses
        # non-RowMajor.
        for indexing in (ColMajor(), Snake())
            lat = build_lattice(Triangular, 4, 4; indexing=indexing)
            @test LatticeCore._has_known_grid(lat) == false
        end
    end

    @testset "Single-sublattice opt-in matches LC FFTW path" begin
        # For the single-sublattice topologies (`Square`,
        # `Triangular`) we register `LatticeCore._has_known_grid` /
        # `LatticeCore._reshape_state` and let `LatticeCoreFFTWExt`
        # drive the FFT. Lock that opt-in down.
        for Topo in (Square, Triangular)
            lat = build_lattice(Topo, 4, 4)
            @test LatticeCore._has_known_grid(lat) == true
            v = collect(1:num_sites(lat))
            grid = LatticeCore._reshape_state(lat, v, (lat.Lx, lat.Ly))
            @test size(grid) == (lat.Lx, lat.Ly)
        end
    end
end
