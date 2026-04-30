using Plots
using LatticeCore

# `runtests.jl` already sets `ENV["GKSwstype"] = "100"` so GR runs in
# a headless mode. We only build the figure objects and inspect their
# series — nothing is written to disk — but the env var still applies
# in case any kwarg path triggers a backend init.

@testset "Lattice2DPlotsExt — plot_state" begin
    @testset "Extension is loaded once Plots is in scope" begin
        ext = Base.get_extension(Lattice2D, :Lattice2DPlotsExt)
        @test ext !== nothing
    end

    @testset "Continuous state on Square (PBC) returns a Plots.Plot" begin
        lat = build_lattice(Square, 4, 4)
        N = num_sites(lat)
        state = [Float64(i) for i in 1:N]

        p = plot_state(lat, state)
        @test p isa Plots.Plot

        # Bonds (lines) + scatter series should both be present.
        nseries = length(p.series_list)
        @test nseries >= 2

        # The continuous code path sets `marker_z` on the scatter series.
        scatter_series = filter(s -> s[:seriestype] === :scatter, p.series_list)
        @test !isempty(scatter_series)
        @test any(!isnothing(s[:marker_z]) for s in scatter_series)
    end

    @testset "Bool state on Honeycomb uses discrete palette" begin
        lat = build_lattice(Honeycomb, 3, 3)
        N = num_sites(lat)
        state = [iseven(i) for i in 1:N]

        p = plot_state(lat, state; marker_size=8, title="bool")
        @test p isa Plots.Plot

        # Discrete path: one scatter series per unique label, none of
        # which carries a marker_z (those would force a colorbar).
        scatter_series = filter(s -> s[:seriestype] === :scatter, p.series_list)
        @test length(scatter_series) == length(unique(state))
        @test all(isnothing(s[:marker_z]) for s in scatter_series)
    end

    @testset "Integer (small distinct count) state is discrete" begin
        lat = build_lattice(Triangular, 3, 3)
        N = num_sites(lat)
        # 3-state clock model: labels 1, 2, 3.
        state = [mod1(i, 3) for i in 1:N]

        p = plot_state(lat, state)
        @test p isa Plots.Plot

        scatter_series = filter(s -> s[:seriestype] === :scatter, p.series_list)
        @test length(scatter_series) == 3
    end

    @testset "show_bonds=false drops the bond series" begin
        lat = build_lattice(Square, 3, 3)
        state = randn(num_sites(lat))

        p_with = plot_state(lat, state)
        p_without = plot_state(lat, state; show_bonds=false)

        @test length(p_without.series_list) < length(p_with.series_list)
    end

    @testset "Length mismatch throws DimensionMismatch" begin
        lat = build_lattice(Square, 3, 3)
        @test_throws DimensionMismatch plot_state(lat, randn(num_sites(lat) + 1))
    end

    @testset "Works on multi-sublattice topology (Kagome)" begin
        lat = build_lattice(Kagome, 3, 3)
        state = randn(num_sites(lat))
        p = plot_state(lat, state; colormap=:plasma)
        @test p isa Plots.Plot
    end
end
