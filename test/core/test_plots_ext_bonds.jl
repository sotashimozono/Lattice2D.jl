using Plots
using LatticeCore

# `runtests.jl` already sets `ENV["GKSwstype"] = "100"` so GR runs in
# a headless mode. We only build the figure objects and inspect them
# — nothing is written to disk — but the env var still applies in
# case any kwarg path triggers a backend init.

@testset "Lattice2DPlotsExt — plot_bonds" begin
    @testset "Extension is loaded once Plots is in scope" begin
        ext = Base.get_extension(Lattice2D, :Lattice2DPlotsExt)
        @test ext !== nothing
    end

    @testset "Standalone plot_bonds returns a Plots.Plot" begin
        # ShastrySutherland has two distinct `Bond.type` tags, so
        # `color_by=:type` exercises the multi-group code path.
        lat = build_lattice(ShastrySutherland, 4, 4)

        p = plot_bonds(lat)
        @test p isa Plots.Plot

        p2 = plot_bonds(lat; color_by=:type, with_sites=true)
        @test p2 isa Plots.Plot

        p3 = plot_bonds(lat; color_by=:direction)
        @test p3 isa Plots.Plot
    end

    @testset "Filtering by bond_types" begin
        lat = build_lattice(ShastrySutherland, 4, 4)

        # Single-symbol selector.
        p = plot_bonds(lat; bond_types=:type_1)
        @test p isa Plots.Plot

        # Vector selector.
        p2 = plot_bonds(lat; bond_types=[:type_1, :type_2])
        @test p2 isa Plots.Plot

        # Tuple selector.
        p3 = plot_bonds(lat; bond_types=(:type_1,))
        @test p3 isa Plots.Plot

        # Empty selection — exercise the empty-bond fast return.
        p4 = plot_bonds(lat; bond_types=:does_not_exist)
        @test p4 isa Plots.Plot
    end

    @testset "Overlay plot_bonds! on plot_lattice" begin
        lat = build_lattice(Honeycomb, 4, 4)

        # `plot_lattice` is provided by LatticeCorePlotsExt; with
        # `show_bonds=false` we only get the site scatter, then
        # overlay the Lattice2D-aware bond renderer on top.
        p = plot_lattice(lat; show_bonds=false)
        plot_bonds!(p, lat; color_by=:direction, lw=2.0)
        @test p isa Plots.Plot
    end

    @testset "Argument validation" begin
        lat = build_lattice(Square, 3, 3)
        @test_throws ArgumentError plot_bonds(lat; color_by=:bogus)
        @test_throws ArgumentError plot_bonds(lat; bond_types=42)
    end
end
