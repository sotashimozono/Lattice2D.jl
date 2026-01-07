using Lattice2D, Test, Plots, Random

@testset "Visualization Tests" begin
    rng = MersenneTwister(42)
    
    # Disable plot display for testing
    ENV["GKSwstype"] = "100"
    
    @testset "Plot Percolation" begin
        lat = build_lattice(Square, 10, 10)
        perc_bond = bond_percolation(lat, 0.5, rng=rng)
        perc_site = site_percolation(lat, 0.6, rng=rng)
        
        # Test bond percolation plot
        p1 = plot_percolation(perc_bond)
        @test p1 isa Plots.Plot
        
        p2 = plot_percolation(perc_bond, show_clusters=false)
        @test p2 isa Plots.Plot
        
        # Test site percolation plot
        p3 = plot_percolation(perc_site)
        @test p3 isa Plots.Plot
        
        p4 = plot_percolation(perc_site, show_clusters=false, markersize=5)
        @test p4 isa Plots.Plot
    end
    
    @testset "Plot Random Walk" begin
        lat = build_lattice(Triangular, 15, 15)
        rw = random_walk(lat)
        walk!(rw, 100, rng=rng)
        
        p1 = plot_random_walk(rw)
        @test p1 isa Plots.Plot
        
        p2 = plot_random_walk(rw, show_lattice=false)
        @test p2 isa Plots.Plot
        
        p3 = plot_random_walk(rw, linewidth=3, markersize=6)
        @test p3 isa Plots.Plot
    end
    
    @testset "Plot DLA" begin
        lat = build_lattice(Square, 30, 30)
        dla = dla_cluster(lat)
        grow!(dla, 50, rng=rng)
        
        p1 = plot_dla(dla)
        @test p1 isa Plots.Plot
        
        p2 = plot_dla(dla, color_by_distance=false)
        @test p2 isa Plots.Plot
        
        p3 = plot_dla(dla, markersize=5)
        @test p3 isa Plots.Plot
    end
    
    @testset "Plot Spanning Tree" begin
        lat = build_lattice(Honeycomb, 8, 8)
        tree = wilson_algorithm(lat, rng=rng)
        
        p1 = plot_spanning_tree(tree)
        @test p1 isa Plots.Plot
        
        p2 = plot_spanning_tree(tree, color_by_depth=false)
        @test p2 isa Plots.Plot
        
        p3 = plot_spanning_tree(tree, linewidth=3, markersize=4)
        @test p3 isa Plots.Plot
    end
    
    @testset "Plot Cluster Size Distribution" begin
        lat = build_lattice(Square, 15, 15)
        perc = bond_percolation(lat, 0.5, rng=rng)
        
        p1 = plot_cluster_size_distribution(perc)
        @test p1 isa Plots.Plot
        
        p2 = plot_cluster_size_distribution(perc, logscale=false)
        @test p2 isa Plots.Plot
    end
    
    @testset "Plot MSD vs Time" begin
        lat = build_lattice(Square, 20, 20)
        
        p = plot_msd_vs_time(lat, 50, 5, rng=rng)
        @test p isa Plots.Plot
    end
    
    @testset "Plot Radial Density" begin
        lat = build_lattice(Square, 30, 30)
        dla = dla_cluster(lat)
        grow!(dla, 50, rng=rng)
        
        p1 = plot_radial_density(dla)
        @test p1 isa Plots.Plot
        
        p2 = plot_radial_density(dla, n_bins=15)
        @test p2 isa Plots.Plot
    end
    
    @testset "Animation" begin
        # Test that animation function exists and can be called
        lat = build_lattice(Square, 10, 10)
        rw = random_walk(lat)
        walk!(rw, 20, rng=rng)
        
        # Create the animation in temp directory
        tmpfile = joinpath(mktempdir(), "test_walk.gif")
        animate_random_walk(rw, tmpfile, fps=10)
        
        # Check the file was created
        @test isfile(tmpfile)
        
        # Cleanup
        rm(tmpfile, force=true)
    end
    
    @testset "Different Lattice Types" begin
        for LatticeType in [Square, Triangular, Honeycomb, Kagome]
            lat = build_lattice(LatticeType, 8, 8)
            
            # Test percolation plot
            perc = bond_percolation(lat, 0.5, rng=rng)
            p1 = plot_percolation(perc)
            @test p1 isa Plots.Plot
            
            # Test random walk plot
            rw = random_walk(lat)
            walk!(rw, 30, rng=rng)
            p2 = plot_random_walk(rw)
            @test p2 isa Plots.Plot
            
            # Test spanning tree plot
            tree = wilson_algorithm(lat, rng=rng)
            p3 = plot_spanning_tree(tree)
            @test p3 isa Plots.Plot
        end
    end
end
