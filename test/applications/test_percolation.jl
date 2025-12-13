using Lattice2D, Test, Random

@testset "Percolation Tests" begin
    # Set random seed for reproducibility
    rng = MersenneTwister(42)
    
    @testset "Bond Percolation" begin
        lat = build_lattice(Square, 10, 10)
        
        # Test basic creation
        perc = bond_percolation(lat, 0.5, rng=rng)
        @test perc isa BondPercolation
        @test 0 <= length(perc.occupied_bonds) <= length(lat.bonds)
        @test perc.probability == 0.5
        
        # Test edge cases
        perc_empty = bond_percolation(lat, 0.0, rng=rng)
        @test length(perc_empty.occupied_bonds) == 0
        @test number_of_clusters(perc_empty) == lat.N  # Each site is its own cluster
        
        perc_full = bond_percolation(lat, 1.0, rng=rng)
        @test length(perc_full.occupied_bonds) == length(lat.bonds)
        
        # Test probability bounds
        @test_throws AssertionError bond_percolation(lat, -0.1, rng=rng)
        @test_throws AssertionError bond_percolation(lat, 1.1, rng=rng)
    end
    
    @testset "Site Percolation" begin
        lat = build_lattice(Triangular, 10, 10)
        
        # Test basic creation
        perc = site_percolation(lat, 0.6, rng=rng)
        @test perc isa SitePercolation
        @test 0 <= length(perc.occupied_sites) <= lat.N
        @test perc.probability == 0.6
        
        # Test edge cases
        perc_empty = site_percolation(lat, 0.0, rng=rng)
        @test length(perc_empty.occupied_sites) == 0
        @test number_of_clusters(perc_empty) == 0
        
        perc_full = site_percolation(lat, 1.0, rng=rng)
        @test length(perc_full.occupied_sites) == lat.N
    end
    
    @testset "Cluster Analysis" begin
        lat = build_lattice(Square, 10, 10)
        perc = bond_percolation(lat, 0.5, rng=rng)
        
        # Test cluster functions
        n_clusters = number_of_clusters(perc)
        @test n_clusters > 0
        @test n_clusters <= lat.N
        
        max_size = largest_cluster_size(perc)
        @test max_size > 0
        @test max_size <= lat.N
        
        # Test cluster size distribution
        dist = cluster_size_distribution(perc)
        @test sum(values(dist)) == n_clusters
        @test sum(k * v for (k, v) in dist) == lat.N  # Total sites
        
        # Test percolation strength
        strength = percolation_strength(perc)
        @test 0.0 <= strength <= 1.0
        @test strength == max_size / lat.N
    end
    
    @testset "Spanning" begin
        lat = build_lattice(Square, 5, 5)
        
        # Full percolation should span
        perc_full = bond_percolation(lat, 1.0, rng=rng)
        @test spans_lattice(perc_full, direction=:horizontal) || spans_lattice(perc_full, direction=:vertical)
        
        # Empty percolation should not span
        perc_empty = bond_percolation(lat, 0.0, rng=rng)
        @test !spans_lattice(perc_empty, direction=:horizontal)
        @test !spans_lattice(perc_empty, direction=:vertical)
        @test !spans_lattice(perc_empty, direction=:both)
    end
    
    @testset "Different Lattices" begin
        # Test on different lattice types
        for LatticeType in [Square, Triangular, Honeycomb, Kagome]
            lat = build_lattice(LatticeType, 8, 8)
            perc_bond = bond_percolation(lat, 0.5, rng=rng)
            perc_site = site_percolation(lat, 0.5, rng=rng)
            
            @test number_of_clusters(perc_bond) > 0
            @test number_of_clusters(perc_site) >= 0
        end
    end
end
