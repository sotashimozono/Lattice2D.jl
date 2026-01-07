using Lattice2D, Test, Random

@testset "DLA Tests" begin
    rng = MersenneTwister(42)
    
    @testset "DLA Cluster Creation" begin
        lat = build_lattice(Square, 30, 30)
        
        # Test default creation
        dla = dla_cluster(lat)
        @test dla isa DLACluster
        @test dla.seed_site == div(lat.N, 2)
        @test dla.n_particles == 1
        @test length(dla.occupied_sites) == 1
        @test dla.seed_site in dla.occupied_sites
        
        # Test custom seed
        dla_custom = dla_cluster(lat, 50)
        @test dla_custom.seed_site == 50
        @test 50 in dla_custom.occupied_sites
        
        # Test bounds
        @test_throws AssertionError dla_cluster(lat, 0)
        @test_throws AssertionError dla_cluster(lat, lat.N + 1)
    end
    
    @testset "Growing DLA" begin
        lat = build_lattice(Square, 30, 30)
        dla = dla_cluster(lat)
        
        # Add single particle
        success = add_particle!(dla, rng=rng)
        # Note: success might be false if particle escapes
        if success
            @test dla.n_particles == 2
            @test length(dla.occupied_sites) == 2
        end
        
        # Grow multiple particles
        initial_size = dla.n_particles
        n_added = grow!(dla, 10, rng=rng)
        @test n_added >= 0
        @test n_added <= 10
        @test dla.n_particles == initial_size + n_added
        @test length(dla.occupied_sites) == dla.n_particles
    end
    
    @testset "Cluster Properties" begin
        lat = build_lattice(Square, 40, 40)
        dla = dla_cluster(lat)
        grow!(dla, 30, rng=rng)
        
        # Test radius
        R = cluster_radius(dla)
        @test R >= 0.0
        @test isfinite(R)
        
        # Test mass center
        cm = mass_center(dla)
        @test length(cm) == 2
        @test all(isfinite.(cm))
        
        # Seed should be near center for small clusters
        seed_pos = lat.positions[dla.seed_site]
        dist_to_cm = sqrt(sum((seed_pos - cm) .^ 2))
        @test dist_to_cm < R * 2  # Reasonable check
    end
    
    @testset "Radial Density" begin
        lat = build_lattice(Square, 40, 40)
        dla = dla_cluster(lat)
        grow!(dla, 50, rng=rng)
        
        bin_edges, counts = radial_density(dla, 10)
        @test length(bin_edges) == 11  # n_bins + 1
        @test length(counts) == 10
        @test sum(counts) == dla.n_particles
        @test all(counts .>= 0)
    end
    
    @testset "Fractal Dimension" begin
        lat = build_lattice(Square, 50, 50)
        dla = dla_cluster(lat)
        
        # Small cluster - may return NaN
        D_small = fractal_dimension(dla)
        @test isnan(D_small) || (0.0 < D_small < 3.0)
        
        # Larger cluster
        grow!(dla, 100, rng=rng)
        if dla.n_particles >= 10
            D = fractal_dimension(dla)
            # DLA in 2D typically has D ≈ 1.7
            # But we just check it's reasonable
            @test isnan(D) || (1.0 < D < 2.5)
        end
    end
    
    @testset "Different Lattices" begin
        for LatticeType in [Square, Triangular, Honeycomb]
            lat = build_lattice(LatticeType, 30, 30)
            dla = dla_cluster(lat)
            
            @test dla.n_particles == 1
            
            # Try to grow
            n_added = grow!(dla, 20, rng=rng)
            @test n_added >= 0
            @test dla.n_particles >= 1
        end
    end
    
    @testset "Radius Updates" begin
        lat = build_lattice(Square, 40, 40)
        dla = dla_cluster(lat)
        
        initial_spawn = dla.spawn_radius
        initial_kill = dla.kill_radius
        
        # Grow cluster
        grow!(dla, 50, rng=rng)
        
        # Radii should have grown with cluster
        @test dla.spawn_radius >= initial_spawn
        @test dla.kill_radius >= initial_kill
        @test dla.kill_radius > dla.spawn_radius
    end
end
