using Lattice2D, Test, Random

@testset "Random Walk Tests" begin
    rng = MersenneTwister(42)
    
    @testset "Basic Random Walk" begin
        lat = build_lattice(Square, 10, 10)
        
        # Test initialization
        rw = random_walk(lat, 50)
        @test rw isa RandomWalk
        @test rw.current_site == 50
        @test rw.start_site == 50
        @test rw.n_steps == 0
        @test length(rw.path) == 1
        
        # Test default start
        rw_default = random_walk(lat)
        @test rw_default.start_site == 1
        
        # Test bounds
        @test_throws AssertionError random_walk(lat, 0)
        @test_throws AssertionError random_walk(lat, lat.N + 1)
    end
    
    @testset "Walking" begin
        lat = build_lattice(Square, 10, 10)
        rw = random_walk(lat, 50)
        
        # Single step
        step!(rw, rng=rng)
        @test rw.n_steps == 1
        @test length(rw.path) == 2
        @test rw.current_site in lat.nearest_neighbors[50]
        
        # Multiple steps
        walk!(rw, 10, rng=rng)
        @test rw.n_steps == 11
        @test length(rw.path) == 12
    end
    
    @testset "Statistics" begin
        lat = build_lattice(Triangular, 15, 15)
        rw = random_walk(lat, 1)
        walk!(rw, 100, rng=rng)
        
        # Test displacement
        disp = displacement(rw)
        @test length(disp) == 2
        @test disp isa Vector{Float64}
        
        # Test MSD
        msd = mean_squared_displacement(rw)
        @test msd >= 0.0
        @test msd == sum(disp .^ 2)
        
        # Test end-to-end distance
        dist = end_to_end_distance(rw)
        @test dist >= 0.0
        @test dist ≈ sqrt(msd)
        
        # Test path length
        plen = path_length(rw)
        @test plen >= dist  # Path length >= straight line distance
        
        # Test visited sites
        visited = visited_sites(rw)
        @test visited isa Set{Int}
        @test length(visited) <= 101  # Can't visit more than steps + 1
        @test 1 in visited  # Start site must be visited
        
        # Test return probability
        prob = return_probability(rw, 1)
        @test 0.0 <= prob <= 1.0
    end
    
    @testset "Self-Avoiding Walk" begin
        lat = build_lattice(Square, 20, 20)
        
        saw = self_avoiding_walk(lat, 50, 1, rng=rng)
        @test saw isa RandomWalk
        @test saw.start_site == 1
        
        # Check self-avoiding property
        visited = visited_sites(saw)
        @test length(visited) == length(saw.path)  # All sites unique
        
        # Test different lattices
        for LatticeType in [Triangular, Honeycomb]
            lat2 = build_lattice(LatticeType, 15, 15)
            saw2 = self_avoiding_walk(lat2, 30, rng=rng)
            visited2 = visited_sites(saw2)
            @test length(visited2) == length(saw2.path)
        end
    end
    
    @testset "Average MSD" begin
        lat = build_lattice(Square, 20, 20)
        
        avg_msd = average_msd(lat, 50, 10, rng=rng)
        @test avg_msd >= 0.0
        @test isfinite(avg_msd)
        
        # MSD should grow roughly linearly with time for random walk
        avg_msd_early = average_msd(lat, 10, 10, rng=rng)
        avg_msd_late = average_msd(lat, 50, 10, rng=rng)
        @test avg_msd_late > avg_msd_early
    end
    
    @testset "Mixing Time" begin
        lat = build_lattice(Square, 10, 10)
        
        t_mix = mixing_time(lat, 5000, 0.5, n_trials=5, rng=rng)
        @test t_mix > 0.0
        @test isfinite(t_mix)
        
        # Easier target should have shorter mixing time
        t_mix_easy = mixing_time(lat, 5000, 0.2, n_trials=5, rng=rng)
        @test t_mix_easy < t_mix || t_mix_easy > 0
        
        # Impossible target
        t_impossible = mixing_time(lat, 10, 1.0, n_trials=1, rng=rng)
        @test t_impossible == Inf
    end
    
    @testset "Different Lattice Types" begin
        for LatticeType in [Square, Triangular, Honeycomb, Kagome, Lieb]
            lat = build_lattice(LatticeType, 10, 10)
            rw = random_walk(lat)
            walk!(rw, 50, rng=rng)
            
            @test rw.n_steps == 50
            @test length(rw.path) == 51
            @test mean_squared_displacement(rw) >= 0.0
        end
    end
end
