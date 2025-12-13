using Lattice2D, Test, Random

@testset "Spanning Tree Tests" begin
    rng = MersenneTwister(42)
    
    @testset "Wilson Algorithm" begin
        lat = build_lattice(Square, 10, 10)
        
        # Test basic creation
        tree = wilson_algorithm(lat, 1, rng=rng)
        @test tree isa SpanningTree
        @test tree.root == 1
        @test length(tree.parent) == lat.N
        @test tree.parent[1] == -1  # Root has no parent
        
        # Test connectivity
        @test is_connected(tree)
        
        # Spanning tree should have N-1 edges
        @test length(tree.tree_edges) == lat.N - 1
        
        # Test different root
        tree2 = wilson_algorithm(lat, 50, rng=rng)
        @test tree2.root == 50
        @test tree2.parent[50] == -1
        @test is_connected(tree2)
    end
    
    @testset "Tree Properties" begin
        lat = build_lattice(Square, 8, 8)
        tree = wilson_algorithm(lat, rng=rng)
        
        # Test depth
        depth1 = tree_depth(tree, 1)
        @test depth1 == 0  # Root has depth 0
        
        depth_last = tree_depth(tree, lat.N)
        @test depth_last >= 0
        
        # Test height
        h = tree_height(tree)
        @test h >= 0
        @test h < lat.N  # Height must be less than number of vertices
        
        # Height should be max depth
        max_depth = maximum(tree_depth(tree, v) for v in 1:lat.N)
        @test h == max_depth
    end
    
    @testset "Paths" begin
        lat = build_lattice(Triangular, 10, 10)
        tree = wilson_algorithm(lat, 1, rng=rng)
        
        # Path to root from root
        path1 = path_to_root(tree, 1)
        @test length(path1) == 1
        @test path1[1] == 1
        
        # Path from arbitrary vertex
        path50 = path_to_root(tree, 50)
        @test path50[1] == 50
        @test path50[end] == 1  # Ends at root
        @test length(path50) == tree_depth(tree, 50) + 1
        
        # Check path validity
        for i in 1:(length(path50)-1)
            @test tree.parent[path50[i]] == path50[i+1]
        end
    end
    
    @testset "Tree Distance" begin
        lat = build_lattice(Square, 10, 10)
        tree = wilson_algorithm(lat, rng=rng)
        
        # Distance from vertex to itself
        @test tree_distance(tree, 1, 1) == 0
        @test tree_distance(tree, 50, 50) == 0
        
        # Distance is symmetric
        d_12 = tree_distance(tree, 10, 20)
        d_21 = tree_distance(tree, 20, 10)
        @test d_12 == d_21
        
        # Distance through root
        d_1_to_v = tree_distance(tree, 1, 50)
        depth_v = tree_depth(tree, 50)
        @test d_1_to_v == depth_v
    end
    
    @testset "Subtree Size" begin
        lat = build_lattice(Square, 10, 10)
        tree = wilson_algorithm(lat, rng=rng)
        
        # Root subtree includes all vertices
        @test subtree_size(tree, 1) == lat.N
        
        # Leaf has subtree size 1
        leaves = leaf_vertices(tree)
        @test !isempty(leaves)
        for leaf in leaves
            @test subtree_size(tree, leaf) == 1
        end
        
        # Sum of all subtree sizes >= N (with overcounting)
        total = sum(subtree_size(tree, v) for v in 1:lat.N)
        @test total >= lat.N
    end
    
    @testset "Leaf Vertices" begin
        lat = build_lattice(Square, 10, 10)
        tree = wilson_algorithm(lat, rng=rng)
        
        leaves = leaf_vertices(tree)
        @test !isempty(leaves)
        @test all(1 .<= leaves .<= lat.N)
        
        # Check leaf property: no vertex has this as parent
        for leaf in leaves
            has_child = any(tree.parent[i] == leaf for i in 1:lat.N)
            @test !has_child
        end
    end
    
    @testset "Degree Sequence" begin
        lat = build_lattice(Triangular, 10, 10)
        tree = wilson_algorithm(lat, rng=rng)
        
        degrees = degree_sequence(tree)
        @test length(degrees) == lat.N
        @test all(degrees .>= 0)
        
        # Sum of degrees = 2 * number of edges
        @test sum(degrees) == 2 * length(tree.tree_edges)
        
        # Tree has leaves (degree 1)
        @test any(degrees .== 1)
    end
    
    @testset "Average Path Length" begin
        lat = build_lattice(Square, 10, 10)
        tree = wilson_algorithm(lat, rng=rng)
        
        avg_len = average_path_length(tree, n_samples=50, rng=rng)
        @test avg_len > 0.0
        @test isfinite(avg_len)
        @test avg_len < lat.N  # Average path can't exceed N
    end
    
    @testset "Maze Generation" begin
        lat = build_lattice(Square, 15, 15, boundary=OBC())
        maze = generate_maze(lat, rng=rng)
        
        @test maze isa SpanningTree
        @test is_connected(maze)
        @test maze.root == 1  # Corner
    end
    
    @testset "Critical Path" begin
        lat = build_lattice(Square, 10, 10)
        tree = wilson_algorithm(lat, rng=rng)
        
        path, length = critical_path(tree)
        @test length >= 0
        @test length == Base.length(path) - 1
        
        # Path should be non-empty for non-trivial tree
        if lat.N > 1
            @test !isempty(path)
        end
    end
    
    @testset "Different Lattice Types" begin
        for LatticeType in [Square, Triangular, Honeycomb, Kagome]
            lat = build_lattice(LatticeType, 8, 8)
            tree = wilson_algorithm(lat, rng=rng)
            
            @test is_connected(tree)
            @test length(tree.tree_edges) == lat.N - 1
            @test tree_height(tree) >= 0
        end
    end
end
