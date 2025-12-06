@testset "Triangular Lattice Tests" begin
    @testset "Unit Cell Definition" begin
        uc = get_unit_cell(Triangular)
        # Basis Vectors are separated by 60 degrees and unit length
        a1, a2 = uc.basis
        @test dot(a1, a2) ≈ cos(pi/3) atol=1e-10
        @test norm(a1) ≈ 1.0 atol=1e-10
        @test norm(a2) ≈ 1.0 atol=1e-10
        # Single sublattice at origin
        @test length(uc.sublattice_positions) == 1
        @test uc.sublattice_positions[1] == [0.0, 0.0]
    end

    @testset "Geometry & Reciprocal Vectors" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Triangular, Lx, Ly)

        # a is the basis vectors, b is the reciprocal vectors
        # Verify a_i · b_j = 2π δ_ij
        a = lat.basis_vectors
        b = lat.reciprocal_vectors

        @test dot(a[1], b[1]) ≈ 2π atol=1e-10
        @test dot(a[2], b[2]) ≈ 2π atol=1e-10
        @test abs(dot(a[1], b[2])) < 1e-10
        @test abs(dot(a[2], b[1])) < 1e-10
    end

    @testset "Topology & Connectivity" begin
        Lx, Ly = 4, 4
        @testset "periodic boundary condition" begin
            lat_pbc = build_lattice(Triangular, Lx, Ly; boundary=PBC())
            
            # each site should have 6 neighbors in PBC
            degrees = length.(lat_pbc.nearest_neighbors)
            @test all(d -> d == 6, degrees)
            
            # Check specific connections: is the left neighbor of (1,1) connected to (Lx,1)?
            # Use site_map[x, y] to get IDs for verification
            id_1_1 = lat_pbc.site_map[1, 1]
            id_L_1 = lat_pbc.site_map[Lx, 1] # left edge
            id_1_L = lat_pbc.site_map[1, Ly] # bottom edge
            
            # Check if the neighbors contain the corresponding IDs
            @test id_L_1 in lat_pbc.nearest_neighbors[id_1_1]
            @test id_1_L in lat_pbc.nearest_neighbors[id_1_1]

            # Bipartite check
            @test lat_pbc.is_bipartite == false
        end
        @testset "open boundary condition" begin
            lat_obc = build_lattice(Triangular, Lx, Ly; boundary=OBC())
            
            # Corner: coordination number 2
            corner_id = lat_obc.site_map[1, 1]
            @test length(lat_obc.nearest_neighbors[corner_id]) == 2
            
            # Edge: coordination number 4 (1 < x < Lx, y=1)
            edge_id = lat_obc.site_map[2, 1]
            @test length(lat_obc.nearest_neighbors[edge_id]) == 4
            
            # Bulk: coordination number 6
            bulk_id = lat_obc.site_map[2, 2]
            @test length(lat_obc.nearest_neighbors[bulk_id]) == 6
        end
    end

    @testset "Index Consistency" begin
        Lx, Ly = 4, 3
        lat = build_lattice(Triangular, Lx, Ly)
        
        for x in 1:Lx, y in 1:Ly
            expected_idx = (x - 1) + (y - 1) * Lx + 1
            @test lat.site_map[x, y] == expected_idx
            
            pos = lat.positions[expected_idx]
            expected_pos = (x-1)*lat.basis_vectors[1] + (y-1)*lat.basis_vectors[2]
            @test pos ≈ expected_pos
        end
    end
end