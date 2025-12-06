@testset "Honeycomb Lattice Tests" begin
    @testset "Unit Cell Definition" begin
        uc = get_unit_cell(Honeycomb)
        # Basis Vectors are separated by 60 degrees and unit length
        a1, a2 = uc.basis
        @test dot(a1, a2) ≈ 1.5 atol=1e-10
        @test norm(a1) ≈ sqrt(3) atol=1e-10
        @test norm(a2) ≈ sqrt(3) atol=1e-10
        # Single sublattice at origin
        @test length(uc.sublattice_positions) == 2
        @test uc.sublattice_positions[1] == [0.0, 0.0]
        @test uc.sublattice_positions[2] == [0.0, 1.0]
    end

    @testset "Geometry & Reciprocal Vectors" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Honeycomb, Lx, Ly)

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
            lat_pbc = build_lattice(Honeycomb, Lx, Ly; boundary=PBC())

            # total number of sites will be 2 * Lx * Ly
            @test lat_pbc.N == 2 * Lx * Ly
            # each site should have 3 neighbors in PBC
            degrees = length.(lat_pbc.nearest_neighbors)
            @test all(d -> d == 3, degrees)
            
            
            id_1_1_A = lat_pbc.site_map[1, 1]     # (1,1) A
            id_1_1_B = lat_pbc.site_map[1, 1] + 1 # (1,1) B
            
            # Check if the neighbors contain the corresponding IDs
            @test id_1_1_B in lat_pbc.nearest_neighbors[id_1_1_A]
            @test id_1_1_A in lat_pbc.nearest_neighbors[id_1_1_B]

            # Bipartite check
            @test lat_pbc.is_bipartite == true
        end
        @testset "open boundary condition" begin
            lat_obc = build_lattice(Honeycomb, Lx, Ly; boundary=OBC())
            
            degrees = length.(lat_obc.nearest_neighbors)
            @test maximum(degrees) == 3
            @test minimum(degrees) >= 1

            # Bulk: coordination number 6
            bulk_id = lat_obc.site_map[1, 1]
            @test length(lat_obc.nearest_neighbors[bulk_id]) <= 3
        end
    end

    @testset "Index Consistency" begin
        Lx, Ly = 3, 3
        n_sub = 2 # Honeycomb has 2 sublattices
        lat = build_lattice(Honeycomb, Lx, Ly)
        
        for x in 1:Lx, y in 1:Ly
            base_idx = lat.site_map[x, y]
            
            # 手動計算との一致: ((x-1) + (y-1)*Lx) * n_sub + 1
            expected_base = ((x - 1) + (y - 1) * Lx) * n_sub + 1
            @test base_idx == expected_base
            
            # 各サブ格子の座標チェック
            for s in 1:n_sub
                idx = base_idx + (s - 1)
                
                # lattice position
                pos = lat.positions[idx]
                
                # expected position: cell origin + sublattice displacement
                cell_origin = (x-1)*lat.basis_vectors[1] + (y-1)*lat.basis_vectors[2]
                expected_pos = cell_origin + lat.unit_cell.sublattice_positions[s]
                
                @test pos ≈ expected_pos
            end
        end
    end
end