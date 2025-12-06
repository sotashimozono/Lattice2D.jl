@testset "Kagome Lattice Tests" begin
    @testset "Unit Cell Definition" begin
        uc = get_unit_cell(Kagome)
        # Basis Vectors are separated by 60 degrees and unit length
        a1, a2 = uc.basis
        @test dot(a1, a2) ≈ cos(π/3) atol=1e-10
        @test norm(a1) ≈ sqrt(1) atol=1e-10
        @test norm(a2) ≈ sqrt(1) atol=1e-10
        # Single sublattice at origin
        @test length(uc.sublattice_positions) == 3

        d_A = uc.sublattice_positions[1]
        d_B = uc.sublattice_positions[2]
        d_C = uc.sublattice_positions[3]
        @test d_A == [0.0, 0.0]
        @test d_B == 0.5 * a1
        @test d_C == 0.5 * a2
    end

    @testset "Geometry & Reciprocal Vectors" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Kagome, Lx, Ly)

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
            lat_pbc = build_lattice(Kagome, Lx, Ly; boundary=PBC())

            # total number of sites will be 3 * Lx * Ly
            @test lat_pbc.N == 3 * Lx * Ly
            # each site should have 4 neighbors in PBC
            degrees = length.(lat_pbc.nearest_neighbors)
            @test all(d -> d == 4, degrees)
            
            id_A = lat_pbc.site_map[1, 1]
            id_B = id_A + 1
            id_C = id_A + 2
            # Check if the neighbors contain the corresponding IDs
            @test id_B in lat_pbc.nearest_neighbors[id_A]
            @test id_A in lat_pbc.nearest_neighbors[id_B]

            # Bipartite check
            @test lat_pbc.is_bipartite == false
        end
        @testset "open boundary condition" begin
            lat_obc = build_lattice(Kagome, Lx, Ly; boundary=OBC())
            
            degrees = length.(lat_obc.nearest_neighbors)
            @test maximum(degrees) == 4
            @test minimum(degrees) >= 2

            # Bulk: coordination number 2
            bulk_id = lat_obc.site_map[1, 1]
            @test length(lat_obc.nearest_neighbors[bulk_id]) <= 2
        end
    end
@testset "Index Consistency" begin
        Lx, Ly = 3, 3
        n_sub = 3 # Kagome is 3
        lat = build_lattice(Kagome, Lx, Ly)
        
        for x in 1:Lx, y in 1:Ly
            base_idx = lat.site_map[x, y]
            
            # Base index calc check
            expected_base = ((x - 1) + (y - 1) * Lx) * n_sub + 1
            @test base_idx == expected_base
            
            # Sublattice position check
            cell_origin = (x-1)*lat.basis_vectors[1] + (y-1)*lat.basis_vectors[2]
            
            # Check A, B, C positions
            for s in 1:n_sub
                pos = lat.positions[base_idx + s - 1]
                expected_pos = cell_origin + lat.unit_cell.sublattice_positions[s]
                @test pos ≈ expected_pos
            end
        end
    end
end