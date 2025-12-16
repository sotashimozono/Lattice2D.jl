@testset "UnionJack Lattice Tests" begin
    @testset "Unit Cell Definition" begin
        uc = get_unit_cell(UnionJack)
        # Basis Vectors are orthogonal and unit length
        a1, a2 = uc.basis
        @test dot(a1, a2) ≈ 0.0 atol=1e-10
        @test norm(a1) ≈ 1.0 atol=1e-10
        @test norm(a2) ≈ 1.0 atol=1e-10
        # 2 sublattices: Corner(A) and Center(B)
        @test length(uc.sublattice_positions) == 2

        d_A = uc.sublattice_positions[1]
        d_B = uc.sublattice_positions[2]
        @test d_A == [0.0, 0.0]
        @test d_B ≈ [0.5, 0.5]
    end

    @testset "Geometry & Reciprocal Vectors" begin
        Lx, Ly = 4, 4
        lat = build_lattice(UnionJack, Lx, Ly)

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
            lat_pbc = build_lattice(UnionJack, Lx, Ly; boundary=PBC())

            # total number of sites will be 2 * Lx * Ly
            @test lat_pbc.N == 2 * Lx * Ly

            # UnionJack: Corner has coordination 8 (4 square + 4 diagonal), Center has coordination 4
            degrees = length.(lat_pbc.nearest_neighbors)
            sub_ids = lat_pbc.sublattice_ids

            # Corner (sublattice 1) should have degree 8
            @test all(degrees[sub_ids .== 1] .== 8)
            # Center (sublattice 2) should have degree 4
            @test all(degrees[sub_ids .== 2] .== 4)

            id_A = lat_pbc.site_map[1, 1]
            id_B = id_A + 1
            # Check if the neighbors contain the corresponding IDs
            @test id_B in lat_pbc.nearest_neighbors[id_A]
            @test id_A in lat_pbc.nearest_neighbors[id_B]

            # Bipartite check - UnionJack is NOT bipartite (Corner sites connect to each other)
            @test lat_pbc.is_bipartite == false
        end
        @testset "open boundary condition" begin
            lat_obc = build_lattice(UnionJack, Lx, Ly; boundary=OBC())

            degrees = length.(lat_obc.nearest_neighbors)
            @test maximum(degrees) <= 8
            @test minimum(degrees) >= 1

            # Corner site at (1,1): has Right, Up, and Center connections
            corner_A = lat_obc.site_map[1, 1]
            # It connects to B(1,1), A(2,1), A(1,2) -> at least 3
            @test length(lat_obc.nearest_neighbors[corner_A]) >= 3

            # Bulk Corner site should have coordination close to 8
            bulk_A = lat_obc.site_map[2, 2]
            @test length(lat_obc.nearest_neighbors[bulk_A]) <= 8
        end
    end
    @testset "Index Consistency" begin
        Lx, Ly = 3, 3
        n_sub = 2 # UnionJack has 2 sublattices
        lat = build_lattice(UnionJack, Lx, Ly)
        uc = get_unit_cell(UnionJack)

        for x in 1:Lx, y in 1:Ly
            base_idx = lat.site_map[x, y]

            # Base index calc check
            expected_base = ((x - 1) + (y - 1) * Lx) * n_sub + 1
            @test base_idx == expected_base

            # Sublattice position check
            cell_origin = (x-1)*lat.basis_vectors[1] + (y-1)*lat.basis_vectors[2]

            # Check Corner and Center positions
            for s in 1:n_sub
                pos = lat.positions[base_idx + s - 1]
                expected_pos = cell_origin + uc.sublattice_positions[s]
                @test pos ≈ expected_pos
            end
        end
    end
end
