@testset "Lieb Lattice Tests" begin
    @testset "Unit Cell Definition" begin
        uc = get_unit_cell(Lieb)
        
        # Lieb has doubled unit cell size: a1=[2.0, 0.0], a2=[0.0, 2.0]
        a1, a2 = uc.basis
        @test norm(a1) ≈ 2.0 atol=1e-10
        @test norm(a2) ≈ 2.0 atol=1e-10
        @test dot(a1, a2) ≈ 0.0 atol=1e-10

        # 2. Sublattices check (3 sites: Corner, Right, Up)
        @test length(uc.sublattice_positions) == 3
        
        d_A = uc.sublattice_positions[1] # [0,0]
        d_B = uc.sublattice_positions[2] # [1,0]
        d_C = uc.sublattice_positions[3] # [0,1]

        @test d_A == [0.0, 0.0]
        @test d_B ≈ 0.5 * a1
        @test d_C ≈ 0.5 * a2
    end

    @testset "Geometry & Reciprocal Vectors" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Lieb, Lx, Ly)
        
        a = lat.basis_vectors
        b = lat.reciprocal_vectors

        # a_i ⋅ b_j = 2π δ_ij
        @test dot(a[1], b[1]) ≈ 2π atol=1e-10
        @test dot(a[2], b[2]) ≈ 2π atol=1e-10
        @test abs(dot(a[1], b[2])) < 1e-10
        @test abs(dot(a[2], b[1])) < 1e-10
    end

    @testset "Topology & Connectivity" begin
        Lx, Ly = 4, 4
        
        @testset "periodic boundary condition" begin
            lat_pbc = build_lattice(Lieb, Lx, Ly; boundary=PBC())
            
            @test lat_pbc.N == 3 * Lx * Ly

            # Sublattice 1 (A): Degree 4
            # Sublattice 2,3 (B, C): Degree 2
            degrees = length.(lat_pbc.nearest_neighbors)
            sub_ids = lat_pbc.sublattice_ids
            
            @test all(degrees[sub_ids .== 1] .== 4)
            @test all(degrees[sub_ids .== 2] .== 2)
            @test all(degrees[sub_ids .== 3] .== 2)
            
            # A connects only to B/C. B/C connects only to A.
            # No B-C connection.
            @test lat_pbc.is_bipartite == true

            # Check A(1,1) connections
            id_A = lat_pbc.site_map[1, 1]
            id_B = id_A + 1
            id_C = id_A + 2
            
            # A should connect to B and C in the same cell
            @test id_B in lat_pbc.nearest_neighbors[id_A]
            @test id_C in lat_pbc.nearest_neighbors[id_A]
        end

        @testset "open boundary condition" begin
            lat_obc = build_lattice(Lieb, Lx, Ly; boundary=OBC())
            
            degrees = length.(lat_obc.nearest_neighbors)
            @test maximum(degrees) == 4
            
            # Corner A-site (1,1): Connects to B(1,1) and C(1,1) -> Degree 2
            # Corner B/C sites might have degree 1 (e.g., B at right edge)
            @test minimum(degrees) >= 1
            
            # Specific Corner Check (1,1)
            # A(1,1) -> B(1,1), C(1,1). (Left and Down are open)
            corner_A = lat_obc.site_map[1, 1]
            @test length(lat_obc.nearest_neighbors[corner_A]) == 2
        end
    end

    @testset "Index Consistency" begin
        Lx, Ly = 3, 3
        n_sub = 3
        lat = build_lattice(Lieb, Lx, Ly)
        
        for x in 1:Lx, y in 1:Ly
            base_idx = lat.site_map[x, y]
            
            # Base index check
            expected_base = ((x - 1) + (y - 1) * Lx) * n_sub + 1
            @test base_idx == expected_base
            
            # Position check
            cell_origin = (x-1)*lat.basis_vectors[1] + (y-1)*lat.basis_vectors[2]
            for s in 1:n_sub
                pos = lat.positions[base_idx + s - 1]
                expected_pos = cell_origin + lat.unit_cell.sublattice_positions[s]
                @test pos ≈ expected_pos
            end
        end
    end
end