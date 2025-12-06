@testset "Shastry-Sutherland Lattice Tests" begin
    @testset "Unit Cell Definition" begin
        uc = get_unit_cell(ShastrySutherland)
        
        # Based on 2x2 square unit cell
        a1, a2 = uc.basis
        @test norm(a1) ≈ 2.0 atol=1e-10
        @test norm(a2) ≈ 2.0 atol=1e-10
        @test dot(a1, a2) ≈ 0.0 atol=1e-10 # Orthogonal

        @test length(uc.sublattice_positions) == 4
        
        # Check positions relative to basis (assuming standard definition)
        # d1=[0,0], d2=[1,0], d3=[0,1], d4=[1,1]
        @test uc.sublattice_positions[1] == [0.0, 0.0]
        @test uc.sublattice_positions[4] == [1.0, 1.0]
    end

    @testset "Geometry & Reciprocal Vectors" begin
        Lx, Ly = 4, 4
        lat = build_lattice(ShastrySutherland, Lx, Ly)
        
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
            lat_pbc = build_lattice(ShastrySutherland, Lx, Ly; boundary=PBC())
            
            @test lat_pbc.N == 4 * Lx * Ly

            # Square Lattice (4 neighbors) + 1 Diagonal Dimer = 5
            degrees = length.(lat_pbc.nearest_neighbors)
            @test all(d -> d == 5, degrees)
            
            # Diagonal bonds create triangles (e.g., 1-2-4), breaking bipartiteness
            @test lat_pbc.is_bipartite == false

            # Check if Site 1 connects to Site 4 (Intra-cell diagonal)
            # This depends on site_map layout, but typically site 1 and 4 in the same cell form a dimer
            id_1 = lat_pbc.site_map[1, 1]
            id_4 = id_1 + 3 # 4th sublattice
            @test id_4 in lat_pbc.nearest_neighbors[id_1]
        end

        @testset "open boundary condition" begin
            lat_obc = build_lattice(ShastrySutherland, Lx, Ly; boundary=OBC())
            
            degrees = length.(lat_obc.nearest_neighbors)
            @test maximum(degrees) == 5
            
            corner_id = lat_obc.site_map[1, 1]
            @test length(lat_obc.nearest_neighbors[corner_id]) == 3
        end
    end

    @testset "Index Consistency" begin
        Lx, Ly = 3, 3
        n_sub = 4 # Shastry-Sutherland is 4
        lat = build_lattice(ShastrySutherland, Lx, Ly)
        
        for x in 1:Lx, y in 1:Ly
            base_idx = lat.site_map[x, y]
            
            expected_base = ((x - 1) + (y - 1) * Lx) * n_sub + 1
            @test base_idx == expected_base
            
            cell_origin = (x-1)*lat.basis_vectors[1] + (y-1)*lat.basis_vectors[2]
            for s in 1:n_sub
                pos = lat.positions[base_idx + s - 1]
                expected_pos = cell_origin + lat.unit_cell.sublattice_positions[s]
                @test pos ≈ expected_pos
            end
        end
    end
end