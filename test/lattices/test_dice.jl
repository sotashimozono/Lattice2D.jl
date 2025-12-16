@testset "Dice Lattice Tests" begin
    @testset "Unit Cell Definition" begin
        uc = get_unit_cell(Dice)
        # Basis Vectors are separated by 60 degrees and unit length (triangular)
        a1, a2 = uc.basis
        @test dot(a1, a2) ≈ cos(π/3) atol=1e-10
        @test norm(a1) ≈ 1.0 atol=1e-10
        @test norm(a2) ≈ 1.0 atol=1e-10
        # 3 sublattices: Hub(1), RimA(2), RimB(3)
        @test length(uc.sublattice_positions) == 3

        d_1 = uc.sublattice_positions[1]
        d_2 = uc.sublattice_positions[2]
        d_3 = uc.sublattice_positions[3]
        @test d_1 == [0.0, 0.0]
        # Check relative positions (Should be center of triangles)
        @test d_2 ≈ (a1 .+ a2) ./ 3
        @test d_3 ≈ (a1 .+ a2) .* (2 / 3)
    end

    @testset "Geometry & Reciprocal Vectors" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Dice, Lx, Ly)

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
            lat_pbc = build_lattice(Dice, Lx, Ly; boundary=PBC())

            # total number of sites will be 3 * Lx * Ly
            @test lat_pbc.N == 3 * Lx * Ly

            # Dice lattice: Hub has coordination 6, Rims have coordination 3
            degrees = length.(lat_pbc.nearest_neighbors)
            sub_ids = lat_pbc.sublattice_ids

            # Hub (sublattice 1) should have degree 6
            @test all(degrees[sub_ids .== 1] .== 6)
            # Rim A and B (sublattices 2, 3) should have degree 3
            @test all(degrees[sub_ids .== 2] .== 3)
            @test all(degrees[sub_ids .== 3] .== 3)

            # --- Connectivity Check (Corrected for proper Dice Geometry) ---
            # Get IDs for cell (1,1)
            id_Hub_11 = lat_pbc.site_map[1, 1]
            id_RimA_11 = id_Hub_11 + 1 # Sublattice 2
            id_RimB_11 = id_Hub_11 + 2 # Sublattice 3

            # 1. RimA (1,1) IS connected to Hub (1,1) [Intra-cell bond]
            @test id_RimA_11 in lat_pbc.nearest_neighbors[id_Hub_11]
            @test id_Hub_11 in lat_pbc.nearest_neighbors[id_RimA_11]

            # 2. RimB (1,1) is NOT connected to Hub (1,1) in correct geometry
            #    (It connects to Hubs in neighbor cells)
            @test !(id_RimB_11 in lat_pbc.nearest_neighbors[id_Hub_11])

            # 3. Check correct neighbor for RimB (1,1)
            #    Based on the fix: RimB(x,y) connects to Hub(x+1, y+1) due to PBC
            #    So RimB(1,1) should connect to Hub(2,2)
            id_Hub_22 = lat_pbc.site_map[2, 2]
            @test id_Hub_22 in lat_pbc.nearest_neighbors[id_RimB_11]

            # Bipartite check - Dice is bipartite (Hub vs Rims)
            @test lat_pbc.is_bipartite == true
        end

        @testset "open boundary condition" begin
            lat_obc = build_lattice(Dice, Lx, Ly; boundary=OBC())

            degrees = length.(lat_obc.nearest_neighbors)
            @test maximum(degrees) <= 6
            @test minimum(degrees) >= 0

            bulk_hub = lat_obc.site_map[2, 2]
            @test length(lat_obc.nearest_neighbors[bulk_hub]) <= 6

            bulk_rim = lat_obc.site_map[2, 2] + 1
            @test length(lat_obc.nearest_neighbors[bulk_rim]) > 0
        end
    end

    @testset "Index Consistency" begin
        Lx, Ly = 3, 3
        n_sub = 3 # Dice has 3 sublattices
        lat = build_lattice(Dice, Lx, Ly)
        uc = get_unit_cell(Dice)

        for x in 1:Lx, y in 1:Ly
            base_idx = lat.site_map[x, y]

            # Base index calc check (assuming RowMajor-like or standard construction)
            expected_base = ((x - 1) + (y - 1) * Lx) * n_sub + 1
            @test base_idx == expected_base

            # Sublattice position check
            cell_origin = (x-1)*lat.basis_vectors[1] + (y-1)*lat.basis_vectors[2]

            # Check Hub, RimA, RimB positions
            for s in 1:n_sub
                pos = lat.positions[base_idx + s - 1]
                expected_pos = cell_origin + uc.sublattice_positions[s]
                @test pos ≈ expected_pos atol=1e-10
            end
        end
    end
end
