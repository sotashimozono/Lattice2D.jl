@testset "Kagome lattice" begin
    @testset "unit cell" begin
        uc = get_unit_cell(Kagome)
        a1, a2 = uc.basis
        @test a1 ≈ [1.0, 0.0]
        @test a2 ≈ [0.5, sqrt(3) / 2]
        @test length(uc.sublattice_positions) == 3   # A, B, C
        @test length(uc.connections) == 6            # 3 intra + 3 inter
    end

    @testset "PBC connectivity: 4 neighbours per site" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Kagome, Lx, Ly)
        @test num_sublattices(lat) == 3
        @test num_sites(lat) == Lx * Ly * 3
        # Every site of Kagome under PBC has 4 nearest neighbours
        @test all(length(neighbors(lat, i)) == 4 for i in 1:num_sites(lat))
        @test count(_ -> true, bonds(lat)) == Lx * Ly * 6
    end

    @testset "triangles prevent bipartiteness" begin
        @test is_bipartite(build_lattice(Kagome, 4, 4)) == false
    end
end
