@testset "Triangular lattice" begin
    @testset "unit cell" begin
        uc = get_unit_cell(Triangular)
        a1, a2 = uc.basis
        @test a1 ≈ [1.0, 0.0]
        @test a2 ≈ [0.5, sqrt(3) / 2]
        @test length(uc.sublattice_positions) == 1
        @test length(uc.connections) == 3   # a1, a2, a2-a1
    end

    @testset "PBC connectivity: 6 neighbours per site" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Triangular, Lx, Ly)
        @test num_sublattices(lat) == 1
        @test all(length(neighbors(lat, i)) == 6 for i in 1:num_sites(lat))
        @test length(bonds(lat)) == Lx * Ly * 3
    end

    @testset "never bipartite under PBC" begin
        @test is_bipartite(build_lattice(Triangular, 4, 4)) == false
        @test is_bipartite(build_lattice(Triangular, 6, 6)) == false
    end
end
