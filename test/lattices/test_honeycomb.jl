@testset "Honeycomb lattice" begin
    @testset "unit cell" begin
        uc = get_unit_cell(Honeycomb)
        a1, a2 = uc.basis
        @test a1 ≈ [sqrt(3), 0.0]
        @test a2 ≈ [sqrt(3) / 2, 1.5]
        @test length(uc.sublattice_positions) == 2   # A, B
        @test length(uc.connections) == 3
    end

    @testset "PBC connectivity: 3 neighbours per site" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Honeycomb, Lx, Ly)
        @test num_sublattices(lat) == 2
        @test num_sites(lat) == Lx * Ly * 2
        @test all(length(neighbors(lat, i)) == 3 for i in 1:num_sites(lat))
        @test length(bonds(lat)) == Lx * Ly * 3
    end

    @testset "bipartite A / B split" begin
        lat = build_lattice(Honeycomb, 4, 4)
        @test is_bipartite(lat) == true
        # A and B sublattices partition the sites cleanly
        a_sites = [i for i in 1:num_sites(lat) if sublattice(lat, i) == 1]
        b_sites = [i for i in 1:num_sites(lat) if sublattice(lat, i) == 2]
        @test length(a_sites) == length(b_sites) == num_sites(lat) ÷ 2
        # Every bond crosses the A / B boundary
        for b in bonds(lat)
            @test sublattice(lat, b.i) != sublattice(lat, b.j)
        end
    end
end
