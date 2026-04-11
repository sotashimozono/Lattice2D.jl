@testset "Dice (T3) lattice" begin
    @testset "unit cell" begin
        uc = get_unit_cell(Dice)
        a1, a2 = uc.basis
        @test a1 ≈ [1.0, 0.0]
        @test a2 ≈ [0.5, sqrt(3) / 2]
        @test length(uc.sublattice_positions) == 3   # hub + 2 rims
        @test length(uc.connections) == 6
    end

    @testset "hub / rim coordination numbers under PBC" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Dice, Lx, Ly)
        @test num_sublattices(lat) == 3
        @test num_sites(lat) == Lx * Ly * 3

        hub_sites = [i for i in 1:num_sites(lat) if sublattice(lat, i) == 1]
        rim_sites = [i for i in 1:num_sites(lat) if sublattice(lat, i) in (2, 3)]

        # Hub site has 6 nearest neighbours (3 from each rim sublattice)
        for i in hub_sites
            @test length(neighbors(lat, i)) == 6
        end
        # Rim sites have 3 nearest neighbours (all hubs)
        for i in rim_sites
            @test length(neighbors(lat, i)) == 3
        end
    end

    @testset "bipartite: hub vs rim" begin
        @test is_bipartite(build_lattice(Dice, 4, 4)) == true
    end
end
