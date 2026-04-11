@testset "Union Jack lattice" begin
    @testset "unit cell" begin
        uc = get_unit_cell(UnionJack)
        a1, a2 = uc.basis
        @test a1 ≈ [1.0, 0.0]
        @test a2 ≈ [0.0, 1.0]
        @test length(uc.sublattice_positions) == 2   # corner + body centre
        @test length(uc.connections) == 6
    end

    @testset "corner / centre coordination under PBC" begin
        Lx, Ly = 4, 4
        lat = build_lattice(UnionJack, Lx, Ly)
        @test num_sublattices(lat) == 2

        corner_sites = [i for i in 1:num_sites(lat) if sublattice(lat, i) == 1]
        body_sites = [i for i in 1:num_sites(lat) if sublattice(lat, i) == 2]

        # Corner (A) sites: 4 square-NN + 4 body-centre neighbours = 8
        for i in corner_sites
            @test length(neighbors(lat, i)) == 8
        end

        # Body-centre (B) sites: 4 corner neighbours
        for i in body_sites
            @test length(neighbors(lat, i)) == 4
        end
    end
end
