@testset "Shastry–Sutherland lattice" begin
    @testset "unit cell" begin
        uc = get_unit_cell(ShastrySutherland)
        @test length(uc.sublattice_positions) == 4
        # Eight square-lattice bonds (four intra-cell + four inter-cell)
        # plus two diagonal dimer bonds.
        @test length(uc.connections) == 10
    end

    @testset "bond count under PBC" begin
        Lx, Ly = 4, 4
        lat = build_lattice(ShastrySutherland, Lx, Ly)
        @test num_sublattices(lat) == 4
        @test num_sites(lat) == Lx * Ly * 4
        @test count(_ -> true, bonds(lat)) ==
            Lx * Ly * length(get_unit_cell(ShastrySutherland).connections)
    end

    @testset "every site has 5 neighbours (4 NN + 1 dimer partner)" begin
        lat = build_lattice(ShastrySutherland, 4, 4)
        @test all(length(neighbors(lat, i)) == 5 for i in 1:num_sites(lat))
    end
end
