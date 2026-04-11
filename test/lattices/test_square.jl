@testset "Square lattice" begin
    @testset "unit cell" begin
        uc = get_unit_cell(Square)
        a1, a2 = uc.basis
        @test dot(a1, a2) ≈ 0.0 atol = 1e-10
        @test norm(a1) ≈ 1.0
        @test norm(a2) ≈ 1.0
        @test length(uc.sublattice_positions) == 1
        @test length(uc.connections) == 2   # +x, +y
    end

    @testset "PBC connectivity" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Square, Lx, Ly)
        @test num_sublattices(lat) == 1
        # Every site has 4 neighbours
        @test all(length(neighbors(lat, i)) == 4 for i in 1:num_sites(lat))
        # Total bonds = Lx * Ly * (connections per cell)
        @test length(bonds(lat)) == Lx * Ly * 2
    end

    @testset "OBC corner / edge / bulk degree" begin
        lat = build_lattice(Square, 4, 4; boundary=OpenAxis())
        corner = site_index(RowMajor(), (4, 4), 1, LatticeCoord((1, 1)))
        edge = site_index(RowMajor(), (4, 4), 1, LatticeCoord((2, 1)))
        bulk = site_index(RowMajor(), (4, 4), 1, LatticeCoord((2, 2)))
        @test length(neighbors(lat, corner)) == 2
        @test length(neighbors(lat, edge)) == 3
        @test length(neighbors(lat, bulk)) == 4
    end

    @testset "positions match the expected row-major cell layout" begin
        Lx, Ly = 4, 3
        lat = build_lattice(Square, Lx, Ly)
        for x in 1:Lx, y in 1:Ly
            idx = site_index(RowMajor(), (Lx, Ly), 1, LatticeCoord((x, y)))
            @test position(lat, idx) == SVector(Float64(x - 1), Float64(y - 1))
        end
    end
end
