@testset "Lieb lattice" begin
    @testset "unit cell" begin
        uc = get_unit_cell(Lieb)
        a1, a2 = uc.basis
        @test a1 ≈ [2.0, 0.0]
        @test a2 ≈ [0.0, 2.0]
        @test length(uc.sublattice_positions) == 3   # A corner, B x-edge, C y-edge
        @test length(uc.connections) == 4
    end

    @testset "PBC connectivity: corner has 4 neighbours, edges have 2" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Lieb, Lx, Ly)
        @test num_sublattices(lat) == 3
        @test num_sites(lat) == Lx * Ly * 3

        # A (corner) sites connect to 4 edge sites (2 B, 2 C)
        a_sites = [i for i in 1:num_sites(lat) if sublattice(lat, i) == 1]
        for i in a_sites
            @test length(neighbors(lat, i)) == 4
        end

        # B, C (edge) sites each connect to 2 A sites
        edge_sites = [i for i in 1:num_sites(lat) if sublattice(lat, i) in (2, 3)]
        for i in edge_sites
            @test length(neighbors(lat, i)) == 2
        end

        @test count(_ -> true, bonds(lat)) == Lx * Ly * 4
    end

    @testset "bipartite (A corners vs B/C edges)" begin
        @test is_bipartite(build_lattice(Lieb, 4, 4)) == true
    end
end
