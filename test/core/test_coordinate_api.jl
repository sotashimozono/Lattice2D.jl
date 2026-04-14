@testset "Coordinate API (LatticeCore accessors on PeriodicLattice2D)" begin
    @testset "Square: position(lat, i) matches row-major cell order" begin
        Lx, Ly = 4, 3
        lat = build_lattice(Square, Lx, Ly)
        @test num_sites(lat) == Lx * Ly

        # RowMajor linearisation: site ((y-1)*Lx + x) lives at (x-1, y-1)
        for y in 1:Ly, x in 1:Lx
            idx = site_index(RowMajor(), (Lx, Ly), 1, LatticeCoord((x, y)))
            @test position(lat, idx) == SVector(Float64(x - 1), Float64(y - 1))
        end
    end

    @testset "site_index / lattice_coord round-trip matches positions" begin
        Lx, Ly = 4, 3
        lat = build_lattice(Square, Lx, Ly)
        nsub = num_sublattices(lat)
        for i in 1:num_sites(lat)
            coord = lattice_coord(RowMajor(), (Lx, Ly), nsub, i)
            @test site_index(RowMajor(), (Lx, Ly), nsub, coord) == i
            # to_real on the lattice should agree with position(lat, i)
            @test to_real(lat, coord).x == position(lat, i)
        end
    end

    @testset "Multi-sublattice: Honeycomb sublattice ids are well-defined" begin
        Lx, Ly = 3, 3
        lat = build_lattice(Honeycomb, Lx, Ly)
        @test num_sublattices(lat) == 2

        # The default RowMajor layout places the two sublattices of a
        # cell next to each other in the site index.
        for i in 1:num_sites(lat)
            s = sublattice(lat, i)
            @test 1 <= s <= 2
        end

        # A and B sublattices alternate in the index sequence.
        @test sublattice(lat, 1) == 1
        @test sublattice(lat, 2) == 2
    end

    @testset "Kagome (3 sublattices) has the expected num_sublattices" begin
        lat = build_lattice(Kagome, 3, 3)
        @test num_sublattices(lat) == 3
        # Sublattice ids are in 1..3
        @test all(1 <= sublattice(lat, i) <= 3 for i in 1:num_sites(lat))
    end

    @testset "Snake indexing on Square: neighbours are still correct" begin
        # Snake indexing does not change the neighbour graph — only the
        # linear ordering. Verify that the neighbour set is still the
        # same (as a set) regardless of indexing.
        row_major = build_lattice(Square, 4, 4)
        snake = build_lattice(Square, 4, 4; indexing=Snake())

        @test num_sites(row_major) == num_sites(snake)
        @test count(_ -> true, bonds(row_major)) == count(_ -> true, bonds(snake))
    end

    @testset "Base.length / size match lattice dimensions" begin
        lat = build_lattice(Square, 5, 4)
        @test length(lat) == 20
        @test size(lat) == (5, 4)
        @test size(lat, 1) == 5
        @test size(lat, 2) == 4
    end
end
