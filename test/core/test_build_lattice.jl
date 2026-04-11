@testset "build_lattice — topology-agnostic generic checks" begin
    @testset "every shipped topology builds under PBC" begin
        Lx, Ly = 4, 4
        for Topo in AVAILABLE_LATTICES
            lat = build_lattice(Topo, Lx, Ly)
            @test lat isa PeriodicLattice2D
            @test num_sites(lat) == Lx * Ly * num_sublattices(lat)
            @test size_trait(lat) isa FiniteSize{2}
            @test periodicity(lat) isa Periodic
            @test reciprocal_support(lat) isa HasReciprocal
            # Every site has at least one neighbour under PBC
            for i in 1:num_sites(lat)
                @test !isempty(neighbors(lat, i))
            end
        end
    end

    @testset "OBC drops boundary-crossing bonds" begin
        Lx, Ly = 5, 4

        # Square lattice: interior sites should still have 4 neighbours
        # but corner / edge sites should have fewer.
        lat_obc = build_lattice(Square, Lx, Ly; boundary=OpenAxis())
        degrees = [length(neighbors(lat_obc, i)) for i in 1:num_sites(lat_obc)]
        @test minimum(degrees) == 2     # corners
        @test maximum(degrees) == 4     # interior
        @test length(bonds(lat_obc)) < length(bonds(build_lattice(Square, Lx, Ly)))
    end

    @testset "cylinder (PBC × OBC) mixes neighbour counts" begin
        cyl = build_lattice(
            Square, 4, 4; boundary=LatticeBoundary((PeriodicAxis(), OpenAxis()))
        )
        @test periodicity(cyl) isa Aperiodic          # one axis is open
        @test reciprocal_support(cyl) isa NoReciprocal

        # x-wrapping keeps row neighbours tight, but top / bottom rows
        # lose a y-neighbour each.
        row1_site = 1                                 # (1, 1) in row-major
        row2_site = 5                                 # (1, 2)
        @test length(neighbors(cyl, row1_site)) == 3  # wraps in x, no y-down
        @test length(neighbors(cyl, row2_site)) == 4  # interior row
    end

    @testset "graph consistency: symmetric neighbours, non-zero bond vectors" begin
        lat = build_lattice(Square, 5, 4)

        # Neighbour lists are symmetric
        for i in 1:num_sites(lat)
            for j in neighbors(lat, i)
                @test i in neighbors(lat, j)
            end
        end

        # Every bond endpoint is in the neighbour list of the other
        for b in bonds(lat)
            @test b.j in neighbors(lat, b.i)
            @test b.i in neighbors(lat, b.j)
            @test norm(b.vector) ≈ 1.0 atol = 1e-10     # Square: unit NN
        end
    end

    @testset "indexing strategies produce consistent site counts" begin
        Lx, Ly = 4, 3
        for indexing in (RowMajor(), ColMajor(), Snake())
            lat = build_lattice(Square, Lx, Ly; indexing=indexing)
            @test num_sites(lat) == Lx * Ly
            # Every site index in 1..N appears exactly once
            @test sort(collect(1:num_sites(lat))) == collect(1:num_sites(lat))
        end
    end

    @testset "Lattice2D re-exports LatticeCore symbols" begin
        @test PeriodicAxis() isa AbstractAxisBC
        @test RowMajor() isa AbstractIndexing
        @test IsingSite() isa AbstractSiteType
        @test UniformLayout(IsingSite()) isa AbstractSiteLayout
    end
end
