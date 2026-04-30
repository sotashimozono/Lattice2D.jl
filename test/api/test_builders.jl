@testset "Convenience builders" begin
    # (builder, Topology, num_sublattices) — covers all 8 shipped topologies.
    cases = [
        (square, Square, 1),
        (triangular, Triangular, 1),
        (honeycomb, Honeycomb, 2),
        (kagome, Kagome, 3),
        (lieb, Lieb, 3),
        (union_jack, UnionJack, 2),
        (dice, Dice, 3),
        (shastry_sutherland, ShastrySutherland, 4),
    ]

    @testset "covers every shipped topology" begin
        # The set of topologies wired through the convenience layer must
        # exactly match `AVAILABLE_LATTICES` so that adding a new topology
        # to the package without a builder fails this test.
        @test Set(c[2] for c in cases) == Set(AVAILABLE_LATTICES)
    end

    @testset "default M = L: $(nameof(builder))" for (builder, Topology, nsub) in cases
        L = 4
        lat = builder(L)
        @test lat isa Lattice{Topology}
        @test num_sublattices(lat) == nsub
        @test num_sites(lat) == L * L * nsub
        # Default boundary is fully periodic on both axes.
        @test periodicity(lat) isa Periodic
    end

    @testset "explicit M ≠ L: $(nameof(builder))" for (builder, Topology, nsub) in cases
        L, M = 3, 5
        lat = builder(L, M)
        @test lat isa Lattice{Topology}
        @test num_sites(lat) == L * M * nsub
    end

    @testset "kw forwarding: $(nameof(builder))" for (builder, Topology, _) in cases
        # `boundary` and `layout` keywords must be passed through to
        # `build_lattice` unchanged.
        lat = builder(3, 3; boundary=OpenAxis(), layout=UniformLayout(XYSite()))
        ref = build_lattice(
            Topology, 3, 3; boundary=OpenAxis(), layout=UniformLayout(XYSite())
        )
        @test lat isa Lattice{Topology}
        @test num_sites(lat) == num_sites(ref)
        @test site_type(lat, 1) === XYSite()
        @test boundary(lat) == boundary(ref)
    end

    @testset "delegation equivalence: $(nameof(builder))" for (builder, Topology, _) in
                                                              cases
        # The shortcut and the long form must produce lattices with the
        # same observable topology / size / connectivity.
        lat = builder(4, 4)
        ref = build_lattice(Topology, 4, 4)
        @test lat isa Lattice{Topology}
        @test num_sites(lat) == num_sites(ref)
        @test num_sublattices(lat) == num_sublattices(ref)
        @test count(_ -> true, bonds(lat)) == count(_ -> true, bonds(ref))
    end
end
