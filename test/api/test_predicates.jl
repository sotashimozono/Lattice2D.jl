using LinearAlgebra: norm

@testset "lattice predicates" begin
    @testset "is_bipartite" begin
        # Bipartite topologies (PBC, even Lx/Ly so the colouring closes)
        @test is_bipartite(build_lattice(Square, 4, 4)) == true
        @test is_bipartite(build_lattice(Honeycomb, 4, 4)) == true
        @test is_bipartite(build_lattice(Honeycomb, 6, 6)) == true
        @test is_bipartite(build_lattice(Lieb, 4, 4)) == true
        @test is_bipartite(build_lattice(Dice, 4, 4)) == true

        # Non-bipartite topologies
        @test is_bipartite(build_lattice(Triangular, 4, 4)) == false
        @test is_bipartite(build_lattice(Triangular, 6, 6)) == false
        @test is_bipartite(build_lattice(Kagome, 4, 4)) == false

        # Square with an odd-cycle wraparound (PBC, odd Lx)
        @test is_bipartite(build_lattice(Square, 3, 4)) == false
        # OBC removes the wraparound, restoring bipartiteness
        @test is_bipartite(build_lattice(Square, 3, 3; boundary=OpenAxis())) == true
    end

    @testset "coordination_number — bulk values, PBC" begin
        # Bulk coordination on uniform lattices (PBC: all sites equivalent
        # for Bravais; sublattice lattices have one number per sublattice
        # but share the same set, which we summarise via `unique` here).
        @test all(==(4), coordination_number(build_lattice(Square, 6, 6)))
        @test all(==(6), coordination_number(build_lattice(Triangular, 6, 6)))
        @test all(==(3), coordination_number(build_lattice(Honeycomb, 4, 4)))
        @test all(==(4), coordination_number(build_lattice(Kagome, 4, 4)))

        # Single-site form agrees with the vector form
        lat = build_lattice(Square, 5, 5)
        v = coordination_number(lat)
        @test length(v) == num_sites(lat)
        @test v[1] == coordination_number(lat, 1)
        @test v[end] == coordination_number(lat, num_sites(lat))
    end

    @testset "coordination_number — OBC drops at the boundary" begin
        lat = build_lattice(Square, 4, 4; boundary=OpenAxis())
        v = coordination_number(lat)
        # Corners have 2 neighbours, edges 3, interior 4.
        @test 2 in v
        @test 3 in v
        @test 4 in v
        @test minimum(v) == 2
        @test maximum(v) == 4
    end

    @testset "mean_coordination" begin
        # PBC: mean equals the bulk coordination.
        @test mean_coordination(build_lattice(Square, 6, 6)) ≈ 4.0
        @test mean_coordination(build_lattice(Triangular, 6, 6)) ≈ 6.0
        @test mean_coordination(build_lattice(Honeycomb, 4, 4)) ≈ 3.0
        @test mean_coordination(build_lattice(Kagome, 4, 4)) ≈ 4.0

        # OBC: strictly below the bulk value.
        lat_obc = build_lattice(Square, 4, 4; boundary=OpenAxis())
        @test mean_coordination(lat_obc) < 4.0

        # Consistency with the vector form.
        v = coordination_number(lat_obc)
        @test mean_coordination(lat_obc) ≈ sum(v) / length(v)
    end

    @testset "bond_distances" begin
        lat = build_lattice(Square, 4, 4)
        d = bond_distances(lat)
        @test length(d) == length(bonds(lat))
        # Square NN bond length is 1.
        @test all(x -> isapprox(x, 1.0; atol=1e-12), d)

        # Honeycomb has a single NN distance (1/√3 in our convention).
        lat_h = build_lattice(Honeycomb, 4, 4)
        dh = bond_distances(lat_h)
        @test length(dh) == length(bonds(lat_h))
        @test all(x -> isapprox(x, dh[1]; atol=1e-10), dh)

        # Cross-check against the underlying Bond.vector.
        for (b, x) in zip(bonds(lat_h), dh)
            @test x ≈ Float64(norm(b.vector))
        end
    end

    @testset "shells: agrees with neighbors(...; shell=k)" begin
        lat = build_lattice(Triangular, 6, 6)
        sh = shells(lat, 1; n_shells=3)
        @test length(sh) == 3
        @test sort(sh[1]) == sort(neighbors(lat, 1; shell=1))
        @test sort(sh[2]) == sort(neighbors(lat, 1; shell=2))
        @test sort(sh[3]) == sort(neighbors(lat, 1; shell=3))
        # Shells are mutually disjoint.
        @test isempty(intersect(Set(sh[1]), Set(sh[2])))
        @test isempty(intersect(Set(sh[1]), Set(sh[3])))
        @test isempty(intersect(Set(sh[2]), Set(sh[3])))
    end

    @testset "shells: default n_shells = 3" begin
        lat = build_lattice(Square, 6, 6)
        @test length(shells(lat, 1)) == 3
    end

    @testset "shells: trailing empty for over-deep request on small OBC" begin
        lat = build_lattice(Square, 3, 3; boundary=OpenAxis())
        sh = shells(lat, 1; n_shells=10)
        @test length(sh) == 10
        # The deepest shells must be empty on a 3×3 OBC sample.
        @test isempty(sh[end])
    end

    @testset "shells: n_shells must be ≥ 1" begin
        lat = build_lattice(Square, 3, 3)
        @test_throws ArgumentError shells(lat, 1; n_shells=0)
        @test_throws ArgumentError shells(lat, 1; n_shells=-1)
    end
end
