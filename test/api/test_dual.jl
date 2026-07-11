@testset "dual lattice construction (dual.jl)" begin
    @testset "dual_topology involution" begin
        pairs = [
            (Lattice2D.Square, Lattice2D.Square),
            (Lattice2D.Triangular, Lattice2D.Honeycomb),
            (Lattice2D.Honeycomb, Lattice2D.Triangular),
            (Lattice2D.Kagome, Lattice2D.Dice),
            (Lattice2D.Dice, Lattice2D.Kagome),
        ]
        for (T, D) in pairs
            @test dual_topology(T) === D
            @test dual_topology(dual_topology(T)) === T        # involution
        end
        # instance form
        @test dual_topology(Lattice2D.Triangular()) === Lattice2D.Honeycomb()
    end

    @testset "unsupported topologies throw" begin
        for T in (Lattice2D.Lieb, Lattice2D.ShastrySutherland, Lattice2D.UnionJack)
            @test_throws ArgumentError dual_topology(T)
        end
        @test_throws ArgumentError dual_lattice(lieb(4, 4))
    end

    @testset "triangular <-> honeycomb" begin
        tri = triangular(6, 6)
        hexy = dual_lattice(tri)
        @test topology(hexy) == topology(honeycomb(6, 6))
        @test (hexy.Lx, hexy.Ly) == (6, 6)
        @test num_sublattices(hexy) == 2
        @test num_sites(hexy) == 6 * 6 * 2
        # independent physical check: periodic coordination numbers
        @test mean_coordination(tri) == 6.0
        @test mean_coordination(hexy) == 3.0
        # dual of dual returns to the original topology / size
        back = dual_lattice(hexy)
        @test topology(back) == topology(tri)
        @test num_sites(back) == num_sites(tri)
    end

    @testset "square is self-dual" begin
        sq = square(5, 7)
        d = dual_lattice(sq)
        @test topology(d) == topology(sq)
        @test (d.Lx, d.Ly) == (5, 7)
        @test num_sites(d) == num_sites(sq)
    end

    @testset "kagome <-> dice" begin
        kag = kagome(4, 4)
        di = dual_lattice(kag)
        @test topology(di) == topology(dice(4, 4))
        @test num_sublattices(di) == 3
        # dice mean coordination = (6 + 3 + 3)/3 = 4, matching kagome's 4
        @test mean_coordination(kag) == 4.0
        @test mean_coordination(di) == 4.0
        @test topology(dual_lattice(di)) == topology(kag)
    end

    @testset "boundary and indexing are preserved" begin
        tri = triangular(5, 5; boundary=OpenAxis())
        hexy = dual_lattice(tri)
        @test boundary(hexy) == boundary(tri)
        # open boundary => the dual has under-coordinated (boundary) sites
        @test minimum(coordination_number(hexy)) < 3
        # override still works: fully periodic honeycomb is uniformly 3-coordinated
        per = dual_lattice(tri; boundary=PeriodicAxis())
        @test minimum(coordination_number(per)) == 3
    end
end
