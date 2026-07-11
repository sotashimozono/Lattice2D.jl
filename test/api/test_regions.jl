@testset "edge/bulk region API (regions.jl)" begin
    @testset "open 1D chain: ends are the edge" begin
        chain = square(10, 1; boundary=OpenAxis())
        @test edge_sites(chain) == [1, 10]
        @test edge_sites(chain; depth=2) == [1, 2, 9, 10]
        @test edge_sites(chain; depth=3) == [1, 2, 3, 8, 9, 10]
        # depth large enough to swallow the whole chain
        @test edge_sites(chain; depth=5) == collect(1:10)
        # bulk is the exact complement
        @test bulk_sites(chain) == collect(2:9)
        @test sort(vcat(edge_sites(chain), bulk_sites(chain))) == collect(1:10)
        @test isempty(intersect(edge_sites(chain), bulk_sites(chain)))
    end

    @testset "open square LxL: boundary ring" begin
        L = 5
        lat = square(L, L; boundary=OpenAxis())
        N = num_sites(lat)
        edge = edge_sites(lat)
        # geometric expectation: the outer ring of an LxL grid has
        # L^2 - (L-2)^2 = 4L - 4 sites.
        @test length(edge) == 4L - 4
        # every edge site is under-coordinated (degree < 4 in the bulk)
        @test all(coordination_number(lat, i) < 4 for i in edge)
        # every bulk site has full coordination 4
        @test all(coordination_number(lat, i) == 4 for i in bulk_sites(lat))
        @test length(bulk_sites(lat)) == (L - 2)^2
        # depth=2 peels one more ring: 4L-4 + 4(L-2)-4 sites
        @test length(edge_sites(lat; depth=2)) == (4L - 4) + (4 * (L - 2) - 4)
    end

    @testset "fully periodic sample: no edge" begin
        for lat in (square(6, 6), triangular(6, 6), honeycomb(4, 4))
            @test isempty(edge_sites(lat))
            @test isempty(edge_sites(lat; depth=3))
            @test bulk_sites(lat) == collect(1:num_sites(lat))
            @test isempty(edge_bonds(lat))
        end
    end

    @testset "edge_bonds: incident to an edge site" begin
        chain = square(6, 1; boundary=OpenAxis())
        edge = Set(edge_sites(chain))            # {1, 6}
        eb = edge_bonds(chain)
        @test !isempty(eb)
        # every returned bond touches the edge
        @test all(b.i ∈ edge || b.j ∈ edge for b in eb)
        # a purely-interior bond (e.g. 3-4) must not appear
        @test all(!(b.i == 3 && b.j == 4) && !(b.i == 4 && b.j == 3) for b in eb)
        # bonds touching the two ends: (1-2) and (5-6) => exactly 2
        @test length(eb) == 2
    end

    @testset "cylinder (one periodic, one open axis)" begin
        # periodic along x, open along y: the two open rows are the edge.
        lat = square(6, 4; boundary=LatticeBoundary((PeriodicAxis(), OpenAxis())))
        edge = edge_sites(lat)
        # sites on the two open-boundary rows are under-coordinated (deg 3),
        # interior rows keep full coordination 4.
        @test all(coordination_number(lat, i) < 4 for i in edge)
        @test all(coordination_number(lat, i) == 4 for i in bulk_sites(lat))
        @test !isempty(edge)
    end

    @testset "argument validation" begin
        lat = square(4, 4; boundary=OpenAxis())
        @test_throws ArgumentError edge_sites(lat; depth=0)
        @test_throws ArgumentError bulk_sites(lat; depth=-1)
    end
end
