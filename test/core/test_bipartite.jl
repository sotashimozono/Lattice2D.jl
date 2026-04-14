@testset "Bipartiteness of shipped topologies" begin
    @testset "Square: bipartite iff both axis lengths are even" begin
        @test is_bipartite(build_lattice(Square, 4, 4)) == true
        @test is_bipartite(build_lattice(Square, 3, 4)) == false
        @test is_bipartite(build_lattice(Square, 3, 3)) == false

        # Open chain is always bipartite
        @test is_bipartite(build_lattice(Square, 3, 3; boundary=OpenAxis())) == true
    end

    @testset "Triangular is never bipartite under PBC" begin
        @test is_bipartite(build_lattice(Triangular, 4, 4)) == false
        @test is_bipartite(build_lattice(Triangular, 6, 6)) == false
    end

    @testset "Honeycomb is bipartite (A/B sublattices)" begin
        @test is_bipartite(build_lattice(Honeycomb, 4, 4)) == true
        @test is_bipartite(build_lattice(Honeycomb, 6, 6)) == true
    end

    @testset "Kagome contains triangles: not bipartite" begin
        @test is_bipartite(build_lattice(Kagome, 4, 4)) == false
    end

    @testset "Lieb lattice is bipartite" begin
        @test is_bipartite(build_lattice(Lieb, 4, 4)) == true
    end
end
