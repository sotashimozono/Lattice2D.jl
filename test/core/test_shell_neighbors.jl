@testset "multi-shell neighbours" begin
    @testset "Square PBC: shells are axial / diagonal / 2-axial" begin
        lat = build_lattice(Square, 6, 6)
        # Shell 1: 4 axial neighbours at distance 1
        @test length(neighbors(lat, 1; shell=1)) == 4
        # Shell 2: 4 diagonal neighbours at distance √2
        @test length(neighbors(lat, 1; shell=2)) == 4
        # Shell 3: 4 next-axial at distance 2
        @test length(neighbors(lat, 1; shell=3)) == 4

        # Declared = shell=1 for Square (NN only)
        @test sort(neighbors(lat, 1)) == sort(neighbors(lat, 1; shell=1))
    end

    @testset "Triangular PBC: six NN / six 2nd-NN / six 3rd-NN" begin
        lat = build_lattice(Triangular, 6, 6)
        @test length(neighbors(lat, 1; shell=1)) == 6
        @test length(neighbors(lat, 1; shell=2)) == 6
        @test length(neighbors(lat, 1; shell=3)) == 6
    end

    @testset "Honeycomb PBC: 3 / 6 / 3 interleaving" begin
        lat = build_lattice(Honeycomb, 4, 4)
        # From an A site: NN are the 3 B partners
        @test length(neighbors(lat, 1; shell=1)) == 3
        # 2nd shell: 6 same-sublattice (triangular) neighbours
        @test length(neighbors(lat, 1; shell=2)) == 6
        # 3rd shell: 3 more B partners at the next Honeycomb distance
        @test length(neighbors(lat, 1; shell=3)) == 3
    end

    @testset "Shastry–Sutherland: declared ≠ geometric shell 1" begin
        lat = build_lattice(ShastrySutherland, 4, 4)
        # Declared: 4 square NN + 1 dimer partner per site
        @test length(neighbors(lat, 1)) == 5
        # Geometric shell 1: only the 4 square NN (distance 1)
        @test length(neighbors(lat, 1; shell=1)) == 4
        # Geometric shell 2: includes the dimer partner and any other
        # sites at the same √2 distance — the dimer partner must be
        # among the returned indices.
        shell2 = neighbors(lat, 1; shell=2)
        declared = neighbors(lat, 1)
        dimer_partner = setdiff(declared, neighbors(lat, 1; shell=1))
        @test length(dimer_partner) == 1
        @test dimer_partner[1] in shell2
    end

    @testset "OBC restricts the shell count for corner sites" begin
        lat = build_lattice(Square, 4, 4; boundary=OpenAxis())
        # Site 1 is the (1, 1) corner: 2 axial neighbours, 1 diagonal.
        @test length(neighbors(lat, 1; shell=1)) == 2
        @test length(neighbors(lat, 1; shell=2)) == 1
    end

    @testset "shell must be ≥ 1" begin
        lat = build_lattice(Square, 3, 3)
        @test_throws ArgumentError neighbors(lat, 1; shell=0)
        @test_throws ArgumentError neighbors(lat, 1; shell=-1)
    end

    @testset "over-deep shell returns empty" begin
        lat = build_lattice(Square, 3, 3)
        # Well beyond any reachable distance on a 3×3 lattice.
        @test isempty(neighbors(lat, 1; shell=100))
    end

    @testset "shells are mutually disjoint" begin
        lat = build_lattice(Triangular, 6, 6)
        s1 = Set(neighbors(lat, 1; shell=1))
        s2 = Set(neighbors(lat, 1; shell=2))
        s3 = Set(neighbors(lat, 1; shell=3))
        @test isempty(intersect(s1, s2))
        @test isempty(intersect(s1, s3))
        @test isempty(intersect(s2, s3))
    end

    @testset "indexing-invariant: shell count is a geometric property" begin
        # The set of shell-1 distances from a PBC Square site is the
        # same regardless of RowMajor / ColMajor / Snake indexing.
        counts = Dict{Symbol,Int}()
        for indexing in (RowMajor(), ColMajor(), Snake())
            lat = build_lattice(Square, 4, 4; indexing=indexing)
            counts[Symbol(typeof(indexing).name.name)] = length(neighbors(lat, 1; shell=1))
        end
        @test length(unique(values(counts))) == 1
    end
end
