@testset "Base iteration / size / eltype on PeriodicLattice2D" begin
    Lx, Ly = 3, 3
    lat = build_lattice(Honeycomb, Lx, Ly)

    @testset "length / size / eltype" begin
        @test length(lat) == num_sites(lat)
        @test size(lat) == (Lx, Ly)
        @test size(lat, 1) == Lx
        @test size(lat, 2) == Ly
        @test size(lat, 3) == 1
        @test eltype(lat) == Int
    end

    @testset "collect / sum / iterate cover 1..N" begin
        collected = collect(lat)
        @test collected == collect(1:num_sites(lat))
        @test sum(lat) == num_sites(lat) * (num_sites(lat) + 1) ÷ 2
    end

    @testset "neighbors returns an AbstractVector{Int}" begin
        nbrs = neighbors(lat, 1)
        @test nbrs isa AbstractVector{Int}
        for n in nbrs
            @test n isa Int
            @test 1 <= n <= num_sites(lat)
        end
    end
end
