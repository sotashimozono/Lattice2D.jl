using Test
using Lattices

@testset "Iterator & Indexing Interfaces" begin
    Lx, Ly = 3, 3
    n_sub = 2
    lat = build_lattice(Honeycomb, Lx, Ly)
    @testset "Basic Iteration" begin
        @test length(lat) == lat.N
        @test eltype(lat) == Int
        @test collect(lat) == 1:lat.N
        @test sum(lat) == lat.N * (lat.N + 1) ÷ 2
    end

    @testset "Neighbors Iterator" begin
        i = 1
        nbrs = neighbors(lat, i)
        
        @test nbrs isa AbstractVector{Int}
        
        @test nbrs == lat.nearest_neighbors[i]
        
        count = 0
        for n in neighbors(lat, i)
            count += 1
            @test n isa Int
        end
        @test count == length(lat.nearest_neighbors[i])
    end
end