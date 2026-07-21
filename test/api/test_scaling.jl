using Lattice2D
using LatticeCore
using Test

@testset "scaling_rule" begin
    for T in AVAILABLE_LATTICES
        @test scaling_rule(build_lattice(T, 2, 2)) === LinearScaling(2)
    end
end

@testset "rescale multiplies the cell counts" begin
    lat = build_lattice(Square, 2, 3)
    @test rescale(lat, 0) === lat
    @test (rescale(lat).Lx, rescale(lat).Ly) == (4, 6)
    @test (rescale(lat, 2).Lx, rescale(lat, 2).Ly) == (8, 12)

    # topology, boundary, indexing and layout survive the step
    up = rescale(lat, 2)
    @test up isa typeof(lat)
    @test up.boundary === lat.boundary
    @test up.indexing === lat.indexing
    @test up.layout === lat.layout

    # site counts scale as factor^(2k)
    @test num_sites(rescale(lat, 1)) == 4 * num_sites(lat)
    @test num_sites(rescale(lat, 2)) == 16 * num_sites(lat)
end

@testset "rescale down, and refusing to round" begin
    lat = build_lattice(Square, 8, 4)
    @test (rescale(lat, -1).Lx, rescale(lat, -1).Ly) == (4, 2)
    @test (rescale(lat, -2).Lx, rescale(lat, -2).Ly) == (2, 1)
    # 8x4 cannot be halved three times: an inexact step raises rather than rounding
    @test_throws ArgumentError rescale(lat, -3)
    @test_throws ArgumentError rescale(build_lattice(Square, 3, 3), -1)
end

@testset "size_sequence over every topology" begin
    for T in AVAILABLE_LATTICES
        lat = build_lattice(T, 2, 2)
        seq = size_sequence(lat, 2)
        @test length(seq) == 3
        @test num_sites(seq[1]) == num_sites(lat)
        @test num_sites(seq[2]) == 4 * num_sites(lat)
        @test num_sites(seq[3]) == 16 * num_sites(lat)
        # a multi-sublattice topology keeps its sublattice count under rescaling
        @test num_sublattices(seq[3]) == num_sublattices(lat)
    end
end

@testset "cell_partition really partitions" begin
    for T in AVAILABLE_LATTICES, (Lx, Ly) in ((4, 4), (8, 4))
        lat = build_lattice(T, Lx, Ly)
        n = num_sites(lat)
        nsub = num_sublattices(lat)
        for k in 1:2
            (Lx % 2^k == 0 && Ly % 2^k == 0) || continue
            groups = cell_partition(lat, k)
            # covers every site exactly once
            @test sort(reduce(vcat, groups)) == collect(1:n)
            # one group per cell of the coarser lattice
            coarse = rescale(lat, -k)
            @test length(groups) == coarse.Lx * coarse.Ly
            # each group holds factor^(2k) cells' worth of sites
            @test all(g -> length(g) == 4^k * nsub, groups)
        end
    end
end

@testset "cell_partition argument checks" begin
    lat = build_lattice(Square, 4, 4)
    @test_throws ArgumentError cell_partition(lat, 0)
    @test_throws ArgumentError cell_partition(lat, -1)
    @test_throws ArgumentError cell_partition(lat, 3)   # 4 is not divisible by 8
end

@testset "indexing schemes agree on the grouping" begin
    # the partition is a geometric statement, so the groups must contain the same SITES regardless of
    # how sites are numbered — only the labels may permute.
    ref = nothing
    for idx in (RowMajor(), ColMajor(), Snake())
        lat = build_lattice(Square, 4, 4; indexing=idx)
        groups = cell_partition(lat, 1)
        # convert each group to the set of lattice coordinates it covers
        as_coords = Set(
            Set(lattice_coord(idx, (4, 4), 1, i).cell for i in g) for g in groups
        )
        ref === nothing ? (ref = as_coords) : (@test as_coords == ref)
    end
end
