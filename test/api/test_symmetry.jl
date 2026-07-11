using Lattice2D
using LatticeCore
using LinearAlgebra
using Test

@testset "point-group symmetry & orbits (symmetry.jl)" begin
    @testset "realized group orders match the known point groups" begin
        @test symmetry_group_order(build_lattice(Square, 6, 6)) == 8          # C4v
        @test symmetry_group_order(build_lattice(Triangular, 6, 6)) == 12     # C6v
        @test symmetry_group_order(build_lattice(Honeycomb, 4, 4)) == 6       # site C3v
        @test symmetry_group_order(build_lattice(Square, 5, 5; boundary=OpenAxis())) == 8
    end

    @testset "C4v composition: 4 rotations + 4 reflections, inversion present" begin
        ops = symmetry_operations(build_lattice(Square, 6, 6))
        @test count(is_rotation, ops) == 4
        @test count(is_reflection, ops) == 4
        @test any(is_inversion, ops)                      # 180° ∈ C4v
        # rotation angles present: 0, ±90, 180 (mod 2π)
        angs = sort([round(rotation_angle(op); digits=6) for op in ops if is_rotation(op)])
        @test angs ≈ sort(round.([0.0, π / 2, -π / 2, π]; digits=6))
    end

    @testset "group axioms (identity, bijectivity, closure, inverses)" begin
        for lat in (build_lattice(Square, 6, 6), build_lattice(Triangular, 6, 6))
            ops = symmetry_operations(lat)
            N = num_sites(lat)
            id = 1:N

            # identity present exactly once
            @test count(op -> op.permutation == collect(id), ops) == 1
            # every operation is a genuine permutation of the sites
            @test all(sort(op.permutation) == collect(id) for op in ops)

            permset = Set(op.permutation for op in ops)
            # closure: composing any two operations lands back in the group
            for a in ops, b in ops
                comp = [a.permutation[b.permutation[i]] for i in id]
                @test comp in permset
            end
            # inverse of each operation is in the group
            for a in ops
                inv = Vector{Int}(undef, N)
                for i in id
                    inv[a.permutation[i]] = i
                end
                @test inv in permset
            end
        end
    end

    @testset "orbits: partition + orbit–stabilizer divisibility" begin
        for lat in (
            build_lattice(Square, 6, 6),
            build_lattice(Triangular, 6, 6),
            build_lattice(Honeycomb, 4, 4),
            build_lattice(Square, 5, 5; boundary=OpenAxis()),
        )
            N = num_sites(lat)
            G = symmetry_group_order(lat)
            orbits = symmetry_orbits(lat)
            # partition: disjoint cover of 1:N
            @test sort(vcat(orbits...)) == collect(1:N)
            # each orbit sorted; orbit sizes divide |G| (orbit–stabilizer)
            @test all(issorted(o) for o in orbits)
            @test all(G % length(o) == 0 for o in orbits)
            # site_orbit agrees with the containing orbit
            for i in (1, N ÷ 2 + 1, N)
                o = site_orbit(lat, i)
                @test i in o
                @test o in orbits
            end
        end
    end

    @testset "origin site is a fixed point of the periodic point group" begin
        lat = build_lattice(Square, 6, 6)
        @test position(lat, 1) ≈ [0.0, 0.0]         # site 1 at the rotation centre
        @test site_orbit(lat, 1) == [1]              # fixed by all of C4v
    end

    @testset "honeycomb A/B sublattices never mix under the site group" begin
        lat = build_lattice(Honeycomb, 4, 4)
        for o in symmetry_orbits(lat)
            subs = unique(sublattice(lat, i) for i in o)
            @test length(subs) == 1                  # orbits respect the sublattice
        end
    end
end
