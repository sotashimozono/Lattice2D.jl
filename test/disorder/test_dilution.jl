@testset "Dilution: site/bond random disorder wrappers" begin
    @testset "dilute_sites with p=0 keeps every site" begin
        base = build_lattice(Square, 5, 5)
        rng = MersenneTwister(0)
        lat = dilute_sites(base, 0.0; rng=rng)

        @test lat isa DilutedLattice
        @test num_sites(lat) == num_sites(base)
        @test all(lat.active_sites)
        @test isempty(lat.killed_bonds)
        # No reindexing happened (1:N → 1:N)
        @test lat.new_to_old == collect(1:num_sites(base))
        @test lat.old_to_new == collect(1:num_sites(base))
        # Bond multiset matches the base lattice (with new == old indices)
        @test length(bonds(lat)) == length(bonds(base))
        @test num_elements(lat, BondCenter()) == num_elements(base, BondCenter())
        # Per-site neighbours agree
        for i in 1:num_sites(lat)
            @test sort(neighbors(lat, i)) == sort(neighbors(base, i))
        end
    end

    @testset "dilute_sites with p=1 removes every site" begin
        base = build_lattice(Square, 4, 4)
        lat = dilute_sites(base, 1.0; rng=MersenneTwister(0))
        @test num_sites(lat) == 0
        @test isempty(bonds(lat))
        @test num_elements(lat, VertexCenter()) == 0
        @test num_elements(lat, BondCenter()) == 0
    end

    @testset "dilute_sites with p=0.3 thins by ~30 %" begin
        base = build_lattice(Square, 20, 20)
        N = num_sites(base)
        rng = MersenneTwister(42)
        lat = dilute_sites(base, 0.3; rng=rng)
        # With 400 i.i.d. Bernoulli(0.7) survivals, the count is highly
        # concentrated; allow a generous ±15 % window so the test stays
        # green across RNG implementations.
        target = round(Int, 0.7N)
        @test abs(num_sites(lat) - target) ≤ round(Int, 0.15N)
        # All neighbours of the diluted lattice must be active sites.
        for i in 1:num_sites(lat)
            for j in neighbors(lat, i)
                @test j != i
                @test 1 ≤ j ≤ num_sites(lat)
            end
        end
        # `bonds` must list every alive endpoint pair at most once and
        # respect the canonical i < j convention for finite lattices.
        seen = Set{Tuple{Int,Int}}()
        for b in bonds(lat)
            key = (min(b.i, b.j), max(b.i, b.j))
            @test !(key in seen)
            push!(seen, key)
        end
    end

    @testset "dilute_sites is reproducible with the same RNG" begin
        base = build_lattice(Square, 6, 6)
        a = dilute_sites(base, 0.4; rng=MersenneTwister(123))
        b = dilute_sites(base, 0.4; rng=MersenneTwister(123))
        @test a.active_sites == b.active_sites
        @test num_sites(a) == num_sites(b)
    end

    @testset "dilute_bonds keeps every site, drops bonds" begin
        base = build_lattice(Square, 6, 6)
        N = num_sites(base)
        Bcount = length(collect(bonds(base)))

        rng = MersenneTwister(7)
        lat = dilute_bonds(base, 0.4; rng=rng)

        @test num_sites(lat) == N
        @test all(lat.active_sites)
        @test length(lat.killed_bonds) > 0
        # New index space is unchanged (no remap)
        @test lat.new_to_old == collect(1:N)
        # Position is identical to base for every site
        for i in 1:N
            @test position(lat, i) == position(base, i)
        end
        # Bond count drops by exactly the number of explicitly killed
        # bonds (every base bond either survives or is killed).
        @test length(bonds(lat)) == Bcount - length(lat.killed_bonds)
        # Killed bonds shouldn't appear in the iterator
        for (a, c) in lat.killed_bonds
            for b in bonds(lat)
                key = (min(b.i, b.j), max(b.i, b.j))
                @test key != (a, c)
            end
        end
    end

    @testset "dilute_bonds with p=0 / p=1 endpoints" begin
        base = build_lattice(Square, 3, 4)
        Bcount = length(collect(bonds(base)))

        keep = dilute_bonds(base, 0.0; rng=MersenneTwister(1))
        @test length(bonds(keep)) == Bcount
        @test isempty(keep.killed_bonds)

        cut = dilute_bonds(base, 1.0; rng=MersenneTwister(1))
        @test isempty(bonds(cut))
        @test length(cut.killed_bonds) == Bcount
        # Sites are still all present, just isolated
        for i in 1:num_sites(cut)
            @test isempty(neighbors(cut, i))
            @test isempty(neighbor_bonds(cut, i))
        end
    end

    @testset "neighbor_bonds endpoints carry the new (remapped) index" begin
        base = build_lattice(Square, 5, 5)
        rng = MersenneTwister(99)
        lat = dilute_sites(base, 0.25; rng=rng)
        for i in 1:num_sites(lat)
            for b in neighbor_bonds(lat, i)
                @test b.i == i
                @test 1 ≤ b.j ≤ num_sites(lat)
                @test b.j != i
                @test i in neighbors(lat, b.j) || b.j in neighbors(lat, i)
            end
        end
    end

    @testset "input validation" begin
        base = build_lattice(Square, 3, 3)
        @test_throws ArgumentError dilute_sites(base, -0.1)
        @test_throws ArgumentError dilute_sites(base, 1.1)
        @test_throws ArgumentError dilute_bonds(base, -0.5)
        @test_throws ArgumentError dilute_bonds(base, 2.0)

        # active_sites length mismatch
        @test_throws ArgumentError DilutedLattice(base, BitVector([true, false]))

        # killed_bonds out-of-range / self-loop
        N = num_sites(base)
        @test_throws ArgumentError DilutedLattice(base, trues(N), Set([(1, N + 1)]))
        @test_throws ArgumentError DilutedLattice(base, trues(N), Set([(2, 2)]))
    end

    @testset "size_trait / boundary / sublattice forwarded through wrapper" begin
        base = build_lattice(Honeycomb, 3, 3)
        lat = dilute_sites(base, 0.0; rng=MersenneTwister(0))
        @test size_trait(lat) === size_trait(base)
        @test boundary(lat) === boundary(base)
        @test num_sublattices(lat) == num_sublattices(base)
        for i in 1:num_sites(lat)
            @test sublattice(lat, i) == sublattice(base, lat.new_to_old[i])
        end
        # Diluted lattices report a `:disordered` topology trait so
        # generic Bravais-aware code does not silently misuse them.
        @test topology(lat) === TopologyTrait{:disordered}()
        @test periodicity(lat) === Aperiodic()
        @test reciprocal_support(lat) === NoReciprocal()
    end
end
