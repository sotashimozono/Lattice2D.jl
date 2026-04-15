@testset "sublattice_layout convenience helper" begin
    @testset "Honeycomb: Ising A / XY B" begin
        layout = sublattice_layout(Honeycomb, 4, 4, (IsingSite(), XYSite()))
        lat = build_lattice(Honeycomb, 4, 4; layout=layout)

        for i in 1:num_sites(lat)
            s = sublattice(lat, i)
            st = site_type(site_layout(lat), i)
            if s == 1
                @test st isa IsingSite
            else
                @test st isa XYSite
            end
        end

        # Half the sites are A, half are B.
        a_count = count(i -> site_type(site_layout(lat), i) isa IsingSite,
                        1:num_sites(lat))
        @test a_count == num_sites(lat) ÷ 2
    end

    @testset "Kagome: three distinct site types across sublattices" begin
        layout = sublattice_layout(
            Kagome, 3, 3, (IsingSite(), PottsSite(3), XYSite())
        )
        lat = build_lattice(Kagome, 3, 3; layout=layout)

        # Every site's site type matches its sublattice.
        type_by_sub = (IsingSite(), PottsSite(3), XYSite())
        for i in 1:num_sites(lat)
            s = sublattice(lat, i)
            st = site_type(site_layout(lat), i)
            @test typeof(st) === typeof(type_by_sub[s])
        end

        # Equal counts per sublattice on a 3×3 sample.
        for k in 1:3
            @test count(
                i -> sublattice(lat, i) == k && typeof(site_type(site_layout(lat), i)) === typeof(type_by_sub[k]),
                1:num_sites(lat),
            ) == 9
        end
    end

    @testset "ColMajor indexing: sublattice_layout respects the indexing kwarg" begin
        layout = sublattice_layout(
            Honeycomb, 3, 3, (IsingSite(), XYSite()); indexing=ColMajor()
        )
        lat = build_lattice(Honeycomb, 3, 3; layout=layout, indexing=ColMajor())
        for i in 1:num_sites(lat)
            s = sublattice(lat, i)
            st = site_type(site_layout(lat), i)
            if s == 1
                @test st isa IsingSite
            else
                @test st isa XYSite
            end
        end
    end

    @testset "Wrong arity throws ArgumentError" begin
        @test_throws ArgumentError sublattice_layout(Honeycomb, 3, 3, (IsingSite(),))
        @test_throws ArgumentError sublattice_layout(
            Kagome, 3, 3, (IsingSite(), XYSite())
        )
    end

    @testset "random_state on a mixed-spin Honeycomb layout" begin
        rng = MersenneTwister(42)
        layout = sublattice_layout(Honeycomb, 4, 4, (IsingSite(), XYSite()))
        lat = build_lattice(Honeycomb, 4, 4; layout=layout)

        # Draw one random state per site using the site's declared type.
        for i in 1:num_sites(lat)
            st = site_type(site_layout(lat), i)
            s = random_state(rng, st)
            @test s isa state_type(st)
        end
    end
end
