@testset "LatticeCore API coverage on Lattice" begin
    @testset "positions(lat) returns every site" begin
        lat = build_lattice(Square, 3, 4)
        ps = collect(positions(lat))
        @test length(ps) == num_sites(lat)
        # Every collected position must agree with `position(lat, i)`.
        for (i, p) in enumerate(ps)
            @test p == position(lat, i)
        end
    end

    @testset "bond_center on PBC: wrapped bonds use unwrapped vector" begin
        lat = build_lattice(Square, 4, 4)
        # Site 1 = (1, 1). Its left neighbour wraps to (4, 1) (= site 4
        # under RowMajor). The bond's stored vector is (-1, 0), so the
        # midpoint sits half a step LEFT of (1, 1), i.e. at (-0.5, 0).
        # The lattice uses 0-based positions internally, so the
        # corner site 1 is at SVector(0.0, 0.0).
        nb = neighbor_bonds(lat, 1)
        wrapped = first(b for b in nb if b.vector[1] == -1.0)
        center = bond_center(lat, wrapped)
        @test center == SVector(-0.5, 0.0)

        # An interior bond (site 6 → site 7 on the same row, dx=+1)
        # should match the literal midpoint.
        nb6 = neighbor_bonds(lat, 6)
        right = first(b for b in nb6 if b.vector == SVector(1.0, 0.0))
        @test bond_center(lat, right) == (position(lat, 6) + position(lat, right.j)) / 2
    end

    @testset "bond_center: every PBC bond's midpoint lies on the bond" begin
        # A robust invariant that doesn't depend on knowing wrap details:
        # the midpoint of every bond must be reachable by walking half
        # the bond's `vector` from its source endpoint.
        lat = build_lattice(Triangular, 4, 4)
        for b in collect(bonds(lat))
            @test bond_center(lat, b) == position(lat, b.i) + b.vector / 2
        end
    end

    @testset "site_type via the site_layout fall-through" begin
        lat = build_lattice(Square, 3, 3)
        layout = site_layout(lat)
        @test layout isa UniformLayout
        @test site_type(layout, 1) === IsingSite{Int8}()
        # Same site type for every site under UniformLayout.
        @test all(site_type(layout, i) === site_type(layout, 1) for i in 1:num_sites(lat))
    end

    @testset "is_finite and size_trait" begin
        lat = build_lattice(Honeycomb, 3, 3)
        @test is_finite(lat)
        @test size_trait(lat) isa FiniteSize{2}
    end

    @testset "momentum_lattice is the reciprocal lattice for Bravais cases" begin
        lat = build_lattice(Triangular, 4, 4)
        ml_rec = reciprocal_lattice(lat)
        ml = momentum_lattice(lat)
        @test ml === ml_rec || (
            num_k_points(ml) == num_k_points(ml_rec) &&
            all(k_point(ml, i) == k_point(ml_rec, i) for i in 1:num_k_points(ml))
        )

        # Mixed-BC samples have no reciprocal lattice → momentum_lattice throws.
        cyl = build_lattice(
            Square, 4, 4; boundary=LatticeBoundary((PeriodicAxis(), OpenAxis()))
        )
        @test_throws ArgumentError momentum_lattice(cyl)
    end

    @testset "structure_factor at Γ on a uniform state = N" begin
        # ⟨s_i⟩ uniform over i ⇒ S(k=0) = N (1/N)|N|² = N.
        lat = build_lattice(Square, 4, 4)
        N = num_sites(lat)
        state = ones(Int, N)
        γ = SVector(0.0, 0.0)
        @test structure_factor(lat, state, γ) ≈ Float64(N)
    end

    @testset "structure_factor over a momentum lattice returns one value per k" begin
        lat = build_lattice(Square, 4, 4)
        ml = reciprocal_lattice(lat)
        state = ones(Int, num_sites(lat))
        sks = structure_factor(lat, state, ml)
        @test length(sks) == num_k_points(ml)
        # On a half-shifted Monkhorst–Pack mesh Γ is not sampled, so we
        # don't assert any specific peak — we only require finite, real,
        # non-negative values.
        @test all(isfinite, sks)
        @test all(sks .≥ 0)

        # On a Γ-centred mesh Γ IS sampled, and a uniform state has
        # S(Γ) = N exactly.
        γ_mesh = gamma_centered(reciprocal_basis(ml), (4, 4))
        sks_γ = structure_factor(lat, state, γ_mesh)
        @test maximum(sks_γ) ≈ Float64(num_sites(lat)) atol = 1e-10
    end

    @testset "site type element_type defaults to VertexCenter" begin
        @test element_type(IsingSite()) === VertexCenter()
        @test element_type(XYSite()) === VertexCenter()
        @test element_type(HeisenbergSite()) === VertexCenter()
    end

    @testset "random_state / state_type round-trips for shipped site types" begin
        rng = MersenneTwister(42)
        for st in (IsingSite(), PottsSite(3), XYSite(), HeisenbergSite())
            s = random_state(rng, st)
            @test s isa state_type(st)
        end
    end
end
