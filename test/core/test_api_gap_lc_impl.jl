@testset "LatticeCore API gap implementations on Lattice" begin

    # ---- to_lattice (issue #41) ------------------------------------

    @testset "to_lattice round-trips with to_real on every shipped topology" begin
        for Topo in (
            Square, Triangular, Honeycomb, Kagome, Lieb, ShastrySutherland, Dice, UnionJack
        )
            lat = build_lattice(Topo, 4, 4)
            for i in 1:num_sites(lat)
                r = to_real(lat, lattice_coord(RowMajor(), (4, 4), num_sublattices(lat), i))
                back = to_lattice(lat, r)
                @test back == lattice_coord(RowMajor(), (4, 4), num_sublattices(lat), i)
            end
        end
    end

    @testset "to_lattice resolves sublattice on Honeycomb" begin
        lat = build_lattice(Honeycomb, 3, 3)
        for i in 1:num_sites(lat)
            p = position(lat, i)
            lc = to_lattice(lat, RealSpace{2,Float64}(SVector(p[1], p[2])))
            @test lc.sublattice == sublattice(lat, i)
        end
    end

    @testset "to_lattice on a slightly perturbed point still rounds correctly" begin
        # Each site of a 4 × 4 Square lattice sits on integer coords
        # (0, 0), (1, 0), .... A small Δ smaller than half the basis
        # vector must still resolve to the same site.
        lat = build_lattice(Square, 4, 4)
        Δ = SVector(0.1, -0.1)
        for i in 1:num_sites(lat)
            p = position(lat, i)
            lc = to_lattice(lat, RealSpace{2,Float64}(p + Δ))
            i_back = site_index(RowMajor(), (4, 4), 1, lc)
            @test i_back == i
        end
    end

    @testset "to_lattice wraps under PBC, clips to candidate cells under OBC" begin
        # PBC: a point one full lattice vector to the right of site 1
        # must land back on site 1's cell because of the wrap.
        lat_pbc = build_lattice(Square, 4, 4)
        a1, a2 = basis_vectors(lat_pbc)[:, 1], basis_vectors(lat_pbc)[:, 2]
        p1 = position(lat_pbc, 1)
        lc = to_lattice(lat_pbc, RealSpace{2,Float64}(p1 + 4 * a1))
        @test lc.cell == (1, 1)

        # OBC: same query on an open lattice does NOT wrap — it returns
        # the nearest unwrapped cell index even though it's out of range.
        # We only assert that the query does not throw and that the
        # returned cell agrees with the unwrapped rounding.
        lat_obc = build_lattice(Square, 4, 4; boundary=OpenAxis())
        lc_obc = to_lattice(lat_obc, RealSpace{2,Float64}(p1 + 4 * a1))
        @test lc_obc.cell == (5, 1)
    end

    @testset "to_lattice(::LatticeCoord) is the identity" begin
        # Inherited from LatticeCore, but exercise it through Lattice2D.
        lat = build_lattice(Square, 3, 3)
        c = LatticeCoord{2}((2, 1), 1)
        @test to_lattice(lat, c) === c
    end

    # ---- element_position(::VertexCenter) (issue #42) ---------------

    @testset "element_position(::VertexCenter) matches position(lat, i)" begin
        for Topo in (Square, Honeycomb, Kagome, Lieb, ShastrySutherland)
            lat = build_lattice(Topo, 3, 3)
            for i in 1:num_sites(lat)
                @test element_position(lat, VertexCenter(), i) == position(lat, i)
            end
        end
    end

    # ---- element_positions specialisations (issue #43) --------------

    @testset "element_positions(::VertexCenter) iterates positions(lat)" begin
        lat = build_lattice(Square, 4, 3)
        ps = collect(element_positions(lat, VertexCenter()))
        @test length(ps) == num_sites(lat)
        @test ps == collect(positions(lat))
    end

    @testset "element_positions(::BondCenter) lazy: matches element_position pointwise" begin
        lat = build_lattice(Triangular, 4, 4)
        nb = num_elements(lat, BondCenter())
        ps = collect(element_positions(lat, BondCenter()))
        @test length(ps) == nb
        for i in 1:nb
            @test ps[i] == element_position(lat, BondCenter(), i)
        end
    end

    @testset "element_positions(::PlaquetteCenter) lazy: matches plaquette centers" begin
        lat = build_lattice(Honeycomb, 3, 3)
        np = num_elements(lat, PlaquetteCenter())
        ps = collect(element_positions(lat, PlaquetteCenter()))
        @test length(ps) == np
        for i in 1:np
            @test ps[i] == element_position(lat, PlaquetteCenter(), i)
        end
    end

    @testset "element_positions does not allocate a Vector{...} eagerly" begin
        # The specialised methods return generators; they should not
        # be `<:AbstractVector`. The user-facing contract is
        # "iterator over positions" — concrete tests for that contract
        # live above; here we just confirm the specialisation didn't
        # accidentally regress to eager `collect`.
        lat = build_lattice(Kagome, 3, 3)
        @test !(element_positions(lat, BondCenter()) isa AbstractVector)
        @test !(element_positions(lat, PlaquetteCenter()) isa AbstractVector)
    end
end
