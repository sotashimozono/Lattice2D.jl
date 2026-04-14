@testset "Reciprocal vectors" begin
    @testset "reciprocal_lattice on every topology passes orthogonality" begin
        for Topo in AVAILABLE_LATTICES
            lat = build_lattice(Topo, 4, 4)
            A = basis_vectors(lat)
            ml = reciprocal_lattice(lat)
            B = ml.reciprocal_basis
            # a_i · b_j = 2π δ_ij
            for i in 1:2, j in 1:2
                expected = i == j ? 2π : 0.0
                @test dot(A[:, i], B[:, j]) ≈ expected atol = 1e-10
            end
        end
    end

    @testset "reciprocal_lattice throws under non-periodic boundary" begin
        lat_obc = build_lattice(Square, 4, 4; boundary=OpenAxis())
        @test_throws ArgumentError reciprocal_lattice(lat_obc)

        cyl = build_lattice(
            Square, 4, 4; boundary=LatticeBoundary((PeriodicAxis(), OpenAxis()))
        )
        @test_throws ArgumentError reciprocal_lattice(cyl)
    end
end
