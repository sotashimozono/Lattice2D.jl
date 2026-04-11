@testset "Reciprocal vectors" begin
    @testset "calc_reciprocal_vectors helper" begin
        basis_ortho = [[1.0, 0.0], [0.0, 1.0]]
        recip_ortho = calc_reciprocal_vectors(basis_ortho)
        @test recip_ortho[1] ≈ [2π, 0.0]
        @test recip_ortho[2] ≈ [0.0, 2π]
        @test dot(basis_ortho[1], recip_ortho[1]) ≈ 2π
        @test dot(basis_ortho[2], recip_ortho[2]) ≈ 2π
        @test abs(dot(basis_ortho[1], recip_ortho[2])) < 1e-10
        @test abs(dot(basis_ortho[2], recip_ortho[1])) < 1e-10

        # Non-orthogonal: triangular basis
        a1 = [1.0, 0.0]
        a2 = [0.5, sqrt(3) / 2]
        recip_tri = calc_reciprocal_vectors([a1, a2])
        @test dot(a1, recip_tri[1]) ≈ 2π
        @test dot(a2, recip_tri[2]) ≈ 2π
        @test abs(dot(a1, recip_tri[2])) < 1e-10
        @test abs(dot(a2, recip_tri[1])) < 1e-10
    end

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
