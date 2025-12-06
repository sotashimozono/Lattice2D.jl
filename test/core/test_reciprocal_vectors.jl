@testset "calc_reciprocal_vectors" begin
    basis_ortho = [[1.0, 0.0], [0.0, 1.0]]
    recip_ortho = calc_reciprocal_vectors(basis_ortho)
    @test recip_ortho[1] ≈ [2π, 0.0]
    @test recip_ortho[2] ≈ [0.0, 2π]
    @test dot(basis_ortho[1], recip_ortho[1]) ≈ 2π
    @test dot(basis_ortho[2], recip_ortho[2]) ≈ 2π
    @test abs(dot(basis_ortho[1], recip_ortho[2])) < 1e-10
    @test abs(dot(basis_ortho[2], recip_ortho[1])) < 1e-10
    # Case 2: Non-Orthogonal Basis (Triangular-like)
    # a1 = [1, 0], a2 = [0.5, √3/2]
    a1 = [1.0, 0.0]
    a2 = [0.5, sqrt(3)/2]
    basis_tri = [a1, a2]
    recip_tri = calc_reciprocal_vectors(basis_tri)
    
    b1, b2 = recip_tri
    # Definition Check: a_i ⋅ b_j = 2π δ_ij
    @test dot(a1, b1) ≈ 2π
    @test dot(a2, b2) ≈ 2π
    @test abs(dot(a1, b2)) < 1e-10
    @test abs(dot(a2, b1)) < 1e-10
end
