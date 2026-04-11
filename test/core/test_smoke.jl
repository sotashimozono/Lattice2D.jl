@testset "smoke test: Square lattice builds and passes basic checks" begin
    lat = build_lattice(Square, 4, 4)

    # Basic counts
    @test num_sites(lat) == 16
    @test lat.Lx == 4 && lat.Ly == 4
    @test num_sublattices(lat) == 1

    # Every interior site has 4 neighbours under PBC
    degrees = [length(neighbors(lat, i)) for i in 1:num_sites(lat)]
    @test all(d == 4 for d in degrees)

    # Site positions are unit-spaced
    @test position(lat, 1) == SVector(0.0, 0.0)
    @test position(lat, 2) == SVector(1.0, 0.0)

    # The 4×4 square is bipartite
    @test is_bipartite(lat) == true
    @test periodicity(lat) isa Periodic
    @test reciprocal_support(lat) isa HasReciprocal

    # Size trait
    @test size_trait(lat) isa FiniteSize{2}
    @test size_trait(lat).dims == (4, 4)

    # Reciprocal lattice: orthogonality a_i · b_j = 2π δ_ij
    ml = reciprocal_lattice(lat)
    @test num_sites(ml) == 16
    A = basis_vectors(lat)
    B = ml.reciprocal_basis
    @test dot(A[:, 1], B[:, 1]) ≈ 2π atol = 1e-10
    @test dot(A[:, 2], B[:, 2]) ≈ 2π atol = 1e-10
    @test abs(dot(A[:, 1], B[:, 2])) < 1e-10

    # Default layout is a uniform Ising layout
    @test site_layout(lat) isa UniformLayout
    @test site_type(lat, 1) isa IsingSite
end
