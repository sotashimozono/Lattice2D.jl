@testset "_coord_to_index" begin
    # Logic: ((x - 1) + (y - 1) * Lx) * n_sub + s
    # Test params: Lx=4, n_sub=2
    Lx = 4
    n_sub = 2
    
    # (1, 1, 1) -> 1
    @test Lattices._coord_to_index(1, 1, 1, Lx, n_sub) == 1
    # (1, 1, 2) -> 2
    @test Lattices._coord_to_index(1, 1, 2, Lx, n_sub) == 2
    
    # (2, 1, 1) -> ((1) + 0)*2 + 1 = 3
    @test Lattices._coord_to_index(2, 1, 1, Lx, n_sub) == 3
    
    # (1, 2, 1) -> (0 + (1)*4)*2 + 1 = 9
    @test Lattices._coord_to_index(1, 2, 1, Lx, n_sub) == 9
    
    # Max index: (Lx, Ly, n_sub)
    Ly = 3
    # ((4-1) + (3-1)*4)*2 + 2 = (3 + 8)*2 + 2 = 24
    @test Lattices._coord_to_index(Lx, Ly, n_sub, Lx, n_sub) == Lx * Ly * n_sub
end