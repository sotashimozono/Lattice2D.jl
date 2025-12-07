@testset "_coord_to_index" begin
    @testset "RowMajorIndexing Logic" begin
        Lx = 4
        Ly = 3
        n_sub = 2
        RMI = RowMajorIndexing()
        # Base function signature:
        # Lattice2D._coord_to_index(::AbstractIndexing, x, y, s, Lx, Ly, n_sub)
        # === 1. セル (1, 1) のテスト ===
        # (x, y, s) = (1, 1, 1) -> [ (0*4 + 0) ] * 2 + 1 = 1
        @test Lattice2D._coord_to_index(RMI, 1, 1, 1, Lx, Ly, n_sub) == 1
        # (x, y, s) = (1, 1, 2) -> [ 0 ] * 2 + 2 = 2
        @test Lattice2D._coord_to_index(RMI, 1, 1, 2, Lx, Ly, n_sub) == 2
        
        # === 2. 同じ行 (y=1) のテスト ===
        # (x, y, s) = (2, 1, 1) -> [ (0*4 + 1) ] * 2 + 1 = 3
        @test Lattice2D._coord_to_index(RMI, 2, 1, 1, Lx, Ly, n_sub) == 3
        # (x, y, s) = (4, 1, 2) -> [ (0*4 + 3) ] * 2 + 2 = 8
        @test Lattice2D._coord_to_index(RMI, 4, 1, 2, Lx, Ly, n_sub) == 8 # y=1 の最大インデックス
        
        # === 3. 2行目 (y=2) のテスト (インデックスがジャンプする箇所) ===
        # (x, y, s) = (1, 2, 1) -> [ (1*4 + 0) ] * 2 + 1 = 9
        @test Lattice2D._coord_to_index(RMI, 1, 2, 1, Lx, Ly, n_sub) == 9 
        
        # (x, y, s) = (4, 2, 2) -> [ (1*4 + 3) ] * 2 + 2 = 16
        @test Lattice2D._coord_to_index(RMI, 4, 2, 2, Lx, Ly, n_sub) == 16 
        
        # === 4. 最終行 (y=3) の最大インデックス ===
        # Cell Index max = (3-1)*4 + (4-1) = 11
        # Global Index max = 11 * 2 + 2 = 24
        @test Lattice2D._coord_to_index(RMI, Lx, Ly, n_sub, Lx, Ly, n_sub) == Lx * Ly * n_sub
        
        # === 5. ColMajorとの比較 (Lx=4, Ly=3 の場合、(2, 1)と(1, 2)の順序が異なる) ===
        # RowMajor: (2, 1, 1) -> 3
        # ColMajor: (1, 2, 1) -> 9
    end
    @testset "SnakeIndexing Logic" begin
        Lx = 4
        Ly = 3
        n_sub = 2
        SI = SnakeIndexing()
        # === 1. 奇数行 (y=1): 昇順 ===
        # (1, 1, 1) -> 1 (RowMajorと同じ)
        @test Lattice2D._coord_to_index(SI, 1, 1, 1, Lx, Ly, n_sub) == 1
        # (4, 1, 2) -> 8 (RowMajorと同じ)
        @test Lattice2D._coord_to_index(SI, 4, 1, 2, Lx, Ly, n_sub) == 8
        
        # === 2. 偶数行 (y=2): 降順 (4->1) ===
        # y=2 の Cell Index ベースは (2-1)*4 = 4
        # (x_map: 4 -> x_eff=4, base_x=3) 
        # (x, y, s) = (4, 2, 1) -> (4*2)+1 = 9 (通常RowMajorの(1, 2, 1)に相当)
        @test Lattice2D._coord_to_index(SI, 4, 2, 1, Lx, Ly, n_sub) == 9 
        
        # (x_map: 1 -> x_eff=1, base_x=0)
        # (x, y, s) = (1, 2, 2) -> (7*2)+2 = 16 (通常RowMajorの(4, 2, 2)に相当)
        @test Lattice2D._coord_to_index(SI, 1, 2, 2, Lx, Ly, n_sub) == 16 
        
        # === 3. 最終行 (y=3): 奇数行なので昇順 ===
        # (x, y, s) = (1, 3, 1) -> [ (2*4 + 0) ] * 2 + 1 = 17
        @test Lattice2D._coord_to_index(SI, 1, 3, 1, Lx, Ly, n_sub) == 17
    end
end