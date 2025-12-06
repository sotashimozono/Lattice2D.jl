@testset "check_bipartite_bfs" begin
    # Case 1: Simple Line 1-2-3 (Bipartite)
    # 1 -- 2 -- 3
    # Color: A - B - A (Valid)
    neighbors_line = [
        [2],    # 1 connects to 2
        [1, 3], # 2 connects to 1, 3
        [2]     # 3 connects to 2
    ]
    @test check_bipartite_bfs(3, neighbors_line) == true

    # Case 2: Triangle 1-2-3-1 (Not Bipartite)
    # 1 -- 2
    #  \  /
    #   3
    neighbors_tri = [
        [2, 3], # 1
        [1, 3], # 2
        [1, 2]  # 3
    ]
    @test check_bipartite_bfs(3, neighbors_tri) == false

    # Case 3: Disconnected Graph (Two separate lines)
    # 1-2, 3-4 (Bipartite)
    neighbors_dis = [
        [2], [1], # Component 1
        [4], [3]  # Component 2
    ]
    @test check_bipartite_bfs(4, neighbors_dis) == true
end
