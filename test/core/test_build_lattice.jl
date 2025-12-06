struct TestSimpleSquare <: Lattices.AbstractTopology{2} end

function Lattices.get_unit_cell(::Type{TestSimpleSquare})
    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]
    conns = [
        Connection(1, 1, 1, 0, 1), # 右 (+x)
        Connection(1, 1, 0, 1, 1)  # 上 (+y)
    ]
    return UnitCell{2,Float64}([a1, a2], [[0.0, 0.0]], conns)
end

@testset "Build Lattice Generic Logic" begin
    Lx, Ly = 5, 4 # あえて異なるサイズにする
    
    @testset "Position Generation" begin
        lat = build_lattice(TestSimpleSquare, Lx, Ly)
        for x in 1:Lx, y in 1:Ly
            idx = lat.site_map[x, y]
            pos = lat.positions[idx]
            expected = [(x - 1) * 1.0, (y - 1) * 1.0]
            @test pos ≈ expected atol=1e-10
        end
    end

    @testset "Boundary Conditions" begin
        lat_obc = build_lattice(TestSimpleSquare, Lx, Ly; boundary=OBC())
        
        # 右端 (Lx, y) は右 (Lx+1, y) に繋がってはいけない
        for y in 1:Ly
            idx_right_edge = lat_obc.site_map[Lx, y]
            neighbors = lat_obc.nearest_neighbors[idx_right_edge]
            # 隣接リストの中に「x=1 (左端)」のIDが含まれていないことを確認
            idx_left_edge = lat_obc.site_map[1, y]
            @test !(idx_left_edge in neighbors)
        end

        lat_pbc = build_lattice(TestSimpleSquare, Lx, Ly; boundary=PBC())
        
        # 右端 (Lx, y) は左端 (1, y) に繋がっていなければならない
        for y in 1:Ly
            idx_right_edge = lat_pbc.site_map[Lx, y]
            idx_left_edge = lat_pbc.site_map[1, y]
            
            neighbors = lat_pbc.nearest_neighbors[idx_right_edge]
            @test idx_left_edge in neighbors
        end
    end

    @testset "Graph Consistency" begin
        lat = build_lattice(TestSimpleSquare, Lx, Ly; boundary=PBC())
        
        for i in 1:lat.N
            for j in lat.nearest_neighbors[i]
                @test i in lat.nearest_neighbors[j]
            end
        end

        for b in lat.bonds
            @test b.dst in lat.nearest_neighbors[b.src]
            @test b.src in lat.nearest_neighbors[b.dst]
            
            @test norm(b.vector) ≈ 1.0 atol=1e-10
        end
    end
    
    @testset "Bipartite Logic Check" begin
        lat_even = build_lattice(TestSimpleSquare, 4, 4; boundary=PBC())
        @test lat_even.is_bipartite == true
        
        lat_odd = build_lattice(TestSimpleSquare, 3, 4; boundary=PBC())
        @test lat_odd.is_bipartite == false
        
        lat_odd_obc = build_lattice(TestSimpleSquare, 3, 4; boundary=OBC())
        @test lat_odd_obc.is_bipartite == true
    end
end