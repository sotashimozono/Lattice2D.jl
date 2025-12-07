# TODO : helper functions for Lattice
"""
    calc_reciprocal_vectors(basis)
Calculate reciprocal lattice vectors from the given basis vectors.
"""
function calc_reciprocal_vectors(basis)
    MatA = hcat(basis...)
    MatB = 2π * inv(MatA')
    return [MatB[:, i] for i in 1:size(MatB, 2)]
end
export calc_reciprocal_vectors
"""
    check_bipartite_bfs(N, neighbors)
check whether the given lattice is bipartite using BFS.
"""
function check_bipartite_bfs(N::Int, neighbors::Vector{Vector{Int}})
    colors = zeros(Int, N)
    # 0: unvisited, 1: colorA, -1: colorB
    # is_bipartite = true
    for i in 1:N
        if colors[i] == 0
            colors[i] = 1
            queue = [i]
            while !isempty(queue)
                u = popfirst!(queue)
                for v in neighbors[u]
                    if colors[v] == 0
                        colors[v] = -colors[u]
                        push!(queue, v)
                    elseif colors[v] == colors[u]
                        return false
                    end
                end
            end
        end
    end
    return true
end
export check_bipartite_bfs

"""
    _coord_to_index(lat::Lattice, x::Int, y::Int, s::Int)
this function converts lattice coordinates to a linear index based on the lattice's indexing method.
- ColMajorIndexing
- RowMajorIndexing
- SnakeIndexing
are available indexing methods.
"""
function _coord_to_index(lat::Lattice, x::Int, y::Int, s::Int)
    n_sub = length(lat.unit_cell.sublattice_positions)
    return _coord_to_index(lat.index_method, x, y, s, lat.Lx, lat.Ly, n_sub)
end

"""
    build_lattice(Topology::Type{<:AbstractTopology}, Lx::Int, Ly::Int; boundary::AbstractBoundaryCondition=PBC())
Construct a lattice with the specified topology, size, and boundary conditions.
this function is available if unitcell information is defined.
"""
function build_lattice(Topology::Type{<:AbstractTopology}, Lx::Int, Ly::Int; 
    boundary::AbstractBoundaryCondition=PBC(), index_method::AbstractIndexing=RowMajorIndexing())
    # --- 1. Initialize Lattice Parameters ---
    uc = get_unit_cell(Topology)
    n_sub = length(uc.sublattice_positions)
    N = Lx * Ly * n_sub
    # --- Pre-allocation ---
    positions = Vector{Vector{Float64}}(undef, N)
    site_map = Matrix{Int}(undef, Lx, Ly)
    sublattice_ids = Vector{Int}(undef, N)
    nn_table = [Int[] for _ in 1:N] # Inner vectors will grow, but outer is fixed

    # Bond数の見積もり: (UnitCell内のボンド数) * (セル数)
    est_bonds = length(uc.connections) * Lx * Ly
    bonds = Bond[]
    sizehint!(bonds, est_bonds)

    a1, a2 = uc.basis[1], uc.basis[2]
    idx = 1
    for y in 1:Ly
        for x in 1:Lx
            # site_map には「そのセルの最初のサブ格子(s=1)のID」を保存
            base_id = _coord_to_index(index_method, x, y, 1, Lx, Ly, n_sub)
            site_map[x, y] = base_id
            cell_origin = (x - 1) * a1 + (y - 1) * a2
            for s in 1:n_sub
                i = base_id + (s - 1)
                
                # positions, sublattice_ids の設定
                positions[i] = cell_origin + uc.sublattice_positions[s]
                sublattice_ids[i] = s
            end
        end
    end
    # --- 2. Bond Generation (Connection Rules) ---
    is_OBC = isa(boundary, OBC)
    for y in 1:Ly, x in 1:Lx
        # 現在のセルのベースID (s=1 の ID)
        # base_src = ( (x-1) + (y-1)*Lx ) * n_sub + 1
        # 計算済みの site_map を使うのと算術計算は、メモリアクセス vs ALU のトレードオフ。
        # ここでは可読性とキャッシュ効率のバランスで site_map を使う（シーケンシャルアクセスなので速い）
        base_src = site_map[x, y]
        for conn in uc.connections
            # 始点
            src = base_src + (conn.src_sub - 1)
            # 終点セル座標
            tx = x + conn.dx
            ty = y + conn.dy
            if is_OBC && (tx < 1 || tx > Lx || ty < 1 || ty > Ly)
                continue
            end
            # PBC Wrap (mod1を使わない最適化: 分岐予測が効くならifの方が速い場合も多いが、mod1は安全)
            # cx = mod1(tx, Lx); cy = mod1(ty, Ly)
            cx = (tx - 1) % Lx + 1
            cx = cx < 1 ? cx + Lx : cx
            cy = (ty - 1) % Ly + 1
            cy = cy < 1 ? cy + Ly : cy
            # 終点IDの計算 (site_mapを見ずに算術で出す)
            # dst_base = ((cx-1) + (cy-1)*Lx) * n_sub + 1
            # dst = dst_base + (conn.dst_sub - 1)
            dst = site_map[cx, cy] + (conn.dst_sub - 1)
            # Vector calculation
            d_vec = (conn.dx * a1 + conn.dy * a2) +
                    (uc.sublattice_positions[conn.dst_sub] - uc.sublattice_positions[conn.src_sub])
            push!(nn_table[src], dst)
            push!(nn_table[dst], src)
            push!(bonds, Bond(src, dst, conn.type, d_vec))
        end
    end
    # --- 3. Finalize ---
    recip_vecs = calc_reciprocal_vectors(uc.basis)
    is_bipartite = check_bipartite_bfs(N, nn_table)
    return Lattice{Topology,Float64,typeof(boundary), typeof(index_method)}(
        Lx, Ly, N, positions, nn_table, bonds,
        uc.basis, recip_vecs, sublattice_ids, is_bipartite,
        site_map, boundary, index_method
    )
end
export build_lattice
