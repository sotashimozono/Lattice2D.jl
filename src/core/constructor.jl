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
function build_lattice(
    Topology::Type{<:AbstractTopology},
    Lx::Int,
    Ly::Int;
    boundary::AbstractBoundaryCondition=PBC(),
    index_method::AbstractIndexing=RowMajorIndexing(),
)
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
            d_vec =
                (conn.dx * a1 + conn.dy * a2) + (
                    uc.sublattice_positions[conn.dst_sub] -
                    uc.sublattice_positions[conn.src_sub]
                )
            push!(nn_table[src], dst)
            push!(nn_table[dst], src)
            push!(bonds, Bond(src, dst, conn.type, d_vec))
        end
    end
    # --- 3. Finalize ---
    recip_vecs = calc_reciprocal_vectors(uc.basis)
    is_bipartite = check_bipartite_bfs(N, nn_table)
    return Lattice{Topology,Float64,typeof(boundary),typeof(index_method)}(
        Lx,
        Ly,
        N,
        positions,
        nn_table,
        bonds,
        uc.basis,
        recip_vecs,
        sublattice_ids,
        is_bipartite,
        site_map,
        boundary,
        index_method,
    )
end
export build_lattice

function Lattice(Topology::Type{<:AbstractTopology}, Lx::Int, Ly::Int; kwargs...)
    build_lattice(Topology::Type{<:AbstractTopology}, Lx::Int, Ly::Int; kwargs...)
end
"""
    get_site_index(lat::Lattice, x::Int, y::Int, s::Int=1)

Get the linear site index from lattice coordinates (x, y) and sublattice index s.
For single-sublattice lattices, s defaults to 1.

# Arguments
- `lat::Lattice`: The lattice structure
- `x::Int`: x-coordinate (1 to Lx)
- `y::Int`: y-coordinate (1 to Ly)
- `s::Int`: sublattice index (1 to number of sublattices), defaults to 1

# Returns
- `Int`: The linear site index

# Examples
```julia
lat = build_lattice(Square, 4, 4)
idx = get_site_index(lat, 2, 3)  # Get index at (2,3)

lat = build_lattice(Honeycomb, 4, 4)
idx_A = get_site_index(lat, 1, 1, 1)  # Get sublattice A at (1,1)
idx_B = get_site_index(lat, 1, 1, 2)  # Get sublattice B at (1,1)
```
"""
function get_site_index(lat::Lattice, x::Int, y::Int, s::Int=1)
    # Get number of sublattices from the lattice
    # We can infer this from the site_map structure
    n_sub = length(lat.sublattice_ids) ÷ (lat.Lx * lat.Ly)

    # Validate inputs
    if x < 1 || x > lat.Lx
        throw(ArgumentError("x must be in range [1, $(lat.Lx)], got $x"))
    end
    if y < 1 || y > lat.Ly
        throw(ArgumentError("y must be in range [1, $(lat.Ly)], got $y"))
    end
    if s < 1 || s > n_sub
        throw(ArgumentError("s must be in range [1, $n_sub], got $s"))
    end

    return _coord_to_index(lat.index_method, x, y, s, lat.Lx, lat.Ly, n_sub)
end
export get_site_index

"""
    get_position(lat::Lattice, idx::Int)

Get the spatial position (coordinates) of a site given its linear index.

# Arguments
- `lat::Lattice`: The lattice structure
- `idx::Int`: The linear site index

# Returns
- `Vector{Float64}`: The position vector [x, y] in real space

# Examples
```julia
lat = build_lattice(Square, 4, 4)
pos = get_position(lat, 5)  # Get position of site 5
```
"""
function get_position(lat::Lattice, idx::Int)
    if idx < 1 || idx > lat.N
        throw(ArgumentError("idx must be in range [1, $(lat.N)], got $idx"))
    end
    return lat.positions[idx]
end
export get_position

"""
    get_coordinates(lat::Lattice, idx::Int)

Get the lattice coordinates (x, y, sublattice) from a linear site index.

# Arguments
- `lat::Lattice`: The lattice structure
- `idx::Int`: The linear site index

# Returns
- `Tuple{Int, Int, Int}`: A tuple (x, y, s) where x, y are unit cell coordinates and s is sublattice index

# Examples
```julia
lat = build_lattice(Honeycomb, 4, 4)
x, y, s = get_coordinates(lat, 10)
```
"""
function get_coordinates(lat::Lattice, idx::Int)
    if idx < 1 || idx > lat.N
        throw(ArgumentError("idx must be in range [1, $(lat.N)], got $idx"))
    end

    n_sub = length(lat.sublattice_ids) ÷ (lat.Lx * lat.Ly)
    s = lat.sublattice_ids[idx]

    # For RowMajorIndexing, we can reverse the calculation
    if isa(lat.index_method, RowMajorIndexing)
        # idx = cell_index_0based * n_sub + s
        # where cell_index_0based = (y - 1) * Lx + (x - 1)
        cell_index_0based = (idx - s) ÷ n_sub
        y = cell_index_0based ÷ lat.Lx + 1
        x = cell_index_0based % lat.Lx + 1
        return (x, y, s)
    elseif isa(lat.index_method, ColMajorIndexing)
        # idx = cell_index_0based * n_sub + s
        # where cell_index_0based = (x - 1) * Ly + (y - 1)
        cell_index_0based = (idx - s) ÷ n_sub
        x = cell_index_0based ÷ lat.Ly + 1
        y = cell_index_0based % lat.Ly + 1
        return (x, y, s)
    elseif isa(lat.index_method, SnakeIndexing)
        # For SnakeIndexing, we need to reverse the snake pattern
        cell_index_0based = (idx - s) ÷ n_sub
        y = cell_index_0based ÷ lat.Lx + 1
        x_map = cell_index_0based % lat.Lx + 1
        # Reverse the snake pattern
        x = isodd(y) ? x_map : lat.Lx + 1 - x_map
        return (x, y, s)
    else
        error("Unknown indexing method: $(typeof(lat.index_method))")
    end
end
export get_coordinates
