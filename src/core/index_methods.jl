abstract type AbstractIndexing end
struct RowMajorIndexing <: AbstractIndexing end
export RowMajorIndexing
struct ColMajorIndexing <: AbstractIndexing end
export ColMajorIndexing
struct SnakeIndexing <: AbstractIndexing end
export SnakeIndexing

function _coord_to_index(::RowMajorIndexing, x::Int, y::Int, s::Int, Lx::Int, Ly::Int, n_sub::Int)
    cell_index_0based = (y - 1) * Lx + (x - 1) 
    return cell_index_0based * n_sub + s
end
function _coord_to_index(::ColMajorIndexing, x::Int, y::Int, s::Int, Lx::Int, Ly::Int, n_sub::Int)
    cell_index_0based = (x - 1) * Ly + (y - 1) 
    return cell_index_0based * n_sub + s
end

function _coord_to_index(::SnakeIndexing, x::Int, y::Int, s::Int, Lx::Int, Ly::Int, n_sub::Int)
    # 1. 偶数行か奇数行かで x 座標の順序を反転
    x_map = isodd(y) ? x : Lx + 1 - x
    
    # 2. その y の行が持つ Cell Index のベース（y-1行目までのサイト数）
    y_offset = (y - 1) * Lx
    
    # 3. 0-based の Cell Index を計算: オフセット + マップされた x 座標
    cell_index_0based = y_offset + (x_map - 1)
    return cell_index_0based * n_sub + s
end