"""
    AbstractLattice{D}
格子の抽象型。Dは空間次元。
ising model などの格子モデルはこの型を継承して実装する。
"""
abstract type AbstractLattice{D} end
"""
    AbstractLatticeConnection
格子の接続（辺、ボンド）を表す抽象型。
"""
abstract type AbstractLatticeConnection end

"""
    Bond
格子の辺を表す型。
- `src::Int`: 辺の始点サイトのインデックス
- `dst::Int`: 辺の終点サイトのインデックス
- `type::Int`: 辺の種類（異なる結合定数を持つ辺を区別するため）
- `vector::Vector{Float64}`: 辺のベクトル表現
"""
struct Bond <: AbstractLatticeConnection
    src::Int
    dst::Int
    type::Int
    vector::Vector{Float64}
end
export Bond

"""
    Connection
単位胞内または単位胞間の接続ルール。
- `src_sub`: 始点のサブ格子インデックス (1, 2, ...)
- `dst_sub`: 終点のサブ格子インデックス
- `dx`, `dy`: 終点がどの相対セルにあるか (0,0 なら同じ単位胞内)
- `type`: 結合の種類
"""
struct Connection <: AbstractLatticeConnection
    src_sub::Int
    dst_sub::Int
    dx::Int
    dy::Int
    type::Int
end
export Connection

"""
    UnitCell{D, T}
格子の幾何学的定義データ。基本的に、この情報をもとに格子を構築する。
`get_unit_cell(::Type{T})` 関数で各トポロジーに対応する単位胞データを取得する。
"""
struct UnitCell{D,T} <: AbstractLattice{D}
    basis::Vector{Vector{T}}
    sublattice_positions::Vector{Vector{T}}
    connections::Vector{Connection}
end
export UnitCell
"""
    AbstractTopology
格子のトポロジーの抽象型。
- Square: 正方格子
- Triangular: 三角格子
- Honeycomb: ハニカム格子
- Kagome: カゴメ格子
"""
abstract type AbstractTopology{D} <: AbstractLattice{D} end
include("lattices.jl")

"""
Lattice2D{Topology<:AbstractTopology, T, B<:AbstractBoundaryCondition}
2次元格子を表す型。
- `Lx::Int`: x方向の格子サイズ
- `Ly::Int`: y方向の格子サイズ
- `N::Int`: サイト総数
- `positions::Vector{Vector{T}}`: 各サイトの位置ベクトル
- `nearest_neighbors::Vector{Vector{Int}}`: 各サイトの最近接サイトのインデックスリスト
- `bonds::Vector{Bond}`: 格子の辺のリスト
- `basis_vectors::Vector{Vector{T}}`: 格子の基底ベクトル
- `reciprocal_vectors::Union{Vector{Vector{T}}, Nothing}`: 逆格子ベクトル
- `sublattice_ids::Vector{Int}`: 各サイトのサブ格子ID
- `is_bipartite::Bool`: 格子が二部グラフかどうか
- `site_map::Union{Matrix{Int}, Nothing}`: 格子上のサイトインデックスのマッピング
- `translation_x::Vector{Int}`: x方向の平行移動ベクトル
- `translation_y::Vector{Int}`: y方向の平行移動ベクトル
- `boundary::B`: 境界条件
"""
struct Lattice2D{Topology<:AbstractTopology,T,B<:AbstractBoundaryCondition} <: AbstractLattice{2}
    Lx::Int
    Ly::Int
    N::Int
    positions::Vector{Vector{T}}
    # graph representation
    nearest_neighbors::Vector{Vector{Int}}
    bonds::Vector{Bond}
    # topology information
    basis_vectors::Vector{Vector{T}}
    reciprocal_vectors::Union{Vector{Vector{T}},Nothing}
    sublattice_ids::Vector{Int}
    is_bipartite::Bool
    site_map::Union{Matrix{Int},Nothing}
    boundary::B
end
export Lattice2D

