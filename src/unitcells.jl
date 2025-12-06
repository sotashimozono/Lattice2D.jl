"""
    get_unit_cell(::Type{T}) where T <: AbstractTopology
this function returns the UnitCell associated with the given Topology type.
If the Topology type is not recognized, it throws an error.
"""
function get_unit_cell(::Type{T}) where T<:AbstractTopology
    error("UnitCell not defined for $T")
end
export get_unit_cell
"""
    Square <: AbstractTopology{2}
struct which represents Square lattice
"""
struct Square <: AbstractTopology{2} end
function get_unit_cell(::Type{Square})
    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]
    conns = [
        Connection(1, 1, 1, 0, 1),
        Connection(1, 1, 0, 1, 2)
    ]
    return UnitCell{2,Float64}([a1, a2], [[0.0, 0.0]], conns)
end

"""
    Triangular <: AbstractTopology{2}
struct which represents Triangular lattice
"""
struct Triangular <: AbstractTopology{2} end
function get_unit_cell(::Type{Triangular})
    a1 = [1.0, 0.0]
    a2 = [0.5, sqrt(3) / 2]
    conns = [
        Connection(1, 1, 1, 0, 1), # right (a1)
        Connection(1, 1, 0, 1, 2), # up (a2)
        Connection(1, 1, -1, 1, 3) # upper left (a2-a1)
    ]
    return UnitCell{2,Float64}([a1, a2], [[0.0, 0.0]], conns)
end
"""
    Honeycomb <: AbstractTopology{2}
struct which represents Honeycomb lattice
"""
struct Honeycomb <: AbstractTopology{2} end
function get_unit_cell(::Type{Honeycomb})
    a1 = [sqrt(3), 0.0]
    a2 = [sqrt(3) / 2, 1.5]
    d_A = [0.0, 0.0]
    d_B = [0.0, 1.0]
    conns = [
        Connection(1, 2, 0, 0, 1),  # A -> B (same cell: up)
        Connection(1, 2, 1, -1, 2), # A -> B (upper left cell)
        Connection(1, 2, 0, -1, 3)   # A -> B (upper cell)
    ]
    return UnitCell{2,Float64}([a1, a2], [d_A, d_B], conns)
end
"""
    Kagome <: AbstractTopology{2}
struct which represents Kagome lattice
"""
struct Kagome <: AbstractTopology{2} end
function get_unit_cell(::Type{Kagome})
    a1 = [1.0, 0.0]
    a2 = [0.5, sqrt(3) / 2]
    d_A = [0.0, 0.0]
    d_B = 0.5 * a1
    d_C = 0.5 * a2
    conns = [
        # triangle connections within the unit cell
        Connection(1, 2, 0, 0, 1), # A-B
        Connection(1, 3, 0, 0, 1), # A-C
        Connection(2, 3, 0, 0, 1), # B-C
        # connections to neighboring unit cells
        Connection(2, 1, 1, 0, 1), # B -> Next A (right)
        Connection(3, 1, 0, 1, 1), # C -> Next A (up)
        Connection(2, 3, 1, -1, 1) # B -> Next C (down-right: Kagome specific)
    ]
    return UnitCell{2,Float64}([a1, a2], [d_A, d_B, d_C], conns)
end
"""
    Lieb <: AbstractTopology{2}
struct which represents Lieb lattice
"""
struct Lieb <: AbstractTopology{2} end
function get_unit_cell(::Type{Lieb})
    a1 = [2.0, 0.0]
    a2 = [0.0, 2.0]

    # A:Corner, B:Right, C:Up
    d_A = [0.0, 0.0]
    d_B = [1.0, 0.0]
    d_C = [0.0, 1.0]

    conns = [
        Connection(1, 2, 0, 0, 1), # A -> B (Intra)
        Connection(1, 3, 0, 0, 2), # A -> C (Intra)
        Connection(2, 1, 1, 0, 1), # B -> Next A (Right)
        Connection(3, 1, 0, 1, 2)  # C -> Next A (Up)
    ]
    return UnitCell{2,Float64}([a1, a2], [d_A, d_B, d_C], conns)
end

"""
    ShastrySutherland <: AbstractTopology{2}
struct which represents Shastry-Sutherland lattice
"""
struct ShastrySutherland <: AbstractTopology{2} end
function get_unit_cell(::Type{ShastrySutherland})
    # 基本は正方格子だが、ユニットセル内に4サイトある
    a1 = [2.0, 0.0]
    a2 = [0.0, 2.0]
    
    # 4サイト (2x2の正方格子の配置)
    d_1 = [0.0, 0.0]
    d_2 = [1.0, 0.0]
    d_3 = [0.0, 1.0]
    d_4 = [1.0, 1.0]

    conns = [
        # --- Nearest Neighbors (Square bonds) ---
        Connection(1, 2, 0, 0, 1), # 1-2 (Right)
        Connection(3, 4, 0, 0, 1), # 3-4 (Right)
        Connection(1, 3, 0, 0, 1), # 1-3 (Up)
        Connection(2, 4, 0, 0, 1), # 2-4 (Up)
        
        # Inter-cell connections for Square
        Connection(2, 1, 1, 0, 1), # 2->1' (Next Right)
        Connection(4, 3, 1, 0, 1), # 4->3'
        Connection(3, 1, 0, 1, 1), # 3->1' (Next Up)
        Connection(4, 2, 0, 1, 1), # 4->2'

        # --- Dimer Bonds (Diagonal) J' ---
        Connection(1, 4, 0, 0, 2), # 1-4 (Cell内 対角)
        Connection(2, 3, 1, -1, 2) # 2-3 (右下のセルとの対角)
    ]
    return UnitCell{2,Float64}([a1, a2], [d_1, d_2, d_3, d_4], conns)
end

const AVAILABLE_LATTICES = (Square, Triangular, Honeycomb, Kagome, Lieb, ShastrySutherland)
export AVAILABLE_LATTICES
export Square, Triangular, Honeycomb, Kagome, Lieb, ShastrySutherland