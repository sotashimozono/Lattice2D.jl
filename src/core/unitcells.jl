"""
    Square <: AbstractTopology{2}

Standard square lattice: one sublattice, orthonormal primitive
vectors, four nearest neighbours per site.
"""
struct Square <: AbstractTopology{2} end

function get_unit_cell(::Type{Square})
    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]
    conns = [Connection(1, 1, 1, 0, 1), Connection(1, 1, 0, 1, 2)]
    # Unit square with corners at (0,0), (1,0), (1,1), (0,1), walked CCW.
    plaqs = [PlaquetteRule([(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)], :square)]
    return UnitCell{2,Float64}([a1, a2], [[0.0, 0.0]], conns, plaqs)
end

"""
    Triangular <: AbstractTopology{2}

Triangular lattice: one sublattice, 60° primitive vectors, six
nearest neighbours per site.
"""
struct Triangular <: AbstractTopology{2} end

function get_unit_cell(::Type{Triangular})
    a1 = [1.0, 0.0]
    a2 = [0.5, sqrt(3) / 2]
    conns = [
        Connection(1, 1, 1, 0, 1),    # right (a1)
        Connection(1, 1, 0, 1, 2),    # up (a2)
        Connection(1, 1, -1, 1, 3),   # upper-left (a2 − a1)
    ]
    # Two triangles per unit cell:
    #   up-triangle   : (0,0), (1,0), (0,1)       — points up-right
    #   down-triangle : (1,0), (1,1), (0,1)       — points down-left
    plaqs = [
        PlaquetteRule([(1, 0, 0), (1, 1, 0), (1, 0, 1)], :up_triangle),
        PlaquetteRule([(1, 1, 0), (1, 1, 1), (1, 0, 1)], :down_triangle),
    ]
    return UnitCell{2,Float64}([a1, a2], [[0.0, 0.0]], conns, plaqs)
end

"""
    Honeycomb <: AbstractTopology{2}

Honeycomb lattice: two-sublattice (A/B) bipartite triangular
Bravais lattice, three nearest neighbours per site.
"""
struct Honeycomb <: AbstractTopology{2} end

function get_unit_cell(::Type{Honeycomb})
    a1 = [sqrt(3), 0.0]
    a2 = [sqrt(3) / 2, 1.5]
    d_A = [0.0, 0.0]
    d_B = [0.0, 1.0]
    conns = [
        Connection(1, 2, 0, 0, 1),     # A → B (same cell: up)
        Connection(1, 2, 1, -1, 2),    # A → B (upper-left cell)
        Connection(1, 2, 0, -1, 3),    # A → B (upper cell)
    ]
    # Hexagon anchored at A(0,0), walking CCW. The 6 corners are 3 A
    # and 3 B sites in the cells (0,0), (1,-1), (1,0), (1,0), (0,1),
    # (0,0). Verified geometrically: every edge has length 1 and the
    # centroid lands at (sqrt(3)/2, 0.5).
    plaqs = [
        PlaquetteRule(
            [
                (1, 0, 0),     # A(0,0)
                (2, 1, -1),    # B(1,-1)
                (1, 1, 0),     # A(1,0)
                (2, 1, 0),     # B(1,0)
                (1, 0, 1),     # A(0,1)
                (2, 0, 0),
            ],    # B(0,0)
            :hexagon,
        ),
    ]
    return UnitCell{2,Float64}([a1, a2], [d_A, d_B], conns, plaqs)
end

"""
    Kagome <: AbstractTopology{2}

Kagome lattice: three-sublattice (A/B/C) structure on a triangular
Bravais lattice, four nearest neighbours per site.
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
        Connection(1, 2, 0, 0, 1),     # A–B
        Connection(1, 3, 0, 0, 1),     # A–C
        Connection(2, 3, 0, 0, 1),     # B–C
        # connections to neighbouring unit cells
        Connection(2, 1, 1, 0, 1),     # B → next A (right)
        Connection(3, 1, 0, 1, 1),     # C → next A (up)
        Connection(2, 3, 1, -1, 1),    # B → next C (down-right)
    ]
    # The up-triangle is the intra-cell triangle A-B-C. The
    # down-triangle lives between three neighbouring cells: B(0,0),
    # A(1,0), C(1,-1). The hexagon surrounds the hole between the
    # anchor cell and its (+x)/(+y) neighbours — a CCW walk through
    # B(0,0), A(1,0), C(1,0), B(0,1), A(0,1), C(0,0). Verified
    # numerically: all 6 corners lie at distance 0.5 from the hex
    # centre and every edge has unit length.
    plaqs = [
        PlaquetteRule([(1, 0, 0), (2, 0, 0), (3, 0, 0)], :up_triangle),
        PlaquetteRule([(2, 0, 0), (1, 1, 0), (3, 1, -1)], :down_triangle),
        PlaquetteRule(
            [(2, 0, 0), (1, 1, 0), (3, 1, 0), (2, 0, 1), (1, 0, 1), (3, 0, 0)], :hexagon
        ),
    ]
    return UnitCell{2,Float64}([a1, a2], [d_A, d_B, d_C], conns, plaqs)
end

"""
    Lieb <: AbstractTopology{2}

Lieb (line-centred square) lattice: three-sublattice structure
on a square Bravais lattice. Used for flat-band physics.
"""
struct Lieb <: AbstractTopology{2} end

function get_unit_cell(::Type{Lieb})
    a1 = [2.0, 0.0]
    a2 = [0.0, 2.0]

    # A: corner, B: right edge, C: top edge.
    d_A = [0.0, 0.0]
    d_B = [1.0, 0.0]
    d_C = [0.0, 1.0]

    conns = [
        Connection(1, 2, 0, 0, 1),     # A → B (intra-cell)
        Connection(1, 3, 0, 0, 2),     # A → C (intra-cell)
        Connection(2, 1, 1, 0, 1),     # B → next A (right)
        Connection(3, 1, 0, 1, 2),     # C → next A (up)
    ]
    # One plaquette per cell: the 2×2 square with 4 corner A's and 4
    # edge sites (2 B's bottom/top, 2 C's left/right). CCW walk:
    #   A(0,0) → B(0,0) → A(1,0) → C(1,0) → A(1,1) → B(0,1) → A(0,1) → C(0,0)
    plaqs = [PlaquetteRule(
        [
            (1, 0, 0),   # A(0,0) — bottom-left corner
            (2, 0, 0),   # B(0,0) — bottom edge
            (1, 1, 0),   # A(1,0) — bottom-right corner
            (3, 1, 0),   # C(1,0) — right edge
            (1, 1, 1),   # A(1,1) — top-right corner
            (2, 0, 1),   # B(0,1) — top edge
            (1, 0, 1),   # A(0,1) — top-left corner
            (3, 0, 0),   # C(0,0) — left edge
        ],
        :square,
    )]
    return UnitCell{2,Float64}([a1, a2], [d_A, d_B, d_C], conns, plaqs)
end

"""
    ShastrySutherland <: AbstractTopology{2}

Shastry–Sutherland lattice: square Bravais lattice with four
sites per unit cell, both nearest-neighbour (square) bonds and
dimer (J′) bonds.
"""
struct ShastrySutherland <: AbstractTopology{2} end

function get_unit_cell(::Type{ShastrySutherland})
    # Square Bravais lattice with a 2×2 four-site unit cell.
    a1 = [2.0, 0.0]
    a2 = [0.0, 2.0]

    d_1 = [0.0, 0.0]
    d_2 = [1.0, 0.0]
    d_3 = [0.0, 1.0]
    d_4 = [1.0, 1.0]

    conns = [
        # --- Nearest neighbour (square) bonds ---
        Connection(1, 2, 0, 0, 1),     # 1–2 (right)
        Connection(3, 4, 0, 0, 1),     # 3–4 (right)
        Connection(1, 3, 0, 0, 1),     # 1–3 (up)
        Connection(2, 4, 0, 0, 1),     # 2–4 (up)

        # Inter-cell nearest-neighbour bonds
        Connection(2, 1, 1, 0, 1),     # 2 → 1' (next right)
        Connection(4, 3, 1, 0, 1),     # 4 → 3'
        Connection(3, 1, 0, 1, 1),     # 3 → 1' (next up)
        Connection(4, 2, 0, 1, 1),     # 4 → 2'

        # --- Dimer bonds (diagonal J′) ---
        Connection(1, 4, 0, 0, 2),     # 1–4 (same cell diagonal)
        Connection(2, 3, 1, -1, 2),    # 2–3 (down-right diagonal)
    ]
    # Four small unit-squares per cell (the underlying square lattice
    # has 2Lx × 2Ly sites on a 2×2 sublattice unit cell, so 4 small
    # squares per cell). Each small square is walked CCW; the `type`
    # tag distinguishes squares that contain an intra-cell J' dimer
    # from those that don't.
    #   P1 (bottom-left, contains the 1-4 dimer):
    #     (0,0)-(1,0)-(1,1)-(0,1)  = sub 1,2,4,3
    #   P2 (bottom-right, dimer-free):
    #     (1,0)-(2,0)-(2,1)-(1,1)  = sub 2@(0,0), 1@(1,0), 3@(1,0), 4@(0,0)
    #   P3 (top-left, dimer-free):
    #     (0,1)-(1,1)-(1,2)-(0,2)  = sub 3@(0,0), 4@(0,0), 2@(0,1), 1@(0,1)
    #   P4 (top-right, contains the 2-3 dimer under cell shift):
    #     (1,1)-(2,1)-(2,2)-(1,2)  = sub 4@(0,0), 3@(1,0), 1@(1,1), 2@(0,1)
    plaqs = [
        PlaquetteRule(
            [(1, 0, 0), (2, 0, 0), (4, 0, 0), (3, 0, 0)], :dimer_square
        ),
        PlaquetteRule(
            [(2, 0, 0), (1, 1, 0), (3, 1, 0), (4, 0, 0)], :square
        ),
        PlaquetteRule(
            [(3, 0, 0), (4, 0, 0), (2, 0, 1), (1, 0, 1)], :square
        ),
        PlaquetteRule(
            [(4, 0, 0), (3, 1, 0), (1, 1, 1), (2, 0, 1)], :dimer_square
        ),
    ]
    return UnitCell{2,Float64}([a1, a2], [d_1, d_2, d_3, d_4], conns, plaqs)
end

"""
    Dice <: AbstractTopology{2}

Dice (T3) lattice: bipartite triangular-based structure with a
single 6-coordinated hub site and two 3-coordinated rim sites per
unit cell.
"""
struct Dice <: AbstractTopology{2} end

function get_unit_cell(::Type{Dice})
    # Triangular basis.
    a1 = [1.0, 0.0]
    a2 = [0.5, sqrt(3) / 2]

    d_1 = [0.0, 0.0]
    d_2 = (a1 .+ a2) ./ 3
    d_3 = (a1 .+ a2) .* (2 / 3)

    conns = [
        Connection(1, 2, 0, 0, 1),
        Connection(2, 1, 1, 0, 1),
        Connection(2, 1, 0, 1, 1),
        Connection(3, 1, 1, 1, 1),
        Connection(3, 1, 1, 0, 1),
        Connection(3, 1, 0, 1, 1),
    ]

    return UnitCell{2,Float64}([a1, a2], [d_1, d_2, d_3], conns)
end

"""
    UnionJack <: AbstractTopology{2}

Union Jack (centred square) lattice: square primitive cell with
two sublattices — a corner site and a body-centred site — and
eight-coordinated corner sites, four-coordinated body sites.
"""
struct UnionJack <: AbstractTopology{2} end

function get_unit_cell(::Type{UnionJack})
    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]

    d_A = [0.0, 0.0]   # corner
    d_B = [0.5, 0.5]   # body centre

    conns = [
        # Square-lattice bonds on A
        Connection(1, 1, 1, 0, 1),     # right
        Connection(1, 1, 0, 1, 1),     # up

        # Corner-to-centre (B) bonds
        Connection(1, 2, 0, 0, 1),     # A → B (intra-cell)
        Connection(1, 2, -1, 0, 1),    # A → B (left)
        Connection(1, 2, 0, -1, 1),    # A → B (down)
        Connection(1, 2, -1, -1, 1),   # A → B (down-left)
    ]
    # Four small triangles per cell, each sharing the body site B.
    # The four corner A's are walked in CCW order around B(0,0):
    #   T1 (south): A(0,0), A(1,0), B(0,0)
    #   T2 (east) : A(1,0), A(1,1), B(0,0)
    #   T3 (north): A(1,1), A(0,1), B(0,0)
    #   T4 (west) : A(0,1), A(0,0), B(0,0)
    plaqs = [
        PlaquetteRule([(1, 0, 0), (1, 1, 0), (2, 0, 0)], :triangle_south),
        PlaquetteRule([(1, 1, 0), (1, 1, 1), (2, 0, 0)], :triangle_east),
        PlaquetteRule([(1, 1, 1), (1, 0, 1), (2, 0, 0)], :triangle_north),
        PlaquetteRule([(1, 0, 1), (1, 0, 0), (2, 0, 0)], :triangle_west),
    ]
    return UnitCell{2,Float64}([a1, a2], [d_A, d_B], conns, plaqs)
end

"""Tuple listing every topology shipped by Lattice2D."""
const AVAILABLE_LATTICES = (
    Square, Triangular, Honeycomb, Kagome, Lieb, ShastrySutherland, Dice, UnionJack
)
