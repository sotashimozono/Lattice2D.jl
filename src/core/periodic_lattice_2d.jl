"""
    PeriodicLattice2D{Topo, T, B, I, L}

Finite 2D periodic lattice built from a topology
([`get_unit_cell`](@ref)) and a
[`LatticeCore.LatticeBoundary`](@ref LatticeCore.LatticeBoundary).
Subtype of `LatticeCore.AbstractLattice{2, T}`; every
LatticeCore-side method (`num_sites`, `position`, `neighbors`,
`boundary`, `size_trait`, `site_layout`, `reciprocal_lattice`, ...)
is implemented below.

Type parameters

- `Topo<:AbstractTopology{2}` — singleton topology marker
- `T<:AbstractFloat` — numeric type for positions
- `B<:LatticeCore.LatticeBoundary` — composite boundary condition
- `I<:LatticeCore.AbstractIndexing` — linearisation strategy
- `L<:LatticeCore.AbstractSiteLayout` — site-type layout

Fields

- `Lx`, `Ly`, `N` — per-axis cell counts and total site count
- `positions` — cached real-space positions, indexed by site id
- `nearest_neighbors` — adjacency lists (one `Vector{Int}` per site)
- `bonds` — list of `LatticeCore.Bond{2, T}` with wrapped vectors
- `basis_matrix` — the `2 × 2` real-space basis (columns are
  primitive vectors)
- `sublattice_ids` — geometric sublattice id (1-based) of each site
- `is_bipartite_flag` — cached result of a BFS 2-colouring
- `site_map` — `Lx × Ly` matrix of the base site id of each cell
- `topology` — topology singleton
- `boundary` — composite boundary
- `indexing` — `RowMajor`, `ColMajor`, or `Snake`
- `layout` — site-type layout

Prefer constructing via [`build_lattice`](@ref) rather than the
inner constructor.
"""
struct PeriodicLattice2D{
    Topo<:AbstractTopology{2},
    T<:AbstractFloat,
    B<:LatticeBoundary,
    I<:AbstractIndexing,
    L<:AbstractSiteLayout,
} <: AbstractLattice{2,T}
    Lx::Int
    Ly::Int
    N::Int
    positions::Vector{SVector{2,T}}
    nearest_neighbors::Vector{Vector{Int}}
    bonds::Vector{Bond{2,T}}
    basis_matrix::SMatrix{2,2,T,4}
    sublattice_ids::Vector{Int}
    is_bipartite_flag::Bool
    site_map::Matrix{Int}
    topology::Topo
    boundary::B
    indexing::I
    layout::L
end

# ---- LatticeCore required interface ---------------------------------

LatticeCore.num_sites(lat::PeriodicLattice2D) = lat.N

LatticeCore.position(lat::PeriodicLattice2D, i::Int) = lat.positions[i]

LatticeCore.neighbors(lat::PeriodicLattice2D, i::Int) = lat.nearest_neighbors[i]

LatticeCore.boundary(lat::PeriodicLattice2D) = lat.boundary

LatticeCore.site_layout(lat::PeriodicLattice2D) = lat.layout

LatticeCore.size_trait(lat::PeriodicLattice2D) = FiniteSize((lat.Lx, lat.Ly))

# ---- Bond iteration --------------------------------------------------
# Override the default generic `bonds(lat)` fallback with the
# pre-computed, wrap-aware bond list.

LatticeCore.bonds(lat::PeriodicLattice2D) = lat.bonds

function LatticeCore.neighbor_bonds(lat::PeriodicLattice2D, i::Int)
    return (b for b in lat.bonds if b.i == i || b.j == i)
end

# ---- Sublattice / topology traits -----------------------------------

function LatticeCore.num_sublattices(lat::PeriodicLattice2D)
    length(get_unit_cell(typeof(lat.topology)).sublattice_positions)
end

LatticeCore.sublattice(lat::PeriodicLattice2D, i::Int) = lat.sublattice_ids[i]

LatticeCore.is_bipartite(lat::PeriodicLattice2D) = lat.is_bipartite_flag

# Periodicity follows the uniform rule: the lattice as a whole is
# periodic iff none of its axes is open. Mixed axis BCs → aperiodic.
function LatticeCore.periodicity(lat::PeriodicLattice2D)
    return all(!(ax isa OpenAxis) for ax in lat.boundary.axes) ? Periodic() : Aperiodic()
end

# Reciprocal support mirrors `periodicity`: only fully-periodic
# samples get a standard reciprocal lattice.
function LatticeCore.reciprocal_support(lat::PeriodicLattice2D)
    return if all(!(ax isa OpenAxis) for ax in lat.boundary.axes)
        HasReciprocal()
    else
        NoReciprocal()
    end
end

# Expose the topology name via the LatticeCore TopologyTrait.
function LatticeCore.topology(lat::PeriodicLattice2D)
    return TopologyTrait{topology_name(lat.topology)}()
end

topology_name(::Square) = :square
topology_name(::Triangular) = :triangular
topology_name(::Honeycomb) = :honeycomb
topology_name(::Kagome) = :kagome
topology_name(::Lieb) = :lieb
topology_name(::ShastrySutherland) = :shastry_sutherland
topology_name(::Dice) = :dice
topology_name(::UnionJack) = :union_jack

# ---- Basis / reciprocal ---------------------------------------------

"""
    basis_vectors(lat::PeriodicLattice2D) → SMatrix{2, 2, T}

Real-space primitive basis. Columns are the two primitive vectors
`a₁`, `a₂` of the underlying topology.
"""
LatticeCore.basis_vectors(lat::PeriodicLattice2D) = lat.basis_matrix

function LatticeCore.reciprocal_lattice(lat::PeriodicLattice2D{Topo,T}) where {Topo,T}
    reciprocal_support(lat) isa HasReciprocal || throw(
        ArgumentError(
            "PeriodicLattice2D has no reciprocal lattice unless every axis is periodic"
        ),
    )
    A = lat.basis_matrix
    B = SMatrix{2,2,T}(T(2π) * inv(transpose(A)))
    return monkhorst_pack(B, (lat.Lx, lat.Ly))
end

# ---- Coordinate conversions -----------------------------------------

function LatticeCore.to_real(
    lat::PeriodicLattice2D{Topo,T}, coord::LatticeCoord{2}
) where {Topo,T}
    cx, cy = coord.cell
    s = coord.sublattice
    uc = get_unit_cell(Topo)
    a1 = SVector{2,T}(T(uc.basis[1][1]), T(uc.basis[1][2]))
    a2 = SVector{2,T}(T(uc.basis[2][1]), T(uc.basis[2][2]))
    sub = SVector{2,T}(T(uc.sublattice_positions[s][1]), T(uc.sublattice_positions[s][2]))
    return RealSpace{2,T}((cx - 1) * a1 + (cy - 1) * a2 + sub)
end
