"""
    AbstractLattice{D}
abstract type for lattices in D dimensions.
"""
abstract type AbstractLattice{D} end
"""
    AbstractLatticeConnection
abstract type for lattice connections (edges, bonds).
"""
abstract type AbstractLatticeConnection end

"""
    Bond
struct for a bond (edge) in the lattice.
- `src::Int`: start site index
- `dst::Int`: destination site index
- `type::Int`: type of the bond for categorization
- `vector::Vector{Float64}`: expresses the bond vector from src to dst
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
Connection rules within or between unit cells.
- `src_sub`: sublattice index of the start point (1, 2, ...)
- `dst_sub`: sublattice index of the end point
- `dx`, `dy`: relative cell position of the end point (0,0 means within the same unit cell)
- `type`: type of the connection
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
Geometric definition data of the lattice. Basically, the lattice is constructed based on this information.
The `get_unit_cell(::Type{T})` function retrieves the unit cell data corresponding to each topology.
"""
struct UnitCell{D,T} <: AbstractLattice{D}
    basis::Vector{Vector{T}}
    sublattice_positions::Vector{Vector{T}}
    connections::Vector{Connection}
end
export UnitCell
"""
    AbstractTopology
Abstract type for lattice topologies.
"""
abstract type AbstractTopology{D} <: AbstractLattice{D} end

"""
Lattice{Topology<:AbstractTopology, T, B<:AbstractBoundaryCondition, I<:AbstractIndexing}
It mainly represents 2-dimiensional lattice, but it can be used as 1-dimensional lattice as well.
- `Lx::Int`: x direction lattice size
- `Ly::Int`: y direction lattice size
- `N::Int`: total number of sites
- `positions::Vector{Vector{T}}`: position vectors of each site
- `nearest_neighbors::Vector{Vector{Int}}`: nearest neighbor indices for each site
- `bonds::Vector{Bond}`: list of bonds (edges) in the lattice
- `basis_vectors::Vector{Vector{T}}`: lattice basis vectors
- `reciprocal_vectors::Union{Vector{Vector{T}}, Nothing}`: reciprocal lattice vectors
- `sublattice_ids::Vector{Int}`: sublattice IDs of each site
- `is_bipartite::Bool`: whether the lattice is bipartite
- `site_map::Union{Matrix{Int}, Nothing}`: mapping of site indices on the lattice
- `translation_x::Vector{Int}`: x direction translation vector
- `translation_y::Vector{Int}`: y direction translation vector
- `boundary::B`: boundary condition
- `index_method::I`: indexing method
"""
struct Lattice{Topology<:AbstractTopology,T,B<:AbstractBoundaryCondition, I<:AbstractIndexing} <: AbstractLattice{2}
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
    index_method::I
end
export Lattice

