"""
Quasicrystal structures and generation methods.
This file contains the infrastructure for generating quasicrystalline patterns.
"""

# Mathematical constants commonly used in quasicrystals
const GOLDEN_RATIO = (1 + sqrt(5)) / 2
const PHI = GOLDEN_RATIO

"""
    AbstractQuasicrystal{D} <: AbstractTopology{D}
Abstract type for quasicrystal topologies in D dimensions.
Unlike periodic lattices, quasicrystals lack translational symmetry but have long-range order.
"""
abstract type AbstractQuasicrystal{D} <: AbstractTopology{D} end

"""
    AbstractGenerationMethod
Abstract type for quasicrystal generation methods.
Different methods can be used to generate the same quasicrystal pattern.
"""
abstract type AbstractGenerationMethod end

"""
    ProjectionMethod <: AbstractGenerationMethod
Generate quasicrystals via projection from higher-dimensional periodic lattices.
This is the most common theoretical approach.
"""
struct ProjectionMethod <: AbstractGenerationMethod end

"""
    SubstitutionMethod <: AbstractGenerationMethod
Generate quasicrystals via substitution (inflation) rules.
This is an alternative algorithmic approach.
"""
struct SubstitutionMethod <: AbstractGenerationMethod end

"""
    QuasicrystalData{D,T}
Data structure holding the generated quasicrystal pattern.
- `positions::Vector{Vector{T}}`: positions of vertices/sites
- `tiles::Vector{Tile{D,T}}`: list of tiles in the pattern
- `generation_method::AbstractGenerationMethod`: method used to generate
- `parameters::Dict{Symbol,Any}`: generation parameters
"""
struct QuasicrystalData{D,T}
    positions::Vector{Vector{T}}
    tiles::Vector{Any}  # Will be specialized based on quasicrystal type
    generation_method::AbstractGenerationMethod
    parameters::Dict{Symbol,Any}
end

"""
    Tile{D,T}
Represents a single tile in the quasicrystal pattern.
- `vertices::Vector{Vector{T}}`: corner positions of the tile
- `type::Int`: tile type identifier (e.g., fat vs thin rhombus)
- `center::Vector{T}`: center position of the tile
"""
struct Tile{D,T}
    vertices::Vector{Vector{T}}
    type::Int
    center::Vector{T}
end

export AbstractQuasicrystal, AbstractGenerationMethod
export ProjectionMethod, SubstitutionMethod
export QuasicrystalData, Tile
export GOLDEN_RATIO, PHI
