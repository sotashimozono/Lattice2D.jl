module Lattice2D

using LinearAlgebra
using StaticArrays

# ---- LatticeCore imports --------------------------------------------
#
# Everything that was previously duplicated inside Lattice2D (the
# AbstractLattice hierarchy, Bond type, PBC / OBC singletons, and the
# RowMajor/ColMajor/Snake indexing strategies) now lives in
# LatticeCore. We import the symbols we rely on and re-export them so
# downstream users can write `using Lattice2D` and still get the
# vocabulary they need.

using LatticeCore

import LatticeCore:
    AbstractLattice,
    num_sites,
    position,
    positions,
    neighbors,
    boundary,
    bonds,
    neighbor_bonds,
    site_layout,
    site_type,
    sublattice,
    num_sublattices,
    size_trait,
    is_bipartite,
    periodicity,
    reciprocal_support,
    topology,
    basis_vectors,
    reciprocal_lattice,
    to_real,
    to_lattice,
    Bond,
    LatticeBoundary,
    AbstractAxisBC,
    PeriodicAxis,
    OpenAxis,
    TwistedAxis,
    NoModifier,
    SSD,
    apply_axis_bc,
    axis_phase,
    bond_weight,
    AbstractCoordinate,
    RealSpace,
    LatticeCoord,
    AbstractIndexing,
    RowMajor,
    ColMajor,
    Snake,
    site_index,
    lattice_coord,
    AbstractSiteType,
    IsingSite,
    PottsSite,
    XYSite,
    HeisenbergSite,
    EmptySite,
    AbstractSiteLayout,
    UniformLayout,
    SublatticeLayout,
    ExplicitLayout,
    TopologyTrait,
    Periodic,
    Aperiodic,
    HasReciprocal,
    HasFourierModule,
    NoReciprocal,
    FiniteSize,
    monkhorst_pack,
    gamma_centered,
    num_k_points,
    k_point,
    reciprocal_basis,
    momentum_lattice,
    fourier_module,
    structure_factor,
    bond_center,
    is_finite,
    state_type,
    random_state,
    zero_state,
    domain,
    element_type,
    AbstractLatticeElement,
    VertexCenter,
    BondCenter,
    PlaquetteCenter,
    CellCenter,
    num_elements,
    elements,
    element_position,
    element_positions,
    element_neighbors,
    incident,
    PlaquetteRule,
    Plaquette,
    plaquettes,
    neighbor_plaquettes,
    plaquette_center

# ---- Lattice2D source files -----------------------------------------

include("core/topology.jl")
include("core/unitcells.jl")
include("core/lattice.jl")
include("core/constructor.jl")
include("utils/iterator.jl")

# ---- Exports --------------------------------------------------------

# Lattice2D-local types and functions
export AbstractTopology, Connection, UnitCell, get_unit_cell
export Lattice
export build_lattice
export AVAILABLE_LATTICES
export Square, Triangular, Honeycomb, Kagome, Lieb, ShastrySutherland, Dice, UnionJack

# Re-exports from LatticeCore so `using Lattice2D` is self-sufficient
export AbstractLattice
export LatticeBoundary, AbstractAxisBC, PeriodicAxis, OpenAxis, TwistedAxis
export AbstractBoundaryModifier, NoModifier, SSD
export AbstractIndexing, RowMajor, ColMajor, Snake
export AbstractSiteLayout, UniformLayout, SublatticeLayout, ExplicitLayout
export AbstractSiteType, IsingSite, PottsSite, XYSite, HeisenbergSite, EmptySite
export AbstractCoordinate, RealSpace, LatticeCoord, HigherDimCoord
export num_sites, position, positions, neighbors, boundary, bonds, neighbor_bonds
export site_layout, site_type, sublattice, num_sublattices
export site_index, lattice_coord
export is_bipartite, topology
export TopologyTrait
export Periodic, Aperiodic, periodicity
export AbstractReciprocalSupport, HasReciprocal, HasFourierModule, NoReciprocal
export reciprocal_support
export AbstractSizeTrait, FiniteSize, InfiniteSize, QuasiInfiniteSize
export size_trait, is_finite
export basis_vectors, reciprocal_lattice, fourier_module, momentum_lattice
export num_k_points, k_point, reciprocal_basis
export to_real, to_lattice, to_hyper
export monkhorst_pack, gamma_centered, structure_factor
export Bond, bond_center
export apply_axis_bc, axis_phase, bond_weight
export AbstractLatticeElement, VertexCenter, BondCenter, PlaquetteCenter, CellCenter
export element_type, state_type, random_state, zero_state, domain
export num_elements, elements, element_position, element_positions
export element_neighbors, incident
export PlaquetteRule, Plaquette, plaquettes, neighbor_plaquettes, plaquette_center
export materialize, require_finite

end # module Lattice2D
