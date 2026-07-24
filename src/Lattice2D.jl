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
    plaquette_center,
    InfiniteSize,
    CellSite,
    CellBond,
    translation_vectors,
    num_basis_sites,
    basis_position,
    cell_bonds,
    site_orbits,
    bond_orbits,
    plaquette_orbits,
    element_orbits,
    element_orbit_position,
    cell_position,
    incident_cell_bonds,
    neighbors_at

# ---- Lattice2D source files -----------------------------------------

include("utils/edge_key.jl")
include("core/topology.jl")
include("core/unitcells.jl")
include("core/lattice.jl")
include("core/cache.jl")
include("core/element_api.jl")
include("core/constructor.jl")
include("core/translation_cell.jl")
include("api/predicates.jl")
include("utils/iterator.jl")
include("api/builders.jl")
include("api/dual.jl")
include("api/symmetry.jl")
include("api/scaling.jl")
include("api/blocking.jl")
include("disorder/dilution.jl")

# ---- Plots-extension stubs ------------------------------------------
#
# Concrete methods live in `ext/Lattice2DPlotsExt.jl` and are loaded
# automatically once `Plots` is in scope.

"""
    plot_bonds(lat::Lattice; bond_types=:all, color_by=:type, kwargs...)

Lattice2D-flavoured standalone bond plot. Returns a fresh
`Plots.Plot` that draws every bond in `lat` (optionally filtered by
`bond_types`) as a line segment, grouped — and therefore coloured —
by either the bond's `:type` tag or its rounded direction.

This function is a stub in the core module; the concrete method
lives in the `Lattice2DPlotsExt` package extension and is loaded
automatically once `Plots` is in scope.

See the `Lattice2DPlotsExt` module docstring for the full keyword
list and worked examples.
"""
function plot_bonds end

"""
    plot_state(lat::Lattice, state::AbstractVector;
               colormap=:viridis, marker_size=12, kwargs...) → Plots.Plot

Visualise a per-site `state::AbstractVector` (with `length(state) ==
num_sites(lat)`) as a coloured scatter on top of the lattice
geometry. Both continuous fields (energy density, charge density,
expectation values) and discrete labels (`Bool`, small-cardinality
`Int` spin / clock-model configurations) are supported through the
same call.

The concrete method lives in the `Lattice2DPlotsExt` package
extension and is loaded automatically once `Plots` is in scope. See
the `Lattice2DPlotsExt` module docstring for the full keyword list
and worked examples.
"""
function plot_state end

"""
    brillouin_zone(lat::Lattice; shell::Int=2) -> Vector{SVector{2,Float64}}

Compute the Brillouin zone of `lat` as the Wigner-Seitz cell of its
reciprocal lattice. Returns an ordered list of polygon vertices in
Cartesian k-space, walked counter-clockwise around the origin.

Throws `ArgumentError` if `lat` has any open axis (no reciprocal
lattice => no BZ).

The concrete method lives in the `Lattice2DPlotsExt` package
extension and is loaded automatically once `Plots` is in scope. See
the `Lattice2DPlotsExt` module docstring for the algorithm and the
`shell` keyword.
"""
function brillouin_zone end

"""
    high_symmetry_points(lat::Lattice) -> Dict{Symbol,SVector{2,Float64}}

Topology-keyed dictionary of high-symmetry points (Gamma, X, M, K,
...) in Cartesian k-coordinates. Currently populated for `Square`,
`Triangular`, and `Honeycomb`; other topologies fall back to a
singleton `:Gamma` entry.

The concrete method lives in the `Lattice2DPlotsExt` package
extension and is loaded automatically once `Plots` is in scope.
"""
function high_symmetry_points end

"""
    plot_brillouin_zone(lat::Lattice; show_mesh=false, ml=nothing,
                        show_high_symmetry=false, kwargs...) -> Plots.Plot

Plot the Brillouin zone of `lat` as a closed polygon. Optionally
overlays a momentum-lattice mesh on top of the BZ and labels
high-symmetry points.

The concrete method lives in the `Lattice2DPlotsExt` package
extension and is loaded automatically once `Plots` is in scope. See
the `Lattice2DPlotsExt` module docstring for the full keyword list.
"""
function plot_brillouin_zone end

# ---- Exports --------------------------------------------------------

# Lattice2D-local types and functions
export AbstractTopology, Connection, UnitCell, get_unit_cell, get_plaquette_rules
export Lattice
export build_lattice
export sublattice_layout
export num_bonds, num_plaquettes, bond_type
export plot_bonds
export plot_state
export brillouin_zone, high_symmetry_points, plot_brillouin_zone
# Backend-dispatched visualization — hosted in LatticeCore; re-exported so
# `using Lattice2D` exposes `plot_lattice(lat; backend = :plots | :makie)` and
# the Makie-only extras.
export plot_lattice, AbstractPlotBackend, PlotsBackend, MakieBackend, default_plot_backend
export makie_state, makie_structure_factor
export DilutedLattice, dilute_sites, dilute_bonds
export AVAILABLE_LATTICES
export Square, Triangular, Honeycomb, Kagome, Lieb, ShastrySutherland, Dice, UnionJack

# Convenience builders (thin wrappers around `build_lattice`)
export square, triangular, honeycomb, kagome, lieb, union_jack, dice, shastry_sutherland

# Structural predicates (`src/api/predicates.jl`)
export coordination_number, mean_coordination, bond_distances, shells

# Edge / bulk region API — generic over AbstractLattice, hosted in LatticeCore;
# re-exported so `using Lattice2D` keeps exposing it.
export edge_sites, bulk_sites, edge_bonds

# Dual-lattice construction (`src/api/dual.jl`)
export dual_lattice, dual_topology
export blocking, block_sublattice

# Point-group symmetry & site orbits (`src/api/symmetry.jl`)
export SymmetryOperation, symmetry_operations, symmetry_group_order
export symmetry_orbits, site_orbit
export is_rotation, is_reflection, is_inversion, rotation_angle

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

# ---- Translation-cell (unit-cell motif) layer ----
export InfiniteLattice
export CellSite, CellBond
export translation_vectors, num_basis_sites, basis_position, cell_bonds
export site_orbits, bond_orbits, plaquette_orbits
export element_orbits, element_orbit_position
export cell_position, incident_cell_bonds, neighbors_at

end # module Lattice2D
