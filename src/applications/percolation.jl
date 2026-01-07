"""
    Percolation Module

This module implements bond and site percolation on 2D lattices.
Percolation theory studies the behavior of connected clusters in random graphs.
"""

"""
    BondPercolation

Represents a bond percolation configuration on a lattice.

# Fields
- `lattice::Lattice`: The underlying lattice
- `occupied_bonds::Set{Int}`: Indices of occupied bonds
- `probability::Float64`: Bond occupation probability
- `clusters::Vector{Vector{Int}}`: List of connected clusters (cached)
- `cluster_map::Vector{Int}`: Maps each site to its cluster ID
"""
mutable struct BondPercolation{T,B,I}
    lattice::Lattice{T,Float64,B,I}
    occupied_bonds::Set{Int}
    probability::Float64
    clusters::Vector{Vector{Int}}
    cluster_map::Vector{Int}
end

"""
    SitePercolation

Represents a site percolation configuration on a lattice.

# Fields
- `lattice::Lattice`: The underlying lattice
- `occupied_sites::Set{Int}`: Indices of occupied sites
- `probability::Float64`: Site occupation probability
- `clusters::Vector{Vector{Int}}`: List of connected clusters (cached)
- `cluster_map::Vector{Int}`: Maps each site to its cluster ID
"""
mutable struct SitePercolation{T,B,I}
    lattice::Lattice{T,Float64,B,I}
    occupied_sites::Set{Int}
    probability::Float64
    clusters::Vector{Vector{Int}}
    cluster_map::Vector{Int}
end

export BondPercolation, SitePercolation

"""
    bond_percolation(lattice::Lattice, p::Float64; rng=Random.default_rng())

Create a bond percolation configuration with occupation probability `p`.

# Arguments
- `lattice::Lattice`: The lattice to perform percolation on
- `p::Float64`: Bond occupation probability (0 ≤ p ≤ 1)
- `rng`: Random number generator (optional)

# Returns
- `BondPercolation`: The percolation configuration

"""
function bond_percolation(lattice::Lattice{T,F,B,I}, p::Float64; rng=Random.default_rng()) where {T,F,B,I}
    @assert 0.0 <= p <= 1.0 "Probability must be between 0 and 1"
    
    n_bonds = length(lattice.bonds)
    occupied = Set{Int}()
    
    for i in 1:n_bonds
        if rand(rng) < p
            push!(occupied, i)
        end
    end
    
    perc = BondPercolation{T,B,I}(
        lattice,
        occupied,
        p,
        Vector{Int}[],
        zeros(Int, lattice.N)
    )
    
    find_clusters!(perc)
    return perc
end
export bond_percolation

"""
    site_percolation(lattice::Lattice, p::Float64; rng=Random.default_rng())

Create a site percolation configuration with occupation probability `p`.

# Arguments
- `lattice::Lattice`: The lattice to perform percolation on
- `p::Float64`: Site occupation probability (0 ≤ p ≤ 1)
- `rng`: Random number generator (optional)

# Returns
- `SitePercolation`: The percolation configuration

"""
function site_percolation(lattice::Lattice{T,F,B,I}, p::Float64; rng=Random.default_rng()) where {T,F,B,I}
    @assert 0.0 <= p <= 1.0 "Probability must be between 0 and 1"
    
    occupied = Set{Int}()
    for i in 1:lattice.N
        if rand(rng) < p
            push!(occupied, i)
        end
    end
    
    perc = SitePercolation{T,B,I}(
        lattice,
        occupied,
        p,
        Vector{Int}[],
        zeros(Int, lattice.N)
    )
    
    find_clusters!(perc)
    return perc
end
export site_percolation

"""
    find_clusters!(perc::BondPercolation)

Find all connected clusters in a bond percolation configuration using Union-Find algorithm.
Updates the `clusters` and `cluster_map` fields.
"""
function find_clusters!(perc::BondPercolation)
    N = perc.lattice.N
    parent = collect(1:N)
    
    # Union-Find helper functions
    function find_root(x::Int)
        if parent[x] != x
            parent[x] = find_root(parent[x])  # path compression
        end
        return parent[x]
    end
    
    function union!(x::Int, y::Int)
        root_x = find_root(x)
        root_y = find_root(y)
        if root_x != root_y
            parent[root_y] = root_x
        end
    end
    
    # Connect sites through occupied bonds
    for bond_idx in perc.occupied_bonds
        bond = perc.lattice.bonds[bond_idx]
        union!(bond.src, bond.dst)
    end
    
    # Build clusters
    cluster_dict = Dict{Int, Vector{Int}}()
    for i in 1:N
        root = find_root(i)
        if !haskey(cluster_dict, root)
            cluster_dict[root] = Int[]
        end
        push!(cluster_dict[root], i)
        perc.cluster_map[i] = root
    end
    
    perc.clusters = collect(values(cluster_dict))
    return perc
end

"""
    find_clusters!(perc::SitePercolation)

Find all connected clusters in a site percolation configuration.
Only considers connections between occupied sites.
Updates the `clusters` and `cluster_map` fields.
"""
function find_clusters!(perc::SitePercolation)
    N = perc.lattice.N
    parent = collect(1:N)
    
    function find_root(x::Int)
        if parent[x] != x
            parent[x] = find_root(parent[x])
        end
        return parent[x]
    end
    
    function union!(x::Int, y::Int)
        root_x = find_root(x)
        root_y = find_root(y)
        if root_x != root_y
            parent[root_y] = root_x
        end
    end
    
    # Only connect occupied sites that are neighbors
    for site in perc.occupied_sites
        for neighbor in perc.lattice.nearest_neighbors[site]
            if neighbor in perc.occupied_sites
                union!(site, neighbor)
            end
        end
    end
    
    # Build clusters (only from occupied sites)
    cluster_dict = Dict{Int, Vector{Int}}()
    for site in perc.occupied_sites
        root = find_root(site)
        if !haskey(cluster_dict, root)
            cluster_dict[root] = Int[]
        end
        push!(cluster_dict[root], site)
        perc.cluster_map[site] = root
    end
    
    perc.clusters = collect(values(cluster_dict))
    return perc
end

"""
    largest_cluster_size(perc::Union{BondPercolation, SitePercolation})

Return the size of the largest cluster.

"""
function largest_cluster_size(perc::Union{BondPercolation, SitePercolation})
    isempty(perc.clusters) && return 0
    return maximum(length(c) for c in perc.clusters)
end
export largest_cluster_size

"""
    number_of_clusters(perc::Union{BondPercolation, SitePercolation})

Return the total number of clusters.

"""
function number_of_clusters(perc::Union{BondPercolation, SitePercolation})
    return length(perc.clusters)
end
export number_of_clusters

"""
    cluster_size_distribution(perc::Union{BondPercolation, SitePercolation})

Return a dictionary mapping cluster sizes to their frequencies.

# Returns
- `Dict{Int, Int}`: Maps cluster size to number of clusters of that size

"""
function cluster_size_distribution(perc::Union{BondPercolation, SitePercolation})
    dist = Dict{Int, Int}()
    for cluster in perc.clusters
        size = length(cluster)
        dist[size] = get(dist, size, 0) + 1
    end
    return dist
end
export cluster_size_distribution

"""
    spans_lattice(perc::Union{BondPercolation, SitePercolation}; direction=:both)

Check if any cluster spans the lattice in the given direction.

# Arguments
- `perc`: Percolation configuration
- `direction`: `:horizontal`, `:vertical`, or `:both` (default)

# Returns
- `Bool`: true if a spanning cluster exists

"""
function spans_lattice(perc::Union{BondPercolation, SitePercolation}; direction::Symbol=:both)
    lattice = perc.lattice
    Lx, Ly = lattice.Lx, lattice.Ly
    
    check_horizontal = direction in [:horizontal, :both]
    check_vertical = direction in [:vertical, :both]
    
    # Count number of unique sublattice sites per unit cell
    n_sub = maximum(lattice.sublattice_ids)
    
    for cluster in perc.clusters
        if check_horizontal
            x_coords = Set{Int}()
            for site in cluster
                # Find x coordinate from site_map
                for x in 1:Lx, y in 1:Ly
                    if lattice.site_map[x, y] <= site <= lattice.site_map[x, y] + n_sub - 1
                        push!(x_coords, x)
                        break
                    end
                end
            end
            if length(x_coords) == Lx  # spans all x coordinates
                return true
            end
        end
        
        if check_vertical
            y_coords = Set{Int}()
            for site in cluster
                for x in 1:Lx, y in 1:Ly
                    if lattice.site_map[x, y] <= site <= lattice.site_map[x, y] + n_sub - 1
                        push!(y_coords, y)
                        break
                    end
                end
            end
            if length(y_coords) == Ly  # spans all y coordinates
                return true
            end
        end
    end
    
    return false
end
export spans_lattice

"""
    percolation_strength(perc::Union{BondPercolation, SitePercolation})

Return the fraction of sites in the largest cluster.

# Returns
- `Float64`: Percolation strength (0 to 1)

"""
function percolation_strength(perc::Union{BondPercolation, SitePercolation})
    N_total = perc isa BondPercolation ? perc.lattice.N : length(perc.occupied_sites)
    return N_total > 0 ? largest_cluster_size(perc) / N_total : 0.0
end
export percolation_strength
