"""
    Diffusion Limited Aggregation (DLA) Module

This module implements Diffusion Limited Aggregation growth on 2D lattices.
DLA is a process where particles undergoing random walk stick together to form fractal structures.
"""

"""
    DLACluster

Represents a DLA cluster growing on a lattice.

# Fields
- `lattice::Lattice`: The underlying lattice
- `occupied_sites::Set{Int}`: Sites in the cluster
- `seed_site::Int`: Initial seed site
- `n_particles::Int`: Number of particles in cluster
- `spawn_radius::Float64`: Radius for spawning new particles
- `kill_radius::Float64`: Radius for killing escaped particles
"""
mutable struct DLACluster{T,B,I}
    lattice::Lattice{T,Float64,B,I}
    occupied_sites::Set{Int}
    seed_site::Int
    n_particles::Int
    spawn_radius::Float64
    kill_radius::Float64
end

export DLACluster

"""
    dla_cluster(lattice::Lattice, seed_site::Int=div(lattice.N, 2))

Initialize a DLA cluster with a single seed particle.

# Arguments
- `lattice::Lattice`: The lattice for DLA growth
- `seed_site::Int`: Seed site index (default: center of lattice)

# Returns
- `DLACluster`: The initialized DLA cluster

"""
function dla_cluster(lattice::Lattice{T,F,B,I}, seed_site::Int=div(lattice.N, 2)) where {T,F,B,I}
    @assert 1 <= seed_site <= lattice.N "Seed site must be within lattice bounds"
    
    return DLACluster{T,B,I}(
        lattice,
        Set([seed_site]),
        seed_site,
        1,
        5.0,  # initial spawn radius
        20.0  # initial kill radius
    )
end
export dla_cluster

"""
    update_radii!(dla::DLACluster)

Update spawn and kill radii based on cluster size.
"""
function update_radii!(dla::DLACluster)
    seed_pos = dla.lattice.positions[dla.seed_site]
    
    # Find maximum distance from seed to any cluster site
    max_dist = 0.0
    for site in dla.occupied_sites
        pos = dla.lattice.positions[site]
        dist = sqrt(sum((pos - seed_pos) .^ 2))
        max_dist = max(max_dist, dist)
    end
    
    dla.spawn_radius = max(5.0, max_dist + 3.0)
    dla.kill_radius = max(20.0, max_dist + 10.0)
end

"""
    spawn_particle(dla::DLACluster; rng=Random.default_rng())

Spawn a new particle at a random location outside the cluster.

# Returns
- `Int`: Site index of spawned particle (or 0 if no suitable site found)
"""
function spawn_particle(dla::DLACluster; rng=Random.default_rng())
    seed_pos = dla.lattice.positions[dla.seed_site]
    attempts = 0
    max_attempts = 1000
    
    while attempts < max_attempts
        # Random angle
        theta = 2π * rand(rng)
        
        # Spawn at spawn_radius distance
        spawn_pos = seed_pos + dla.spawn_radius * [cos(theta), sin(theta)]
        
        # Find nearest lattice site
        min_dist = Inf
        nearest_site = 0
        
        for site in 1:dla.lattice.N
            if site in dla.occupied_sites
                continue
            end
            
            site_pos = dla.lattice.positions[site]
            dist = sqrt(sum((site_pos - spawn_pos) .^ 2))
            
            if dist < min_dist
                min_dist = dist
                nearest_site = site
            end
        end
        
        if nearest_site > 0
            return nearest_site
        end
        
        attempts += 1
    end
    
    # Fallback: random unoccupied site
    for _ in 1:100
        site = rand(rng, 1:dla.lattice.N)
        if !(site in dla.occupied_sites)
            return site
        end
    end
    
    return 0
end

"""
    is_adjacent_to_cluster(dla::DLACluster, site::Int)

Check if a site is adjacent to the cluster.

# Returns
- `Bool`: true if site has a neighbor in the cluster
"""
function is_adjacent_to_cluster(dla::DLACluster, site::Int)
    for neighbor in dla.lattice.nearest_neighbors[site]
        if neighbor in dla.occupied_sites
            return true
        end
    end
    return false
end

"""
    is_too_far(dla::DLACluster, site::Int)

Check if a particle is too far from cluster (should be killed).

# Returns
- `Bool`: true if particle exceeded kill radius
"""
function is_too_far(dla::DLACluster, site::Int)
    seed_pos = dla.lattice.positions[dla.seed_site]
    site_pos = dla.lattice.positions[site]
    dist = sqrt(sum((site_pos - seed_pos) .^ 2))
    return dist > dla.kill_radius
end

"""
    add_particle!(dla::DLACluster; max_walk_steps::Int=10000, rng=Random.default_rng())

Add one particle to the DLA cluster through random walk.

# Arguments
- `dla::DLACluster`: The DLA cluster
- `max_walk_steps::Int`: Maximum walk steps before killing particle (default: 10000)
- `rng`: Random number generator (optional)

# Returns
- `Bool`: true if particle was successfully added, false if failed

"""
function add_particle!(dla::DLACluster; max_walk_steps::Int=10000, rng=Random.default_rng())
    # Spawn particle
    particle_site = spawn_particle(dla, rng=rng)
    
    if particle_site == 0
        return false
    end
    
    # Random walk until adjacent to cluster or too far
    for _ in 1:max_walk_steps
        # Check if adjacent to cluster
        if is_adjacent_to_cluster(dla, particle_site)
            # Stick to cluster
            push!(dla.occupied_sites, particle_site)
            dla.n_particles += 1
            update_radii!(dla)
            return true
        end
        
        # Check if too far
        if is_too_far(dla, particle_site)
            # Kill particle and respawn
            return false
        end
        
        # Random step
        neighbors = dla.lattice.nearest_neighbors[particle_site]
        # Exclude occupied neighbors
        available = filter(n -> !(n in dla.occupied_sites), neighbors)
        
        if isempty(available)
            # Trapped by cluster - stick here
            push!(dla.occupied_sites, particle_site)
            dla.n_particles += 1
            update_radii!(dla)
            return true
        end
        
        particle_site = rand(rng, available)
    end
    
    # Max steps exceeded
    return false
end
export add_particle!

"""
    grow!(dla::DLACluster, n_particles::Int; max_attempts::Int=10, rng=Random.default_rng())

Grow the DLA cluster by adding multiple particles.

# Arguments
- `dla::DLACluster`: The DLA cluster
- `n_particles::Int`: Target number of particles to add
- `max_attempts::Int`: Max attempts per particle (default: 10)
- `rng`: Random number generator (optional)

# Returns
- `Int`: Actual number of particles added

"""
function grow!(dla::DLACluster, n_particles::Int; max_attempts::Int=10, rng=Random.default_rng())
    added = 0
    
    for _ in 1:n_particles
        success = false
        for attempt in 1:max_attempts
            if add_particle!(dla, rng=rng)
                success = true
                break
            end
        end
        
        if success
            added += 1
        else
            # If we fail too many times, stop
            break
        end
    end
    
    return added
end
export grow!

"""
    cluster_radius(dla::DLACluster)

Calculate the radius of gyration of the cluster.

# Returns
- `Float64`: Radius of gyration

"""
function cluster_radius(dla::DLACluster)
    seed_pos = dla.lattice.positions[dla.seed_site]
    
    sum_sq_dist = 0.0
    for site in dla.occupied_sites
        pos = dla.lattice.positions[site]
        sum_sq_dist += sum((pos - seed_pos) .^ 2)
    end
    
    return sqrt(sum_sq_dist / dla.n_particles)
end
export cluster_radius

"""
    fractal_dimension(dla::DLACluster)

Estimate the fractal dimension using box-counting method.

# Returns
- `Float64`: Estimated fractal dimension

"""
function fractal_dimension(dla::DLACluster)
    if dla.n_particles < 10
        return NaN
    end
    
    # Get cluster positions
    positions = [dla.lattice.positions[site] for site in dla.occupied_sites]
    
    # Find bounding box
    xs = [p[1] for p in positions]
    ys = [p[2] for p in positions]
    x_min, x_max = minimum(xs), maximum(xs)
    y_min, y_max = minimum(ys), maximum(ys)
    
    # Box sizes
    box_sizes = [1.0, 2.0, 4.0, 8.0]
    counts = Float64[]
    sizes = Float64[]
    
    for box_size in box_sizes
        # Count boxes containing cluster points
        boxes = Set{Tuple{Int,Int}}()
        
        for pos in positions
            box_x = floor(Int, (pos[1] - x_min) / box_size)
            box_y = floor(Int, (pos[2] - y_min) / box_size)
            push!(boxes, (box_x, box_y))
        end
        
        if length(boxes) > 0
            push!(counts, Float64(length(boxes)))
            push!(sizes, box_size)
        end
    end
    
    if length(counts) < 2
        return NaN
    end
    
    # Linear regression on log-log plot
    log_sizes = log.(sizes)
    log_counts = log.(counts)
    
    # Fit log(count) = D * log(1/size) + const
    # D = -slope of log(count) vs log(size)
    n = length(log_sizes)
    mean_x = sum(log_sizes) / n
    mean_y = sum(log_counts) / n
    
    numerator = sum((log_sizes .- mean_x) .* (log_counts .- mean_y))
    denominator = sum((log_sizes .- mean_x) .^ 2)
    
    slope = numerator / denominator
    return -slope
end
export fractal_dimension

"""
    mass_center(dla::DLACluster)

Calculate the center of mass of the cluster.

# Returns
- `Vector{Float64}`: Center of mass position

"""
function mass_center(dla::DLACluster)
    if dla.n_particles == 0
        return [0.0, 0.0]
    end
    
    cm = zeros(2)
    for site in dla.occupied_sites
        cm += dla.lattice.positions[site]
    end
    
    return cm / dla.n_particles
end
export mass_center

"""
    radial_density(dla::DLACluster, n_bins::Int=20)

Calculate radial density distribution from seed.

# Arguments
- `dla::DLACluster`: The DLA cluster
- `n_bins::Int`: Number of radial bins (default: 20)

# Returns
- `Tuple{Vector{Float64}, Vector{Int}}`: (bin_edges, counts)

"""
function radial_density(dla::DLACluster, n_bins::Int=20)
    seed_pos = dla.lattice.positions[dla.seed_site]
    
    # Calculate distances
    distances = Float64[]
    for site in dla.occupied_sites
        pos = dla.lattice.positions[site]
        dist = sqrt(sum((pos - seed_pos) .^ 2))
        push!(distances, dist)
    end
    
    if isempty(distances)
        return Float64[], Int[]
    end
    
    # Create bins
    max_dist = maximum(distances)
    bin_edges = range(0, max_dist, length=n_bins+1)
    counts = zeros(Int, n_bins)
    
    # Count particles in each bin
    for dist in distances
        bin_idx = searchsortedfirst(bin_edges, dist) - 1
        bin_idx = clamp(bin_idx, 1, n_bins)
        counts[bin_idx] += 1
    end
    
    return collect(bin_edges), counts
end
export radial_density
