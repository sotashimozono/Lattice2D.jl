"""
    Random Walk Module

This module implements random walks on 2D lattices with various statistics.
"""

"""
    RandomWalk

Represents a random walk on a lattice.

# Fields
- `lattice::Lattice`: The underlying lattice
- `path::Vector{Int}`: Sequence of visited sites
- `current_site::Int`: Current position
- `start_site::Int`: Starting position
- `n_steps::Int`: Number of steps taken
"""
mutable struct RandomWalk{T,B,I}
    lattice::Lattice{T,Float64,B,I}
    path::Vector{Int}
    current_site::Int
    start_site::Int
    n_steps::Int
end

export RandomWalk

"""
    random_walk(lattice::Lattice, start_site::Int=1)

Initialize a random walk on a lattice.

# Arguments
- `lattice::Lattice`: The lattice for the random walk
- `start_site::Int`: Starting site index (default: 1)

# Returns
- `RandomWalk`: The random walk object

"""
function random_walk(lattice::Lattice{T,F,B,I}, start_site::Int=1) where {T,F,B,I}
    @assert 1 <= start_site <= lattice.N "Start site must be within lattice bounds"
    
    return RandomWalk{T,B,I}(
        lattice,
        [start_site],
        start_site,
        start_site,
        0
    )
end
export random_walk

"""
    step!(rw::RandomWalk; rng=Random.default_rng())

Take one random step in the walk.

# Arguments
- `rw::RandomWalk`: The random walk to update
- `rng`: Random number generator (optional)

# Returns
- `Int`: The new current site

"""
function step!(rw::RandomWalk; rng=Random.default_rng())
    neighbors = rw.lattice.nearest_neighbors[rw.current_site]
    @assert !isempty(neighbors) "Current site has no neighbors"
    
    # Choose a random neighbor
    next_site = rand(rng, neighbors)
    rw.current_site = next_site
    push!(rw.path, next_site)
    rw.n_steps += 1
    
    return next_site
end
export step!

"""
    walk!(rw::RandomWalk, n_steps::Int; rng=Random.default_rng())

Take multiple random steps.

# Arguments
- `rw::RandomWalk`: The random walk to update
- `n_steps::Int`: Number of steps to take
- `rng`: Random number generator (optional)

"""
function walk!(rw::RandomWalk, n_steps::Int; rng=Random.default_rng())
    for _ in 1:n_steps
        step!(rw, rng=rng)
    end
    return rw
end
export walk!

"""
    displacement(rw::RandomWalk)

Calculate the displacement vector from start to current position.

# Returns
- `Vector{Float64}`: Displacement vector

"""
function displacement(rw::RandomWalk)
    start_pos = rw.lattice.positions[rw.start_site]
    current_pos = rw.lattice.positions[rw.current_site]
    return current_pos - start_pos
end
export displacement

"""
    mean_squared_displacement(rw::RandomWalk)

Calculate the mean squared displacement (MSD) from the start.

# Returns
- `Float64`: Mean squared displacement

"""
function mean_squared_displacement(rw::RandomWalk)
    disp = displacement(rw)
    return sum(disp .^ 2)
end
export mean_squared_displacement

"""
    end_to_end_distance(rw::RandomWalk)

Calculate the Euclidean distance from start to current position.

# Returns
- `Float64`: End-to-end distance

"""
function end_to_end_distance(rw::RandomWalk)
    return sqrt(mean_squared_displacement(rw))
end
export end_to_end_distance

"""
    path_length(rw::RandomWalk)

Calculate the total path length traveled.

# Returns
- `Float64`: Total path length

"""
function path_length(rw::RandomWalk)
    total_length = 0.0
    for i in 2:length(rw.path)
        prev_site = rw.path[i-1]
        curr_site = rw.path[i]
        prev_pos = rw.lattice.positions[prev_site]
        curr_pos = rw.lattice.positions[curr_site]
        total_length += sqrt(sum((curr_pos - prev_pos) .^ 2))
    end
    return total_length
end
export path_length

"""
    visited_sites(rw::RandomWalk)

Return the set of unique sites visited during the walk.

# Returns
- `Set{Int}`: Set of visited site indices

"""
function visited_sites(rw::RandomWalk)
    return Set(rw.path)
end
export visited_sites

"""
    return_probability(rw::RandomWalk, site::Int)

Calculate the fraction of time spent at a given site.

# Arguments
- `rw::RandomWalk`: The random walk
- `site::Int`: Site index

# Returns
- `Float64`: Return probability (0 to 1)

"""
function return_probability(rw::RandomWalk, site::Int)
    count = sum(rw.path .== site)
    return count / length(rw.path)
end
export return_probability

"""
    self_avoiding_walk(lattice::Lattice, max_steps::Int, start_site::Int=1; rng=Random.default_rng())

Perform a self-avoiding random walk (SAW).

# Arguments
- `lattice::Lattice`: The lattice for the walk
- `max_steps::Int`: Maximum number of steps to attempt
- `start_site::Int`: Starting site index (default: 1)
- `rng`: Random number generator (optional)

# Returns
- `RandomWalk`: The self-avoiding walk (may be shorter than max_steps if trapped)

"""
function self_avoiding_walk(lattice::Lattice{T,F,B,I}, max_steps::Int, start_site::Int=1; 
                           rng=Random.default_rng()) where {T,F,B,I}
    @assert 1 <= start_site <= lattice.N "Start site must be within lattice bounds"
    @assert max_steps > 0 "max_steps must be positive"
    
    rw = RandomWalk{T,B,I}(
        lattice,
        [start_site],
        start_site,
        start_site,
        0
    )
    
    visited = Set([start_site])
    
    for _ in 1:max_steps
        neighbors = lattice.nearest_neighbors[rw.current_site]
        
        # Filter out visited sites
        available = filter(n -> !(n in visited), neighbors)
        
        if isempty(available)
            # Trapped - cannot continue
            break
        end
        
        # Choose random available neighbor
        next_site = rand(rng, available)
        rw.current_site = next_site
        push!(rw.path, next_site)
        push!(visited, next_site)
        rw.n_steps += 1
    end
    
    return rw
end
export self_avoiding_walk

"""
    average_msd(lattice::Lattice, n_steps::Int, n_walks::Int=100; rng=Random.default_rng())

Calculate average mean squared displacement over multiple walks.

# Arguments
- `lattice::Lattice`: The lattice
- `n_steps::Int`: Number of steps per walk
- `n_walks::Int`: Number of walks to average (default: 100)
- `rng`: Random number generator (optional)

# Returns
- `Float64`: Average MSD

"""
function average_msd(lattice::Lattice, n_steps::Int, n_walks::Int=100; rng=Random.default_rng())
    total_msd = 0.0
    
    for _ in 1:n_walks
        start_site = rand(rng, 1:lattice.N)
        rw = random_walk(lattice, start_site)
        walk!(rw, n_steps, rng=rng)
        total_msd += mean_squared_displacement(rw)
    end
    
    return total_msd / n_walks
end
export average_msd

"""
    mixing_time(lattice::Lattice, n_steps::Int, target_coverage::Float64=0.9; 
                n_trials::Int=10, rng=Random.default_rng())

Estimate the mixing time - average steps needed to visit a fraction of sites.

# Arguments
- `lattice::Lattice`: The lattice
- `n_steps::Int`: Maximum steps to try
- `target_coverage::Float64`: Target fraction of sites to visit (default: 0.9)
- `n_trials::Int`: Number of trials to average (default: 10)
- `rng`: Random number generator (optional)

# Returns
- `Float64`: Average number of steps to reach target coverage (or Inf if not reached)

"""
function mixing_time(lattice::Lattice, n_steps::Int, target_coverage::Float64=0.9;
                    n_trials::Int=10, rng=Random.default_rng())
    @assert 0.0 < target_coverage <= 1.0 "target_coverage must be in (0, 1]"
    
    target_sites = ceil(Int, target_coverage * lattice.N)
    times = Float64[]
    
    for _ in 1:n_trials
        start_site = rand(rng, 1:lattice.N)
        rw = random_walk(lattice, start_site)
        visited = Set([start_site])
        
        for step in 1:n_steps
            step!(rw, rng=rng)
            push!(visited, rw.current_site)
            
            if length(visited) >= target_sites
                push!(times, Float64(step))
                break
            end
        end
    end
    
    return isempty(times) ? Inf : sum(times) / length(times)
end
export mixing_time
