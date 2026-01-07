# Applications

Lattice2D.jl provides several application examples that demonstrate how to use 2D lattices for various computational physics simulations.

## Available Applications

### 1. Percolation Theory

Bond and site percolation simulations with cluster analysis.

**Features:**
- Bond percolation
- Site percolation  
- Cluster identification using Union-Find algorithm
- Spanning detection (horizontal, vertical, or both)
- Cluster size distribution
- Percolation strength calculation

**Main Functions:**
- `bond_percolation(lattice, p)` - Create bond percolation configuration
- `site_percolation(lattice, p)` - Create site percolation configuration
- `largest_cluster_size(perc)` - Get largest cluster size
- `spans_lattice(perc)` - Check if any cluster spans the lattice
- `cluster_size_distribution(perc)` - Get cluster size histogram

### 2. Random Walks

Random walk simulations on lattices with statistical analysis.

**Features:**
- Standard random walk
- Self-avoiding random walk (SAW)
- Mean squared displacement (MSD) tracking
- Path length calculation
- Visited sites tracking
- Return probability
- Mixing time estimation

**Main Functions:**
- `random_walk(lattice, start_site)` - Initialize random walk
- `step!(rw)` - Take one random step
- `walk!(rw, n_steps)` - Take multiple steps
- `self_avoiding_walk(lattice, max_steps)` - Generate SAW
- `mean_squared_displacement(rw)` - Calculate MSD
- `average_msd(lattice, n_steps, n_walks)` - Average MSD over multiple walks

### 3. Diffusion Limited Aggregation (DLA)

Fractal growth through diffusion-limited aggregation.

**Features:**
- DLA cluster growth
- Fractal dimension estimation
- Radial density distribution
- Center of mass calculation
- Adaptive spawn and kill radii

**Main Functions:**
- `dla_cluster(lattice, seed_site)` - Initialize DLA cluster
- `add_particle!(dla)` - Add one particle
- `grow!(dla, n_particles)` - Grow cluster by n particles
- `cluster_radius(dla)` - Calculate radius of gyration
- `fractal_dimension(dla)` - Estimate fractal dimension
- `radial_density(dla, n_bins)` - Get radial density profile

### 4. Spanning Trees

Uniform random spanning tree generation using Wilson's algorithm.

**Features:**
- Wilson's algorithm for uniform spanning trees
- Tree distance calculations
- Depth and height calculations
- Leaf vertices identification
- Subtree size calculations
- Maze generation

**Main Functions:**
- `wilson_algorithm(lattice, root)` - Generate uniform random spanning tree
- `tree_depth(tree, vertex)` - Get depth from root
- `tree_distance(tree, u, v)` - Distance between vertices in tree
- `path_to_root(tree, vertex)` - Find path to root
- `leaf_vertices(tree)` - Get all leaf vertices
- `generate_maze(lattice)` - Generate maze from spanning tree

### 5. Visualization

Plotting functions for lattices and all applications.

**Features:**
- Percolation visualization with cluster coloring
- Random walk path visualization
- DLA cluster visualization with distance coloring
- Spanning tree visualization with depth coloring
- Cluster size distribution plots
- MSD vs time plots
- Radial density plots
- Animation support for random walks

**Main Functions:**
- `plot_percolation(perc)` - Visualize percolation configuration
- `plot_random_walk(rw)` - Visualize random walk path
- `plot_dla(dla)` - Visualize DLA cluster
- `plot_spanning_tree(tree)` - Visualize spanning tree
- `plot_cluster_size_distribution(perc)` - Plot cluster sizes
- `plot_msd_vs_time(lattice, n_steps, n_walks)` - Plot MSD evolution
- `animate_random_walk(rw, filename)` - Create animation

## Usage Examples

### Percolation Example

```julia
using Lattice2D

# Create lattice
lat = build_lattice(Square, 50, 50)

# Bond percolation at critical probability
perc = bond_percolation(lat, 0.5)

# Analyze clusters
println("Number of clusters: ", number_of_clusters(perc))
println("Largest cluster: ", largest_cluster_size(perc))
println("Spans lattice: ", spans_lattice(perc))

# Visualize
using Plots
plot_percolation(perc)
```

### Random Walk Example

```julia
using Lattice2D

# Create lattice
lat = build_lattice(Triangular, 30, 30)

# Perform random walk
rw = random_walk(lat)
walk!(rw, 1000)

# Calculate statistics
println("MSD: ", mean_squared_displacement(rw))
println("Visited sites: ", length(visited_sites(rw)))

# Visualize
plot_random_walk(rw)
```

### DLA Example

```julia
using Lattice2D

# Create lattice
lat = build_lattice(Square, 100, 100)

# Grow DLA cluster
dla = dla_cluster(lat)
grow!(dla, 500)

# Analyze
println("Cluster radius: ", cluster_radius(dla))
println("Fractal dimension: ", fractal_dimension(dla))

# Visualize
plot_dla(dla)
```

### Spanning Tree Example

```julia
using Lattice2D

# Create lattice
lat = build_lattice(Honeycomb, 20, 20)

# Generate spanning tree
tree = wilson_algorithm(lat)

# Analyze
println("Tree height: ", tree_height(tree))
println("Number of leaves: ", length(leaf_vertices(tree)))

# Visualize
plot_spanning_tree(tree)
```

## Performance Considerations

- For large lattices (>100x100), DLA growth can be slow. Consider using smaller lattices or fewer particles.
- Random walk statistics are more accurate with more samples (higher `n_walks` parameter).
- Spanning tree generation time is O(N) where N is the number of sites.
- Visualization can be slow for very large structures. Consider reducing marker/line sizes or plotting subsets.

## References

- Percolation: Stauffer & Aharony, "Introduction to Percolation Theory"
- DLA: Witten & Sander, Phys. Rev. Lett. 47, 1400 (1981)
- Spanning Trees: Wilson, JACM 43, 37-82 (1996)
