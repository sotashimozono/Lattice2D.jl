"""
    Visualization Module

This module provides enhanced visualization functions for lattices and applications.
"""

"""
    plot_percolation(perc::Union{BondPercolation, SitePercolation}; 
                    kwargs...)

Visualize a percolation configuration with clusters colored.

# Arguments
- `perc`: Percolation configuration
- `kwargs...`: Additional arguments passed to plot

# Returns
- `Plots.Plot`: The plot object

"""
function plot_percolation(perc::Union{BondPercolation, SitePercolation}; 
                         show_clusters::Bool=true,
                         markersize::Int=3,
                         kwargs...)
    lattice = perc.lattice
    
    # Plot lattice structure
    p = plot(; aspect_ratio=:equal, legend=false, grid=false, 
             framestyle=:box, kwargs...)
    
    if perc isa BondPercolation
        # Plot all bonds lightly
        for bond in lattice.bonds
            pos1 = lattice.positions[bond.src]
            pos2 = lattice.positions[bond.dst]
            plot!(p, [pos1[1], pos2[1]], [pos1[2], pos2[2]], 
                  color=:lightgray, alpha=0.3, linewidth=0.5)
        end
        
        # Plot occupied bonds
        for bond_idx in perc.occupied_bonds
            bond = lattice.bonds[bond_idx]
            pos1 = lattice.positions[bond.src]
            pos2 = lattice.positions[bond.dst]
            
            if show_clusters
                cluster_id = perc.cluster_map[bond.src]
                color_idx = mod(cluster_id, 10) + 1
                plot!(p, [pos1[1], pos2[1]], [pos1[2], pos2[2]], 
                      color=color_idx, linewidth=2, alpha=0.8)
            else
                plot!(p, [pos1[1], pos2[1]], [pos1[2], pos2[2]], 
                      color=:blue, linewidth=2, alpha=0.8)
            end
        end
    else  # SitePercolation
        # Plot all sites lightly
        xs = [pos[1] for pos in lattice.positions]
        ys = [pos[2] for pos in lattice.positions]
        scatter!(p, xs, ys, color=:lightgray, markersize=markersize, 
                alpha=0.3, markerstrokewidth=0)
        
        # Plot occupied sites
        if show_clusters
            for cluster_id in unique(values(perc.cluster_map))
                if cluster_id == 0
                    continue
                end
                cluster_sites = [s for s in perc.occupied_sites if perc.cluster_map[s] == cluster_id]
                xs_cluster = [lattice.positions[s][1] for s in cluster_sites]
                ys_cluster = [lattice.positions[s][2] for s in cluster_sites]
                
                color_idx = mod(cluster_id, 10) + 1
                scatter!(p, xs_cluster, ys_cluster, color=color_idx, 
                        markersize=markersize*2, alpha=0.8, markerstrokewidth=0)
            end
        else
            xs_occ = [lattice.positions[s][1] for s in perc.occupied_sites]
            ys_occ = [lattice.positions[s][2] for s in perc.occupied_sites]
            scatter!(p, xs_occ, ys_occ, color=:blue, markersize=markersize*2, 
                    alpha=0.8, markerstrokewidth=0)
        end
    end
    
    return p
end
export plot_percolation

"""
    plot_random_walk(rw::RandomWalk; kwargs...)

Visualize a random walk path on a lattice.

# Arguments
- `rw::RandomWalk`: The random walk to plot
- `kwargs...`: Additional arguments passed to plot

# Returns
- `Plots.Plot`: The plot object

"""
function plot_random_walk(rw::RandomWalk; 
                         show_lattice::Bool=true,
                         linewidth::Int=2,
                         markersize::Int=4,
                         kwargs...)
    lattice = rw.lattice
    
    p = plot(; aspect_ratio=:equal, legend=false, grid=false, 
             framestyle=:box, kwargs...)
    
    # Plot lattice lightly if requested
    if show_lattice
        xs = [pos[1] for pos in lattice.positions]
        ys = [pos[2] for pos in lattice.positions]
        scatter!(p, xs, ys, color=:lightgray, markersize=2, 
                alpha=0.3, markerstrokewidth=0)
    end
    
    # Plot walk path
    path_xs = [lattice.positions[site][1] for site in rw.path]
    path_ys = [lattice.positions[site][2] for site in rw.path]
    
    plot!(p, path_xs, path_ys, color=:blue, linewidth=linewidth, alpha=0.6)
    
    # Mark start and end
    scatter!(p, [path_xs[1]], [path_ys[1]], color=:green, 
            markersize=markersize*2, label="Start", markerstrokewidth=0)
    scatter!(p, [path_xs[end]], [path_ys[end]], color=:red, 
            markersize=markersize*2, label="End", markerstrokewidth=0)
    
    return p
end
export plot_random_walk

"""
    plot_dla(dla::DLACluster; kwargs...)

Visualize a DLA cluster.

# Arguments
- `dla::DLACluster`: The DLA cluster to plot
- `kwargs...`: Additional arguments passed to plot

# Returns
- `Plots.Plot`: The plot object

"""
function plot_dla(dla::DLACluster; 
                 color_by_distance::Bool=true,
                 markersize::Int=3,
                 kwargs...)
    lattice = dla.lattice
    
    p = plot(; aspect_ratio=:equal, legend=false, grid=false, 
             framestyle=:box, kwargs...)
    
    if color_by_distance
        seed_pos = lattice.positions[dla.seed_site]
        
        # Calculate distances from seed
        distances = Float64[]
        for site in dla.occupied_sites
            pos = lattice.positions[site]
            dist = sqrt(sum((pos - seed_pos) .^ 2))
            push!(distances, dist)
        end
        
        max_dist = maximum(distances)
        
        # Plot each site colored by distance
        for (i, site) in enumerate(dla.occupied_sites)
            pos = lattice.positions[site]
            color_val = distances[i] / max_dist
            scatter!(p, [pos[1]], [pos[2]], 
                    color=cgrad(:viridis)[color_val], 
                    markersize=markersize, markerstrokewidth=0)
        end
    else
        xs = [lattice.positions[s][1] for s in dla.occupied_sites]
        ys = [lattice.positions[s][2] for s in dla.occupied_sites]
        scatter!(p, xs, ys, color=:blue, markersize=markersize, 
                markerstrokewidth=0)
    end
    
    # Mark seed
    seed_pos = lattice.positions[dla.seed_site]
    scatter!(p, [seed_pos[1]], [seed_pos[2]], color=:red, 
            markersize=markersize*3, markershape=:star5, 
            markerstrokewidth=1, markerstrokecolor=:black)
    
    return p
end
export plot_dla

"""
    plot_spanning_tree(tree::SpanningTree; kwargs...)

Visualize a spanning tree.

# Arguments
- `tree::SpanningTree`: The spanning tree to plot
- `kwargs...`: Additional arguments passed to plot

# Returns
- `Plots.Plot`: The plot object

"""
function plot_spanning_tree(tree::SpanningTree; 
                           color_by_depth::Bool=true,
                           linewidth::Int=2,
                           markersize::Int=3,
                           kwargs...)
    lattice = tree.lattice
    
    p = plot(; aspect_ratio=:equal, legend=false, grid=false, 
             framestyle=:box, kwargs...)
    
    # Plot tree edges
    for edge in tree.tree_edges
        u, v = edge
        pos1 = lattice.positions[u]
        pos2 = lattice.positions[v]
        plot!(p, [pos1[1], pos2[1]], [pos1[2], pos2[2]], 
              color=:black, linewidth=linewidth, alpha=0.6)
    end
    
    # Plot vertices colored by depth
    if color_by_depth
        max_depth = tree_height(tree)
        
        for v in 1:lattice.N
            depth = tree_depth(tree, v)
            pos = lattice.positions[v]
            color_val = depth / max(max_depth, 1)
            scatter!(p, [pos[1]], [pos[2]], 
                    color=cgrad(:viridis)[color_val], 
                    markersize=markersize, markerstrokewidth=0)
        end
    else
        xs = [pos[1] for pos in lattice.positions]
        ys = [pos[2] for pos in lattice.positions]
        scatter!(p, xs, ys, color=:blue, markersize=markersize, 
                markerstrokewidth=0)
    end
    
    # Mark root
    root_pos = lattice.positions[tree.root]
    scatter!(p, [root_pos[1]], [root_pos[2]], color=:red, 
            markersize=markersize*3, markershape=:star5,
            markerstrokewidth=1, markerstrokecolor=:black)
    
    return p
end
export plot_spanning_tree

"""
    animate_random_walk(rw::RandomWalk, filename::String; 
                       fps::Int=10, kwargs...)

Create an animation of a random walk.

# Arguments
- `rw::RandomWalk`: The random walk to animate
- `filename::String`: Output filename (e.g., "walk.gif")
- `fps::Int`: Frames per second (default: 10)
- `kwargs...`: Additional arguments passed to plot

"""
function animate_random_walk(rw::RandomWalk, filename::String; 
                            fps::Int=10, kwargs...)
    lattice = rw.lattice
    
    # Create animation
    anim = @animate for i in 1:min(length(rw.path), 500)  # Limit frames
        p = plot(; aspect_ratio=:equal, legend=false, grid=false, 
                 framestyle=:box, kwargs...)
        
        # Plot lattice lightly
        xs = [pos[1] for pos in lattice.positions]
        ys = [pos[2] for pos in lattice.positions]
        scatter!(p, xs, ys, color=:lightgray, markersize=2, 
                alpha=0.3, markerstrokewidth=0)
        
        # Plot walk up to current frame
        path_xs = [lattice.positions[site][1] for site in rw.path[1:i]]
        path_ys = [lattice.positions[site][2] for site in rw.path[1:i]]
        
        plot!(p, path_xs, path_ys, color=:blue, linewidth=2, alpha=0.6)
        
        # Mark current position
        scatter!(p, [path_xs[end]], [path_ys[end]], color=:red, 
                markersize=6, markerstrokewidth=0)
        
        title!("Step $i / $(rw.n_steps)")
    end
    
    gif(anim, filename, fps=fps)
end
export animate_random_walk

"""
    plot_cluster_size_distribution(perc::Union{BondPercolation, SitePercolation}; 
                                   kwargs...)

Plot the cluster size distribution.

# Arguments
- `perc`: Percolation configuration
- `kwargs...`: Additional arguments passed to plot

# Returns
- `Plots.Plot`: The plot object

"""
function plot_cluster_size_distribution(perc::Union{BondPercolation, SitePercolation}; 
                                       logscale::Bool=true, kwargs...)
    dist = cluster_size_distribution(perc)
    
    sizes = sort(collect(keys(dist)))
    counts = [dist[s] for s in sizes]
    
    if logscale
        p = scatter(sizes, counts, xlabel="Cluster Size", ylabel="Count",
                   title="Cluster Size Distribution", 
                   xscale=:log10, yscale=:log10,
                   legend=false, markerstrokewidth=0, kwargs...)
    else
        p = scatter(sizes, counts, xlabel="Cluster Size", ylabel="Count",
                   title="Cluster Size Distribution",
                   legend=false, markerstrokewidth=0, kwargs...)
    end
    
    return p
end
export plot_cluster_size_distribution

"""
    plot_msd_vs_time(lattice::Lattice, n_steps::Int, n_walks::Int=10; kwargs...)

Plot mean squared displacement vs time for random walks.

# Arguments
- `lattice::Lattice`: The lattice
- `n_steps::Int`: Number of steps
- `n_walks::Int`: Number of walks to average (default: 10)
- `kwargs...`: Additional arguments passed to plot

# Returns
- `Plots.Plot`: The plot object

"""
function plot_msd_vs_time(lattice::Lattice, n_steps::Int, n_walks::Int=10; 
                         rng=Random.default_rng(), kwargs...)
    times = 1:n_steps
    msd_values = zeros(n_steps)
    
    for _ in 1:n_walks
        start_site = rand(rng, 1:lattice.N)
        rw = random_walk(lattice, start_site)
        
        for t in times
            step!(rw, rng=rng)
            msd_values[t] += mean_squared_displacement(rw)
        end
    end
    
    msd_values ./= n_walks
    
    p = plot(times, msd_values, xlabel="Time Steps", ylabel="MSD",
             title="Mean Squared Displacement vs Time",
             linewidth=2, legend=false, kwargs...)
    
    # Add theoretical line for 2D random walk (MSD ~ t)
    plot!(p, times, times, linestyle=:dash, color=:red, 
          label="Theory (MSD ~ t)", linewidth=1)
    
    return p
end
export plot_msd_vs_time

"""
    plot_radial_density(dla::DLACluster; kwargs...)

Plot the radial density distribution of a DLA cluster.

# Arguments
- `dla::DLACluster`: The DLA cluster
- `kwargs...`: Additional arguments passed to plot

# Returns
- `Plots.Plot`: The plot object

"""
function plot_radial_density(dla::DLACluster; n_bins::Int=20, kwargs...)
    bin_edges, counts = radial_density(dla, n_bins)
    
    # Use bin centers for plotting
    bin_centers = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in 1:n_bins]
    
    p = plot(bin_centers, counts, xlabel="Distance from Seed", 
             ylabel="Particle Count", title="Radial Density Distribution",
             linewidth=2, marker=:circle, legend=false, kwargs...)
    
    return p
end
export plot_radial_density
