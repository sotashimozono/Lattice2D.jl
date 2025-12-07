"""
Visualization utilities for quasicrystals.
"""

using Plots

"""
    visualize_quasicrystal_positions(qc_data::QuasicrystalData{1,T}; kwargs...) where T
Visualize 1D quasicrystal positions as points on a line.
"""
function visualize_quasicrystal_positions(qc_data::QuasicrystalData{1,T}; kwargs...) where T
    positions_1d = [p[1] for p in qc_data.positions]
    y_vals = zeros(length(positions_1d))
    
    p = scatter(positions_1d, y_vals,
        marker=:circle,
        markersize=6,
        label="Sites",
        ylims=(-0.5, 0.5),
        yticks=[],
        xlabel="Position",
        title="1D Quasicrystal",
        legend=:topright;
        kwargs...
    )
    
    return p
end

"""
    visualize_quasicrystal_positions(qc_data::QuasicrystalData{2,T}; kwargs...) where T
Visualize 2D quasicrystal positions as scattered points.
"""
function visualize_quasicrystal_positions(qc_data::QuasicrystalData{2,T}; kwargs...) where T
    x_vals = [p[1] for p in qc_data.positions]
    y_vals = [p[2] for p in qc_data.positions]
    
    p = scatter(x_vals, y_vals,
        marker=:circle,
        markersize=3,
        label="Vertices",
        aspect_ratio=:equal,
        xlabel="x",
        ylabel="y",
        title="2D Quasicrystal",
        legend=:topright;
        kwargs...
    )
    
    return p
end

"""
    visualize_quasicrystal_tiles(qc_data::QuasicrystalData{2,T}; kwargs...) where T
Visualize 2D quasicrystal tiles.
"""
function visualize_quasicrystal_tiles(qc_data::QuasicrystalData{2,T};  kwargs...) where T
    if isempty(qc_data.tiles)
        @warn "No tiles to visualize. Use visualize_quasicrystal_positions instead."
        return visualize_quasicrystal_positions(qc_data; kwargs...)
    end
    
    p = plot(aspect_ratio=:equal, xlabel="x", ylabel="y", 
             title="2D Quasicrystal Tiling", legend=false; kwargs...)
    
    # Group tiles by type
    tile_types = unique([t.type for t in qc_data.tiles])
    colors = distinguishable_colors(length(tile_types) + 1)[2:end]  # Skip white
    
    for (idx, tile_type) in enumerate(tile_types)
        tiles_of_type = filter(t -> t.type == tile_type, qc_data.tiles)
        
        for tile in tiles_of_type
            # Close the polygon by adding the first vertex at the end
            vertices = vcat(tile.vertices, [tile.vertices[1]])
            x_coords = [v[1] for v in vertices]
            y_coords = [v[2] for v in vertices]
            
            plot!(p, x_coords, y_coords,
                fillcolor=colors[idx],
                fillalpha=0.3,
                linecolor=:black,
                linewidth=1,
                label=(idx == 1 ? "Type $tile_type" : "")
            )
        end
    end
    
    return p
end

export visualize_quasicrystal_positions, visualize_quasicrystal_tiles
