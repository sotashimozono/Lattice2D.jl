ENV["GKSwstype"] = "100"

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Lattice2D
using Plots
using LinearAlgebra

const OUT_DIR = joinpath(pkgdir(Lattice2D), "docs", "src", "assets", "figures", "lattice")
mkpath(OUT_DIR)

function plot_lattice_structure(lat::Lattice, title_str)
    p = plot(title=title_str, aspect_ratio=:equal, grid=false, axis=false, legend=:outertopright)
    
    threshold = 1.5 * maximum(norm.(lat.unit_cell.basis)) 
    
    seg_x, seg_y = Float64[], Float64[]
    for b in lat.bonds
        r1, r2 = lat.positions[b.src], lat.positions[b.dst]
        if norm(r1 - r2) < threshold
            push!(seg_x, r1[1], r2[1], NaN)
            push!(seg_y, r1[2], r2[2], NaN)
        end
    end
    plot!(p, seg_x, seg_y, color=:gray, lw=1.5, label="")

    markers = [:circle, :rect, :utriangle, :diamond]
    colors = [:firebrick, :royalblue, :forestgreen, :darkorange]
    
    sub_ids = sort(unique(lat.sublattice_ids))
    for (i, sid) in enumerate(sub_ids)
        mask = lat.sublattice_ids .== sid
        xs = [r[1] for r in lat.positions[mask]]
        ys = [r[2] for r in lat.positions[mask]]
        
        scatter!(p, xs, ys, label="Sub $sid",
            marker=markers[mod1(i, 4)], 
            color=colors[mod1(i, 4)],
            markersize=5, markerstrokewidth=0)
    end
    return p
end

function plot_unit_cell_geometry(uc::UnitCell{2})
    p = plot(aspect_ratio=:equal, grid=false, frame=:box, title="Unit Cell", legend=:outertopright)
    a1, a2 = uc.basis[1], uc.basis[2]
    O = [0.0, 0.0]
    corners = [O, a1, a1 .+ a2, a2, O]
    plot!(p, [c[1] for c in corners], [c[2] for c in corners], 
          fill=(0, 0.05, :blue), line=(:dash, :gray), label="unit cell area")

    # Basis Vectors
    quiver!(p, [0.0], [0.0], quiver=([a1[1]], [a1[2]]), color=:black, lw=2, label="a1")
    annotate!(p, a1[1]/2, a1[2]/2 - 0.1, text("a1", 10))
    quiver!(p, [0.0], [0.0], quiver=([a2[1]], [a2[2]]), color=:black, lw=2, label="a2")
    annotate!(p, a2[1]/2 - 0.1, a2[2]/2, text("a2", 10))


    # Connections
    for conn in uc.connections
        r_src = uc.sublattice_positions[conn.src_sub]
        shift = conn.dx .* a1 .+ conn.dy .* a2
        r_dst = uc.sublattice_positions[conn.dst_sub] .+ shift
        
        is_internal = (conn.dx == 0 && conn.dy == 0)
        l_style = is_internal ? :solid : :dot
        l_alpha = is_internal ? 1.0 : 0.5
        
        plot!(p, [r_src[1], r_dst[1]], [r_src[2], r_dst[2]],
              line=(l_style, 2.0, :gray), alpha=l_alpha, label="")
        
        if !is_internal
            scatter!(p, [r_dst[1]], [r_dst[2]], 
                    marker=:circle, ms=3, color=:gray, alpha=0.5, label="")
        end
    end
    plot!(p, [NaN], [NaN], line=(:solid, 2.0, :gray), label="Intra-cell Bond")
    plot!(p, [NaN], [NaN], line=(:dot, 1.5, :gray), alpha=0.6, label="Inter-cell Bond")
    scatter!(p, [NaN], [NaN], marker=:circle, ms=4, color=:gray, alpha=0.5, label="Image Site")
    markers = [:circle, :rect, :utriangle, :diamond]
    colors = [:firebrick, :royalblue, :forestgreen, :darkorange]

    for (i, pos) in enumerate(uc.sublattice_positions)
        scatter!(p, [pos[1]], [pos[2]],
            marker=markers[mod1(i, 4)],
            color=colors[mod1(i, 4)],
            markersize=8, markerstrokecolor=:black,
            label="Sub $i")
    end
    
    return p
end

lattices_to_plot = [
    (Square, "Square"),
    (Triangular, "Triangular"),
    (Honeycomb, "Honeycomb"),
    (Kagome, "Kagome"),
    (Lieb, "Lieb"),
    (ShastrySutherland, "Shastry-Sutherland")
]

Lx, Ly = 4, 4 

for (T, name) in lattices_to_plot
    println("Generating plot for: $name")
    lat = build_lattice(T, Lx, Ly)
    
    p1 = plot_lattice_structure(lat, name)
    p2 = plot_unit_cell_geometry(lat.unit_cell)
    
    final_plot = plot(p1, p2, layout=(1, 2), size=(900, 400), margin=5Plots.mm)
    
    fname = lowercase(string(T)) * "_lattice.png"
    save_path = joinpath(OUT_DIR, fname)
    
    savefig(final_plot, save_path)
    println("Saved to $save_path")
end