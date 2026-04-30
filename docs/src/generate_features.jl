using Lattice2D
using Plots
using LinearAlgebra

out_dir = joinpath(@__DIR__, "assets", "figures", "features")
mkpath(out_dir)

# 1. Site specification
l = build_lattice(Square, 4, 4)
p = plot(; aspect_ratio=:equal, legend=false, grid=false, axis=false, size=(400, 400))
# Plot all sites empty
for i in 1:num_sites(l)
    pos = position(l, i)
    scatter!(
        p,
        [pos[1]],
        [pos[2]];
        color=:lightgray,
        markersize=10,
        markerstrokewidth=2,
        markercolor=:white,
    )
end
# Highlight site 6 and its neighbors
target = 6
pos = position(l, target)
scatter!(p, [pos[1]], [pos[2]]; color=:red, markersize=12)
for n in neighbors(l, target)
    npos = position(l, n)
    plot!(p, [pos[1], npos[1]], [pos[2], npos[2]]; color=:red, linewidth=3)
    scatter!(p, [npos[1]], [npos[2]]; color=:orange, markersize=10)
end
savefig(p, joinpath(out_dir, "site_neighbors.png"))

# 2. Bond specification
p2 = plot(; aspect_ratio=:equal, legend=false, grid=false, axis=false, size=(400, 400))
for b in bonds(l)
    p1 = position(l, b.i)
    p2_pos = position(l, b.j)
    plot!(p2, [p1[1], p2_pos[1]], [p1[2], p2_pos[2]]; color=:lightgray, linewidth=2)
end
for i in 1:num_sites(l)
    pos = position(l, i)
    scatter!(p2, [pos[1]], [pos[2]]; color=:darkgray, markersize=6)
end
# Highlight a specific bond
target_b = bonds(l)[10]
p1 = position(l, target_b.i)
p2_pos = position(l, target_b.j)
plot!(p2, [p1[1], p2_pos[1]], [p1[2], p2_pos[2]]; color=:blue, linewidth=6)
center = bond_center(l, target_b)
scatter!(p2, [center[1]], [center[2]]; color=:cyan, markersize=8, marker=:square)
savefig(p2, joinpath(out_dir, "bond_selection.png"))

# 3. Plaquette specification
l_kagome = build_lattice(Kagome, 3, 3)
p3 = plot(; aspect_ratio=:equal, legend=false, grid=false, axis=false, size=(400, 400))
for b in bonds(l_kagome)
    p1 = position(l_kagome, b.i)
    p2_pos = position(l_kagome, b.j)
    plot!(p3, [p1[1], p2_pos[1]], [p1[2], p2_pos[2]]; color=:lightgray, linewidth=1)
end
# highlight a plaquette (hexagon)
plaqs = plaquettes(l_kagome)
target_p = plaqs[3] # pick a hexagon
v_positions = [position(l_kagome, i) for i in target_p.vertices]
push!(v_positions, position(l_kagome, target_p.vertices[1])) # close it
px = [v[1] for v in v_positions]
py = [v[2] for v in v_positions]
plot!(
    p3, px, py; fill=true, color=:lightgreen, fillalpha=0.5, linewidth=2, linecolor=:green
)
center = target_p.center
scatter!(p3, [center[1]], [center[2]]; color=:darkgreen, markersize=8, marker=:star5)
savefig(p3, joinpath(out_dir, "plaquette_selection.png"))

println("Features generated!")
