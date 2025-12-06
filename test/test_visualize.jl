

function plot_lattice(p, lat::Lattice2D)
    # 1. Draw Bonds (PBCの長い線は除外する簡易フィルタ)
    seg_x, seg_y = Float64[], Float64[]
    for b in lat.bonds
        r1, r2 = lat.positions[b.src], lat.positions[b.dst]
        if norm(r1 - r2) < 2.0
            push!(seg_x, r1[1], r2[1], NaN)
            push!(seg_y, r1[2], r2[2], NaN)
        end
    end
    plot!(p, seg_x, seg_y, color=:gray, lw=1.5, label="")
    # 2. Draw Sites (Sublatticeごとに色分け)
    markers = [:circle, :rect, :utriangle]
    colors = [:red, :blue, :green]
    sub_ids = unique(lat.sublattice_ids)
    for (i, sid) in enumerate(sub_ids)
        mask = lat.sublattice_ids .== sid
        xs = [r[1] for r in lat.positions[mask]]
        ys = [r[2] for r in lat.positions[mask]]
        scatter!(p, xs, ys, label="Sub $sid",
            marker=markers[mod1(i, 3)], color=colors[mod1(i, 3)],
            markersize=6, markerstrokewidth=0)
    end
    return p
end

function plot_unit_cell!(p::Plots.Plot, uc::UnitCell{2})
    a1 = uc.basis[1]
    a2 = uc.basis[2]

    # --- 1. 単位胞の枠 (Parallelipiped) ---
    # 頂点: 0 -> a1 -> a1+a2 -> a2 -> 0
    O = [0.0, 0.0]
    corners = [O, a1, a1 .+ a2, a2, O]
    cx = [c[1] for c in corners]
    cy = [c[2] for c in corners]

    # 枠を薄く塗りつぶす
    plot!(p, cx, cy, fill=(0, 0.1, :blue), line=(:dash, :gray), label="Unit Cell Area")

    # --- 2. 基本ベクトル (Basis Vectors) ---
    quiver!(p, [0.0], [0.0], quiver=([a1[1]], [a1[2]]), color=:black, lw=2, label="a1")
    annotate!(p, a1[1], a1[2], text("a1", :bottom, 10))

    quiver!(p, [0.0], [0.0], quiver=([a2[1]], [a2[2]]), color=:black, lw=2, label="a2")
    annotate!(p, a2[1], a2[2], text("a2", :left, 10))

    # --- 3. 結合 (Connections) ---
    for (i, conn) in enumerate(uc.connections)
        r_src = uc.sublattice_positions[conn.src_sub]

        # 接続先の座標計算: 相対位置 + シフト(dx*a1 + dy*a2)
        shift = conn.dx .* a1 .+ conn.dy .* a2
        r_dst = uc.sublattice_positions[conn.dst_sub] .+ shift

        # 単位胞内(dx=0, dy=0)か、外への接続かでスタイルを変える
        is_internal = (conn.dx == 0 && conn.dy == 0)

        l_style = is_internal ? :solid : :dot
        l_alpha = is_internal ? 1.0 : 0.6
        l_width = is_internal ? 2.0 : 1.5

        # 線を引く
        plot!(p, [r_src[1], r_dst[1]], [r_src[2], r_dst[2]],
            line=(l_style, l_width, :gray), alpha=l_alpha, label="")

        # もし外への接続なら、接続先の「ゴースト原子」を描画して分かりやすくする
        if !is_internal
            scatter!(p, [r_dst[1]], [r_dst[2]],
                marker=:circle, markersize=4, color=:gray, alpha=0.5, label="")
        end
    end
    # --- 4. 原子 (Sites) ---
    colors = [:red, :blue, :green, :orange] # sublattice 毎の色
    markers = [:circle, :square, :utriangle, :diamond]

    for (i, pos) in enumerate(uc.sublattice_positions)
        scatter!(p, [pos[1]], [pos[2]],
            marker=markers[mod1(i, 4)],
            color=colors[mod1(i, 4)],
            markersize=8,
            markerstrokecolor=:black,
            label="Sub $i")
    end

    # 見た目の調整
    plot!(p, aspect_ratio=:equal, grid=false, frame=:box)
    return p
end

function plot_unit_cell(uc::UnitCell; kwargs...)
    p = plot(; kwargs...)
    plot_unit_cell!(p, uc)
    return p
end
@testset "Lattice Construction" begin
    Lx, Ly = 4, 4
    @testset "Square" begin
        lat = build_lattice(Square, Lx, Ly)
        @test lat.N == Lx * Ly
        @test lat.is_bipartite == true
        # 配位数 4
        @test all(length(n) == 4 for n in lat.nearest_neighbors)

        p1 = plot(title="Square Lattice", aspect_ratio=:equal, grid=false, axis=false, legend=:outerright)
        plot_lattice(p1, lat)
        p2 = plot_unit_cell(lat.unit_cell)
        p = plot(p1, p2, layout=(1, 2), size=(800, 400))
        savefig(p, joinpath(PATHS[:geometry], "square_lattice.png"))
    end
    @testset "Triangular" begin
        lat = build_lattice(Triangular, Lx, Ly)
        @test lat.N == Lx * Ly
        @test lat.is_bipartite == false
        # 配位数 6
        @test all(length(n) == 6 for n in lat.nearest_neighbors)
        p1 = plot(title="Triangular Lattice", aspect_ratio=:equal, grid=false, axis=false, legend=:outerright)
        plot_lattice(p1, lat)
        p2 = plot_unit_cell(lat.unit_cell)
        p = plot(p1, p2, layout=(1, 2), size=(800, 400))
        savefig(p, joinpath(PATHS[:geometry], "triangular_lattice.png"))
    end
    @testset "Honeycomb" begin
        lat = build_lattice(Honeycomb, Lx, Ly)
        @test lat.N == 2 * Lx * Ly
        @test lat.is_bipartite == true
        # 配位数 3
        @test all(length(n) == 3 for n in lat.nearest_neighbors)
        p1 = plot(title="Honeycomb Lattice", aspect_ratio=:equal, grid=false, axis=false, legend=:outerright)
        plot_lattice(p1, lat)
        p2 = plot_unit_cell(lat.unit_cell)
        p = plot(p1, p2, layout=(1, 2), size=(800, 400))
        savefig(p, joinpath(PATHS[:geometry], "honeycomb_lattice.png"))
    end
    @testset "Kagome" begin
        lat = build_lattice(Kagome, Lx, Ly)
        @test lat.N == 3 * Lx * Ly
        @test lat.is_bipartite == false
        # 配位数 4
        @test all(length(n) == 4 for n in lat.nearest_neighbors)
        p1 = plot(title="Kagome Lattice", aspect_ratio=:equal, grid=false, axis=false, legend=:outerright)
        plot_lattice(p1, lat)
        p2 = plot_unit_cell(lat.unit_cell)
        p = plot(p1, p2, layout=(1, 2), size=(800, 400))
        savefig(p, joinpath(PATHS[:geometry], "kagome_lattice.png"))
    end
    @testset "Lieb" begin
        lat = build_lattice(Lieb, Lx, Ly)
        @test lat.N == 3 * Lx * Ly
        @test lat.is_bipartite == true
        # @test all(length(n) == 4 for n in lat.nearest_neighbors)
        p1 = plot(title="Lieb Lattice", aspect_ratio=:equal, grid=false, axis=false, legend=:outerright)
        plot_lattice(p1, lat)
        p2 = plot_unit_cell(lat.unit_cell)
        p = plot(p1, p2, layout=(1, 2), size=(800, 400))
        savefig(p, joinpath(PATHS[:geometry], "lieb_lattice.png"))
    end
end