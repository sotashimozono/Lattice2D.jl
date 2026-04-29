@testset "_connection_steps type stability" begin
    # `_connection_steps` is on the hot path of `neighbors` /
    # `neighbor_bonds`. Issue #45: previously its eltype was an
    # abstract `Tuple{Int,SVector{2,T},Symbol}` literal and it
    # built bond-type symbols dynamically per iteration via
    # `Symbol("type_", conn.type)`. The new implementation returns a
    # `Vector{Lattice2D.Step{T}}` and looks up symbols from a cached
    # per-topology table, so the return type must be inferrable.
    for Topo in (
        Lattice2D.Square,
        Lattice2D.Triangular,
        Lattice2D.Honeycomb,
        Lattice2D.Kagome,
        Lattice2D.Lieb,
        Lattice2D.ShastrySutherland,
        Lattice2D.Dice,
        Lattice2D.UnionJack,
    )
        lat = build_lattice(Topo, 3, 3)
        steps = @inferred Lattice2D._connection_steps(lat, 1)
        @test eltype(steps) === Lattice2D.Step{Float64}
        # Every emitted symbol matches the cached table value.
        syms = Lattice2D._conn_type_symbols(Topo)
        for st in steps
            @test st.type in syms
        end

        # Downstream `neighbor_bonds` should still be inferrable.
        nb = @inferred neighbor_bonds(lat, 1)
        @test !isempty(nb)
    end
end
