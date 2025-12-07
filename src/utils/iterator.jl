Base.length(lat::Lattice{Topology}) where {Topology} = lat.N

Base.size(lat::Lattice{Topology}) where {Topology} = (lat.Lx, lat.Ly)
Base.size(lat::Lattice, d::Int) = d == 1 ? lat.Lx : (d == 2 ? lat.Ly : 1)

Base.eltype(::Type{<:Lattice}) = Int

function Base.iterate(lat::Lattice, state=1)
    if state > lat.N
        return nothing
    else
        return (state, state + 1)
    end
end

function Base.getproperty(lat::Lattice{Topology, Type, Boundary}, sym::Symbol) where {Topology, Type, Boundary}
    if sym === :unit_cell
        return get_unit_cell(Topology)
    elseif sym === :Topology || sym === :topology
        return Topology
    elseif sym === :Boundary || sym === :boundary
        return Boundary
    else
        return getfield(lat, sym)
    end
end

function Base.propertynames(lat::Lattice)
    return (fieldnames(Lattice)..., :unit_cell)
end

function Base.show(io::IO, lat::Lattice)
    T_str = string(lat.topology)
    B_str = string(lat.boundary)
    str = lat.is_bipartite ? "bipartite" : "not bipartite"
    # 単一の print(io, ...) で全体を出力
    print(io, 
        "\n" *
        "Lattice Shape: $T_str\n" *
        "    Lattice Size: $(lat.Lx) x $(lat.Ly)\n" *
        "    total site length : $(lat.N)\n" *
        "    Boundary Condition: $B_str\n" *
        "    Indexing strategy : $(lat.index_method)\n" *
        "    Lattice is $(str)\n" *
        "Connectivity:\n" *
        "    Total Bonds: $(length(lat.bonds))\n" *
        "    Is Bipartite: $(lat.is_bipartite)\n" *
        "Geometry:\n" *
        "    basis vector: $(lat.basis_vectors)\n" *
        "    reciprocal vector: $(lat.reciprocal_vectors)\n"
    )
end