Base.length(lat::Lattice) = num_sites(lat)

Base.size(lat::Lattice) = (lat.Lx, lat.Ly)
Base.size(lat::Lattice, d::Int) = d == 1 ? lat.Lx : (d == 2 ? lat.Ly : 1)

Base.eltype(::Type{<:Lattice}) = Int

function Base.iterate(lat::Lattice, state=1)
    if state > num_sites(lat)
        return nothing
    else
        return (state, state + 1)
    end
end

function Base.show(io::IO, lat::Lattice{Topo}) where {Topo}
    bc_summary = join([string(typeof(ax).name.name) for ax in lat.boundary.axes], ", ")
    bipartite_str = is_bipartite(lat) ? "bipartite" : "not bipartite"
    print(
        io,
        "\n" *
        "Lattice Shape: $(Topo)\n" *
        "    Lattice Size: $(lat.Lx) x $(lat.Ly)\n" *
        "    total site length : $(num_sites(lat))\n" *
        "    Boundary axes: ($bc_summary)\n" *
        "    Indexing strategy : $(typeof(lat.indexing).name.name)\n" *
        "    Lattice is $(bipartite_str)\n" *
        "Geometry:\n" *
        "    basis matrix (columns are primitive vectors):\n" *
        "      $(basis_vectors(lat))\n",
    )
end
