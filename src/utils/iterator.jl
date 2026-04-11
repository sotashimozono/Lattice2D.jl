Base.length(lat::PeriodicLattice2D) = lat.N

Base.size(lat::PeriodicLattice2D) = (lat.Lx, lat.Ly)
Base.size(lat::PeriodicLattice2D, d::Int) = d == 1 ? lat.Lx : (d == 2 ? lat.Ly : 1)

Base.eltype(::Type{<:PeriodicLattice2D}) = Int

function Base.iterate(lat::PeriodicLattice2D, state=1)
    if state > lat.N
        return nothing
    else
        return (state, state + 1)
    end
end

function Base.show(io::IO, lat::PeriodicLattice2D)
    Topo = string(typeof(lat.topology))
    bc_summary = join(
        [string(typeof(ax).name.name) for ax in lat.boundary.axes], ", "
    )
    bipartite_str = lat.is_bipartite_flag ? "bipartite" : "not bipartite"
    print(
        io,
        "\n" *
        "Lattice Shape: $Topo\n" *
        "    Lattice Size: $(lat.Lx) x $(lat.Ly)\n" *
        "    total site length : $(lat.N)\n" *
        "    Boundary axes: ($bc_summary)\n" *
        "    Indexing strategy : $(typeof(lat.indexing).name.name)\n" *
        "    Lattice is $(bipartite_str)\n" *
        "Connectivity:\n" *
        "    Total Bonds: $(length(lat.bonds))\n" *
        "Geometry:\n" *
        "    basis matrix (columns are primitive vectors):\n" *
        "      $(lat.basis_matrix)\n",
    )
end
