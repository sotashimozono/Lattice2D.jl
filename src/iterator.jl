Base.length(lat::Lattice2D{Topology}) where {Topology} = lat.N

Base.size(lat::Lattice2D{Topology}) where {Topology} = (lat.Lx, lat.Ly)
Base.size(lat::Lattice2D, d::Int) = d == 1 ? lat.Lx : (d == 2 ? lat.Ly : 1)

Base.size(lat::Lattice2D) = (lat.Lx, lat.Ly)
Base.size(lat::Lattice2D, d::Int) = d == 1 ? lat.Lx : (d == 2 ? lat.Ly : 1)
Base.eltype(::Type{<:Lattice2D}) = Int

function Base.iterate(lat::Lattice2D, state=1)
    if state > lat.N
        return nothing
    else
        return (state, state + 1)
    end
end
