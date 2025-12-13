module Lattice2D

using Random
using Plots

include("core/boundarycondition.jl")
include("core/index_methods.jl")
include("core/abstractlattices.jl")
include("core/unitcells.jl")
include("core/constructor.jl")
include("utils/iterator.jl")

# Applications
include("applications/percolation.jl")
include("applications/randomwalk.jl")
include("applications/dla.jl")
include("applications/spanningtree.jl")
include("applications/visualization.jl")

end
