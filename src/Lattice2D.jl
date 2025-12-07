module Lattice2D

using LinearAlgebra

include("core/boundarycondition.jl")
include("core/index_methods.jl")
include("core/abstractlattices.jl")
include("core/quasicrystals.jl")
include("core/unitcells.jl")
include("core/constructor.jl")
include("utils/iterator.jl")

# Quasicrystal implementations
include("quasicrystals/fibonacci.jl")
include("quasicrystals/penrose.jl")
include("quasicrystals/ammann_beenker.jl")

end
