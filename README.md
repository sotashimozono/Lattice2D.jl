# Lattice2D.jl

[![docs: dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sotashimozono.github.io/Lattice2D.jl/dev/)
[![Julia](https://img.shields.io/badge/julia-v1.10+-9558b2.svg)](https://julialang.org)
[![Code Style: Blue](https://img.shields.io/badge/Code%20Style-Blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![codecov](https://codecov.io/gh/sotashimozono/Lattice2D.jl/graph/badge.svg?token=6E7VZ9MJMK)](https://codecov.io/gh/sotashimozono/Lattice2D.jl)
[![Build Status](https://github.com/sotashimozono/Lattice2D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sotashimozono/Lattice2D.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Features

This package provides two dimensional lattices. If you know unit cell information (basic vectors, connection), you can construct arbitrary lattices.

Available lattices by defaults are:

- Square lattice
- Triangular lattice
- Honeycomb lattice
- Kagome lattice
- Lieb lattice
- Shastry-Sutherland lattice

### Application Examples

The package now includes 5 comprehensive application examples:

1. **Percolation** - Bond and site percolation with cluster analysis
2. **Random Walks** - Standard and self-avoiding random walks with statistics
3. **Diffusion Limited Aggregation (DLA)** - Fractal growth simulations
4. **Spanning Trees** - Uniform random spanning tree generation (Wilson's algorithm)
5. **Visualization** - Enhanced plotting functions for all applications

See [Applications Documentation](docs/src/applications.md) for detailed usage.

## Installation

```julia
using Pkg
Pkg.add("Lattice2D")
```

## Example

Here we show Honeycomb lattice as an example. The Module enables us to have unit cell informations and connection infos, reciprocal vectors etc. And more, we can know the graph is bipartite or not, we can set periodic boundary and open boundary.

```julia
using Lattice2D

# Create a 4x4 Honeycomb lattice with Periodic Boundary Conditions (PBC)
lat = build_lattice(Honeycomb, 4, 4; boundary=PBC())

# Check basic properties
println("Total sites: ", lat.N)
println("Is bipartite? ", lat.is_bipartite)
```

![Honeycomb Lattice](docs/src/assets/figures/lattice/honeycomb_lattice.png)

## Contributing

Contributions are welcome! Please feel free to open an Issue or submit a pull requests.
