# Quasicrystal Physics Application

This module adds quantum physics calculations for quasicrystals using tight-binding models.

## Overview

The tight-binding approximation is a standard method in condensed matter physics for studying electronic properties of materials. For quasicrystals, this reveals unique phenomena:

1. **Fractal Energy Spectrum**: Unlike periodic crystals, quasicrystals exhibit complex, hierarchical energy spectra
2. **Anderson Localization**: Wave functions can become localized even without disorder
3. **Critical States**: In 2D, states are neither fully extended nor fully localized
4. **Singular Continuous Spectrum**: A unique feature of quasiperiodic systems

## Features

### Tight-Binding Models
```julia
using Lattice2D

# Build model for Fibonacci lattice
fib = generate_fibonacci_substitution(8)
model = build_fibonacci_tight_binding(fib, t=1.0, onsite_modulation=1.0)

# Solve for energy spectrum
eigenvalues, eigenvectors = solve_eigenspectrum(model)
```

### Localization Analysis
```julia
# Compute Inverse Participation Ratio (IPR)
iprs = compute_all_iprs(eigenvectors)

# IPR ~ 1/N: extended state
# IPR ~ 1: localized state
```

### Density of States
```julia
energies, dos = compute_dos(eigenvalues, n_bins=100, sigma=0.1)
```

## Physical Significance

### Fibonacci Lattice (1D)
- Demonstrates metal-insulator transition
- On-site modulation induces localization
- Shows fractal spectrum characteristic of quasiperiodicity

### Penrose Tiling (2D)  
- Exhibits pseudo-gap in DOS
- Critical states (multifractal wavefunctions)
- Unique to 2D quasicrystals

## Applications

1. **Quasicrystalline Alloys**: Understanding electrical conductivity
2. **Photonic Quasicrystals**: Light localization for optical devices
3. **Cold Atoms**: Optical lattice experiments
4. **Topological Physics**: Edge states and topological invariants

## Example Usage

See `examples/quasicrystal_physics.jl` for a comprehensive demonstration including:
- Spectrum evolution with modulation strength
- IPR analysis and localization visualization
- Wavefunction density plots
- 2D Penrose tiling electronic structure

## API

### Core Types
- `TightBindingModel{T}`: Hamiltonian and site positions
  
### Model Building
- `build_fibonacci_tight_binding(qc_data; t, onsite_modulation)`
- `build_penrose_tight_binding(qc_data; t, cutoff)`

### Analysis
- `solve_eigenspectrum(model; k)`: Eigenvalues and eigenvectors
- `compute_dos(eigenvalues; n_bins, sigma)`: Density of states
- `compute_ipr(eigenvector)`: Inverse Participation Ratio
- `compute_all_iprs(eigenvectors)`: IPR for all states

## References

1. Kohmoto et al., PRL (1983): "Localization Problem in One Dimension"
2. Sutherland, PRB (1986): "Self-similar ground-state wave function for electrons on a two-dimensional Penrose lattice"
3. Janot, "Quasicrystals: A Primer" (1992)
4. Kraus et al., PRL (2012): "Topological States and Adiabatic Pumping in Quasicrystals"
