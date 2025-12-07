"""
Example: Generating and visualizing quasicrystals

This example demonstrates the generation of various quasicrystal structures
using both projection and substitution methods.
"""

ENV["GKSwstype"] = "100"  # For headless plotting

using Lattice2D
using Plots

# Create output directory
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)

println("="^60)
println("Quasicrystal Generation Examples")
println("="^60)

# ============================================================================
# 1. Fibonacci Lattice (1D Quasicrystal)
# ============================================================================
println("\n1. Fibonacci Lattice (1D)")
println("-"^60)

# Projection method
println("  Generating Fibonacci lattice via projection...")
fib_proj = generate_fibonacci_projection(30)
println("    - Generated $(length(fib_proj.positions)) points")

# Substitution method
println("  Generating Fibonacci lattice via substitution (8 generations)...")
fib_subst = generate_fibonacci_substitution(8)
println("    - Generated $(length(fib_subst.positions)) points")
println("    - Sequence length: $(fib_subst.parameters[:sequence_length])")

# Visualize
p1 = visualize_quasicrystal_positions(fib_proj, title="Fibonacci Lattice (Projection)")
savefig(p1, joinpath(output_dir, "fibonacci_projection.png"))

p2 = visualize_quasicrystal_positions(fib_subst, title="Fibonacci Lattice (Substitution)")
savefig(p2, joinpath(output_dir, "fibonacci_substitution.png"))

println("  ✓ Visualizations saved to $output_dir")

# ============================================================================
# 2. Penrose P3 Tiling (2D Quasicrystal with 5-fold symmetry)
# ============================================================================
println("\n2. Penrose P3 Tiling (2D, 5-fold symmetry)")
println("-"^60)

# Projection method
println("  Generating Penrose tiling via projection (radius=4.0)...")
penrose_proj = generate_penrose_projection(4.0)
println("    - Generated $(length(penrose_proj.positions)) vertices")

# Substitution method
println("  Generating Penrose tiling via substitution (3 generations)...")
penrose_subst = generate_penrose_substitution(3)
println("    - Generated $(penrose_subst.parameters[:n_tiles]) tiles")
println("    - Total vertices: $(length(penrose_subst.positions))")

# Visualize
p3 = visualize_quasicrystal_positions(penrose_proj,  
    title="Penrose P3 (Projection)",
    markersize=2)
savefig(p3, joinpath(output_dir, "penrose_projection.png"))

p4 = visualize_quasicrystal_tiles(penrose_subst,
    title="Penrose P3 Tiling (Substitution)")
savefig(p4, joinpath(output_dir, "penrose_substitution.png"))

println("  ✓ Visualizations saved to $output_dir")

# ============================================================================
# 3. Ammann-Beenker Tiling (2D Quasicrystal with 8-fold symmetry)
# ============================================================================
println("\n3. Ammann-Beenker Tiling (2D, 8-fold symmetry)")
println("-"^60)

# Projection method
println("  Generating Ammann-Beenker tiling via projection (radius=4.0)...")
ab_proj = generate_ammann_beenker_projection(4.0)
println("    - Generated $(length(ab_proj.positions)) vertices")

# Substitution method
println("  Generating Ammann-Beenker tiling via substitution (2 generations)...")
ab_subst = generate_ammann_beenker_substitution(2)
println("    - Generated $(ab_subst.parameters[:n_tiles]) tiles")
println("    - Total vertices: $(length(ab_subst.positions))")

# Visualize
p5 = visualize_quasicrystal_positions(ab_proj,
    title="Ammann-Beenker (Projection)",
    markersize=2)
savefig(p5, joinpath(output_dir, "ammann_beenker_projection.png"))

p6 = visualize_quasicrystal_tiles(ab_subst,
    title="Ammann-Beenker Tiling (Substitution)")
savefig(p6, joinpath(output_dir, "ammann_beenker_substitution.png"))

println("  ✓ Visualizations saved to $output_dir")

# ============================================================================
# Summary
# ============================================================================
println("\n" * "="^60)
println("Summary")
println("="^60)
println("Generated quasicrystals:")
println("  • Fibonacci lattice (1D) - 2 methods")
println("  • Penrose P3 tiling (2D, 5-fold) - 2 methods")
println("  • Ammann-Beenker tiling (2D, 8-fold) - 2 methods")
println("\nAll visualizations saved to: $output_dir")
println("="^60)
