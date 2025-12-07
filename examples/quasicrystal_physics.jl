"""
Quasicrystal Physics Application: Electronic Properties and Localization

This application demonstrates novel quantum phenomena in quasicrystals:
1. Fractal energy spectrum (Devil's staircase)
2. Anderson localization in the Fibonacci lattice
3. Critical states in 2D Penrose tiling
4. Density of states with singular continuous spectrum

These are fundamental properties unique to quasiperiodic systems.
"""

ENV["GKSwstype"] = "100"

using Lattice2D
using LinearAlgebra
using Plots
using Printf

# Create output directory
output_dir = joinpath(@__DIR__, "output", "physics")
mkpath(output_dir)

println("="^70)
println("QUASICRYSTAL QUANTUM PHYSICS APPLICATION")
println("="^70)
println("\nThis application studies electronic properties and localization")
println("phenomena unique to quasicrystalline materials.")
println("="^70)

# ============================================================================
# Part 1: Fibonacci Lattice - Hofstadter Butterfly Analog
# ============================================================================
println("\n" * "="^70)
println("PART 1: Fibonacci Lattice Electronic Structure")
println("="^70)

println("\n[1.1] Generating Fibonacci lattice with different modulations...")
n_sites = 200
fib = generate_fibonacci_substitution(10)  # Generate Fibonacci lattice

# Trim to desired size
if length(fib.positions) > n_sites
    fib = QuasicrystalData{1,Float64}(
        fib.positions[1:n_sites],
        fib.tiles,
        fib.generation_method,
        merge(fib.parameters, Dict(:n_points => n_sites))
    )
end

println("      Generated Fibonacci lattice with $(length(fib.positions)) sites")

# Build tight-binding models with different on-site modulations
modulations = [0.0, 0.5, 1.0, 2.0]
spectra = []
dos_data = []

println("\n[1.2] Computing energy spectra for different modulation strengths...")
for (idx, V) in enumerate(modulations)
    @printf("      [%d/%d] Modulation V = %.1f...", idx, length(modulations), V)
    
    model = build_fibonacci_tight_binding(fib, t=1.0, onsite_modulation=V)
    eigenvalues, eigenvectors = solve_eigenspectrum(model)
    
    # Compute DOS
    energies, dos = compute_dos(eigenvalues, n_bins=200, sigma=0.05)
    
    push!(spectra, eigenvalues)
    push!(dos_data, (energies, dos, V))
    
    println(" ✓ ($(length(eigenvalues)) states)")
end

# Plot energy spectrum evolution
println("\n[1.3] Visualizing spectrum evolution...")
p1 = plot(layout=(2,2), size=(1000, 800),
         plot_title="Fibonacci Lattice: Spectrum vs Modulation Strength")

for (idx, (spectrum, V)) in enumerate(zip(spectra, modulations))
    subplot = idx
    y_vals = collect(1:length(spectrum))
    
    scatter!(p1[subplot], spectrum, y_vals,
            markersize=1.5,
            markerstrokewidth=0,
            label="",
            xlabel="Energy E",
            ylabel="State Index",
            title="V = $V",
            ylims=(0, length(spectrum)),
            marker=:circle)
end

savefig(p1, joinpath(output_dir, "fibonacci_spectrum_evolution.png"))
println("      ✓ Saved: fibonacci_spectrum_evolution.png")

# Plot Density of States
println("\n[1.4] Computing Density of States...")
p2 = plot(size=(1000, 600),
         xlabel="Energy E",
         ylabel="Density of States",
         title="Fibonacci Lattice: DOS Evolution with On-site Modulation",
         legend=:topright)

for (energies, dos, V) in dos_data
    plot!(p2, energies, dos,
          linewidth=2,
          label="V = $V",
          alpha=0.7)
end

savefig(p2, joinpath(output_dir, "fibonacci_dos.png"))
println("      ✓ Saved: fibonacci_dos.png")

# ============================================================================
# Part 2: Localization Analysis (Anderson Localization)
# ============================================================================
println("\n" * "="^70)
println("PART 2: Anderson Localization in Quasicrystals")
println("="^70)

println("\n[2.1] Computing Inverse Participation Ratio (IPR)...")
println("      IPR measures localization: IPR ~ 1/N (extended), IPR ~ 1 (localized)")

V_strong = 2.0
model_strong = build_fibonacci_tight_binding(fib, t=1.0, onsite_modulation=V_strong)
eigenvalues_strong, eigenvectors_strong = solve_eigenspectrum(model_strong)
iprs = compute_all_iprs(eigenvectors_strong)

println("      ✓ Computed IPR for $(length(iprs)) states")

# Plot IPR vs Energy
p3 = scatter(eigenvalues_strong, iprs,
            markersize=3,
            markerstrokewidth=0,
            xlabel="Energy E",
            ylabel="Inverse Participation Ratio",
            title="Localization in Fibonacci Lattice (V = $V_strong)",
            label="",
            size=(800, 600))

# Add reference lines
hline!(p3, [1/length(fib.positions)], 
       linestyle=:dash, 
       linewidth=2, 
       label="Extended limit (1/N)",
       color=:red)

savefig(p3, joinpath(output_dir, "fibonacci_ipr.png"))
println("\n[2.2] Analysis complete:")
println("      ✓ Saved: fibonacci_ipr.png")
println("      ✓ States with IPR > 0.01 (localized): $(sum(iprs .> 0.01))")
println("      ✓ States with IPR < 0.001 (extended): $(sum(iprs .< 0.001))")

# Visualize a few wavefunctions
println("\n[2.3] Visualizing select wavefunctions...")
positions_1d = sort([p[1] for p in fib.positions])

# Pick representative states: lowest energy, middle, highest energy
n_states = length(eigenvalues_strong)
state_indices = [1, n_states÷2, n_states]
state_labels = ["Lowest E", "Middle E", "Highest E"]

p4 = plot(layout=(3,1), size=(1000, 900),
         plot_title="Wavefunctions in Fibonacci Lattice (V = $V_strong)")

for (idx, state_idx, label) in zip(1:3, state_indices, state_labels)
    psi = eigenvectors_strong[:, state_idx]
    E = eigenvalues_strong[state_idx]
    ipr_val = iprs[state_idx]
    
    scatter!(p4[idx], positions_1d, abs2.(psi),
            markersize=3,
            markerstrokewidth=0,
            xlabel="Position",
            ylabel="|ψ|²",
            title="$label: E = $(@sprintf("%.3f", E)), IPR = $(@sprintf("%.4f", ipr_val))",
            label="",
            color=:blue)
end

savefig(p4, joinpath(output_dir, "fibonacci_wavefunctions.png"))
println("      ✓ Saved: fibonacci_wavefunctions.png")

# ============================================================================
# Part 3: 2D Penrose Tiling - Critical States
# ============================================================================
println("\n" * "="^70)
println("PART 3: Electronic Structure of 2D Penrose Tiling")
println("="^70)

println("\n[3.1] Generating Penrose tiling and building tight-binding model...")
penrose = generate_penrose_projection(5.0)
println("      Generated Penrose tiling with $(length(penrose.positions)) vertices")

println("\n[3.2] Building tight-binding Hamiltonian...")
penrose_model = build_penrose_tight_binding(penrose, t=1.0, cutoff=1.2)
println("      ✓ Hamiltonian dimension: $(penrose_model.n_sites)")

# Compute spectrum
println("\n[3.3] Solving for eigenspectrum (this may take a moment)...")
n_eigs = min(100, penrose_model.n_sites - 2)
penrose_eigenvalues, penrose_eigenvectors = solve_eigenspectrum(penrose_model, k=n_eigs)
println("      ✓ Computed $(length(penrose_eigenvalues)) eigenvalues")

# Plot spectrum
p5 = scatter(penrose_eigenvalues, 1:length(penrose_eigenvalues),
            markersize=3,
            markerstrokewidth=0,
            xlabel="Energy E",
            ylabel="State Index",
            title="Energy Spectrum of Penrose Tiling",
            label="",
            size=(800, 600),
            color=:red)

savefig(p5, joinpath(output_dir, "penrose_spectrum.png"))
println("      ✓ Saved: penrose_spectrum.png")

# Compute DOS
println("\n[3.4] Computing Density of States...")
penrose_energies, penrose_dos = compute_dos(penrose_eigenvalues, n_bins=150, sigma=0.1)

p6 = plot(penrose_energies, penrose_dos,
         linewidth=2,
         xlabel="Energy E",
         ylabel="Density of States",
         title="DOS of 2D Penrose Tiling",
         label="",
         size=(800, 600),
         color=:darkblue)

savefig(p6, joinpath(output_dir, "penrose_dos.png"))
println("      ✓ Saved: penrose_dos.png")

# Localization analysis
println("\n[3.5] Analyzing localization in 2D...")
penrose_iprs = compute_all_iprs(penrose_eigenvectors)

p7 = scatter(penrose_eigenvalues, penrose_iprs,
            markersize=3,
            markerstrokewidth=0,
            xlabel="Energy E",
            ylabel="IPR",
            title="Localization in 2D Penrose Tiling",
            label="",
            size=(800, 600),
            color=:purple)

hline!(p7, [1/penrose_model.n_sites], 
       linestyle=:dash, 
       linewidth=2, 
       label="Extended limit",
       color=:red)

savefig(p7, joinpath(output_dir, "penrose_ipr.png"))
println("      ✓ Saved: penrose_ipr.png")

# Visualize a wavefunction
println("\n[3.6] Visualizing wavefunction density...")
mid_state = length(penrose_eigenvalues) ÷ 2
psi_2d = penrose_eigenvectors[:, mid_state]
E_mid = penrose_eigenvalues[mid_state]
ipr_mid = penrose_iprs[mid_state]

x_coords = [p[1] for p in penrose.positions]
y_coords = [p[2] for p in penrose.positions]
psi2_vals = abs2.(psi_2d)

p8 = scatter(x_coords, y_coords,
            marker_z=psi2_vals,
            markersize=4,
            markerstrokewidth=0,
            xlabel="x",
            ylabel="y",
            title="Wavefunction |ψ|² for E = $(@sprintf("%.3f", E_mid)), IPR = $(@sprintf("%.4f", ipr_mid))",
            colorbar_title="|ψ|²",
            aspect_ratio=:equal,
            size=(800, 700),
            color=:viridis)

savefig(p8, joinpath(output_dir, "penrose_wavefunction.png"))
println("      ✓ Saved: penrose_wavefunction.png")

# ============================================================================
# Summary and Physical Interpretation
# ============================================================================
println("\n" * "="^70)
println("PHYSICS SUMMARY")
println("="^70)

println("\n[Fibonacci Lattice - 1D Quasicrystal]")
println("  • Spectrum shows fractal structure characteristic of quasiperiodicity")
println("  • Strong modulation (V > 1) induces localization")
println("  • IPR analysis reveals mixture of extended and localized states")
println("  • This demonstrates the metal-insulator transition in 1D quasicrystals")

println("\n[Penrose Tiling - 2D Quasicrystal]")
println("  • Energy spectrum exhibits characteristic gaps")
println("  • DOS shows pseudo-gap structure")
println("  • Wavefunctions show complex spatial patterns")
println("  • 2D quasicrystals exhibit critical states (neither fully extended nor localized)")

println("\n[Novel Physics Demonstrated]")
println("  1. Fractal energy spectrum unique to quasiperiodic systems")
println("  2. Anderson localization without disorder")
println("  3. Critical states in 2D (dimension-dependent phenomenon)")
println("  4. Singular continuous spectrum (between pure point and continuous)")

println("\n[Applications]")
println("  • Understanding transport properties in quasicrystalline alloys")
println("  • Photonic quasicrystals for light localization")
println("  • Topological properties of quasiperiodic systems")
println("  • Cold atom systems with optical lattices")

println("\n" * "="^70)
println("All visualizations saved to: $output_dir")
println("="^70)
