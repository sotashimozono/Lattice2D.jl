# LatticeCore API

`Lattice2D` builds on [LatticeCore.jl](https://github.com/sotashimozono/LatticeCore.jl),
which defines the abstract lattice vocabulary (the `AbstractLattice` hierarchy,
`Bond`, boundary-condition singletons, indexing strategies, site types, …). These
symbols are re-exported by `Lattice2D`, so `using Lattice2D` alone is enough to
reach them.

They are documented here so that the `@ref` cross-references inside `Lattice2D`
docstrings resolve to a real target (see issue #68).

```@autodocs
Modules = [LatticeCore]
# Skip the plotting/diffraction stubs: their docstrings reference the
# `LatticeCorePlotsExt` package-extension module, which is not documented in
# this doc set (the concrete methods live in the extension). Lattice2D does
# not re-export these anyway — its own visualisation lives in Lattice2DPlotsExt.
Filter = obj -> try
    nameof(obj) ∉ (:plot_lattice, :plot_bonds!, :plot_sites!, :diffraction_pattern)
catch
    true
end
```
