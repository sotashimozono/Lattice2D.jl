"""
    Lattice2D convenience builders
    ==============================

This file defines a thin layer of one-liner builders that delegate to
[`build_lattice`](@ref). They exist purely to make call sites read more
naturally:

```julia
lat = honeycomb(6, 6; boundary = OpenAxis())
# is exactly equivalent to
lat = build_lattice(Honeycomb, 6, 6; boundary = OpenAxis())
```

Every builder accepts the same keyword arguments as `build_lattice`
(`boundary`, `indexing`, `layout`) and returns the same `Lattice` value
— no behaviour is added or hidden.

The second linear extent defaults to the first, so `square(8)` builds
an `8 × 8` sample.
"""

"""
    square(L::Int, M::Int = L; kw...) -> Lattice

Shortcut for `build_lattice(Square, L, M; kw...)`.
"""
square(L::Int, M::Int=L; kw...) = build_lattice(Square, L, M; kw...)

"""
    triangular(L::Int, M::Int = L; kw...) -> Lattice

Shortcut for `build_lattice(Triangular, L, M; kw...)`.
"""
triangular(L::Int, M::Int=L; kw...) = build_lattice(Triangular, L, M; kw...)

"""
    honeycomb(L::Int, M::Int = L; kw...) -> Lattice

Shortcut for `build_lattice(Honeycomb, L, M; kw...)`.
"""
honeycomb(L::Int, M::Int=L; kw...) = build_lattice(Honeycomb, L, M; kw...)

"""
    kagome(L::Int, M::Int = L; kw...) -> Lattice

Shortcut for `build_lattice(Kagome, L, M; kw...)`.
"""
kagome(L::Int, M::Int=L; kw...) = build_lattice(Kagome, L, M; kw...)

"""
    lieb(L::Int, M::Int = L; kw...) -> Lattice

Shortcut for `build_lattice(Lieb, L, M; kw...)`.
"""
lieb(L::Int, M::Int=L; kw...) = build_lattice(Lieb, L, M; kw...)

"""
    union_jack(L::Int, M::Int = L; kw...) -> Lattice

Shortcut for `build_lattice(UnionJack, L, M; kw...)`.
"""
union_jack(L::Int, M::Int=L; kw...) = build_lattice(UnionJack, L, M; kw...)

"""
    dice(L::Int, M::Int = L; kw...) -> Lattice

Shortcut for `build_lattice(Dice, L, M; kw...)`.
"""
dice(L::Int, M::Int=L; kw...) = build_lattice(Dice, L, M; kw...)

"""
    shastry_sutherland(L::Int, M::Int = L; kw...) -> Lattice

Shortcut for `build_lattice(ShastrySutherland, L, M; kw...)`.
"""
shastry_sutherland(L::Int, M::Int=L; kw...) =
    build_lattice(ShastrySutherland, L, M; kw...)
