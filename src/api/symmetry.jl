"""
    Lattice2D point-group symmetry API
    ==================================

Realized **point-group** symmetries of a lattice sample and the site orbits
they generate (issue #33). A symmetry is an orthogonal transformation
(rotation, mirror reflection, or inversion) about a fixed centre that maps
every site of `lat` onto a site of `lat`; on a periodic sample the mapping is
taken modulo the supercell.

The group is discovered, not assumed: the 24 crystallographic candidate
operations (rotations in 30° steps, mirror lines in 15° steps — enough to
express every 2-, 3-, 4- and 6-fold point group) are each tested against the
sample, and only those that permute the sites are kept. This returns the
*realized* group of the given sample, which depends on its shape and boundary
conditions:

* a fully periodic square `L×L` realizes `C₄ᵥ` (order 8), triangular `C₆ᵥ`
  (order 12), honeycomb the site group `C₃ᵥ` (order 6);
* an open `L×L` square realizes the sample's `D₄` (order 8).

Accordingly the default rotation/reflection centre is the lattice origin for a
fully periodic sample (a lattice point, so the torus point group is realized)
and the site centroid for a sample with any open axis (the sample-shape point
group). Override it with the `center` keyword.
"""

# ---- symmetry operation ---------------------------------------------

"""
    SymmetryOperation

A realized point-group operation of a lattice: the `2×2` orthogonal `matrix`
acting on real-space positions, together with the induced site `permutation`
(`permutation[i]` is the site that site `i` is mapped to).

Query it with [`is_rotation`](@ref), [`is_reflection`](@ref),
[`is_inversion`](@ref) and [`rotation_angle`](@ref).
"""
struct SymmetryOperation
    matrix::SMatrix{2,2,Float64,4}
    permutation::Vector{Int}
end

"""
    is_rotation(op::SymmetryOperation) -> Bool

`true` if `op` is a proper rotation (`det = +1`), including the identity and
the 180° inversion.
"""
is_rotation(op::SymmetryOperation) = det(op.matrix) > 0

"""
    is_reflection(op::SymmetryOperation) -> Bool

`true` if `op` is an improper operation — a mirror reflection (`det = -1`).
"""
is_reflection(op::SymmetryOperation) = det(op.matrix) < 0

"""
    is_inversion(op::SymmetryOperation) -> Bool

`true` if `op` is the point inversion `r ↦ -r` (in 2D, the 180° rotation).
"""
is_inversion(op::SymmetryOperation) = isapprox(op.matrix, -I; atol=1e-8)

"""
    rotation_angle(op::SymmetryOperation) -> Float64

Rotation angle in radians (in `(-π, π]`) for a proper rotation. For a
reflection it returns the angle of the equivalent orthogonal matrix; use
[`is_rotation`](@ref) to guard.
"""
rotation_angle(op::SymmetryOperation) = atan(op.matrix[2, 1], op.matrix[1, 1])

# ---- candidate crystallographic operations --------------------------

# 12 rotations (30° steps) + 12 mirror lines (15° steps): a superset of every
# 2D crystallographic point-group operation. Non-symmetries are filtered out.
function _candidate_matrices()
    mats = SMatrix{2,2,Float64,4}[]
    for k in 0:11
        θ = k * π / 6
        push!(mats, SMatrix{2,2}(cos(θ), sin(θ), -sin(θ), cos(θ)))
    end
    for m in 0:11
        φ = m * π / 12
        c2, s2 = cos(2φ), sin(2φ)
        push!(mats, SMatrix{2,2}(c2, s2, s2, -c2))
    end
    return mats
end

# ---- realized symmetry operations -----------------------------------

"""
    symmetry_operations(lat::Lattice; center=:auto, tol=1e-7)
        -> Vector{SymmetryOperation}

Realized point-group operations of `lat` — every rotation / reflection /
inversion about `center` that permutes the sites (modulo the supercell on
periodic axes). Always contains the identity.

`center` is `:auto` (origin for a fully periodic sample, centroid otherwise),
`:origin`, `:centroid`, or an explicit 2-vector. `tol` is the position-matching
tolerance.

See also [`symmetry_orbits`](@ref), [`symmetry_group_order`](@ref).
"""
function symmetry_operations(lat::Lattice; center=:auto, tol::Real=1e-7)
    N = num_sites(lat)
    P = [SVector{2,Float64}(position(lat, i)) for i in 1:N]

    ax = boundary(lat).axes
    per1 = ax[1] isa PeriodicAxis
    per2 = ax[2] isa PeriodicAxis
    c = _resolve_center(center, P, per1, per2)

    B = basis_vectors(lat)
    a1 = SVector(B[1, 1], B[2, 1])
    a2 = SVector(B[1, 2], B[2, 2])
    M = hcat(lat.Lx * a1, lat.Ly * a2)
    Minv = inv(M)

    key(p) = _reduce_key(p, Minv, per1, per2, Float64(tol))
    lookup = Dict{Tuple{Int,Int},Int}()
    for j in 1:N
        lookup[key(P[j])] = j
    end

    ops = SymmetryOperation[]
    for R in _candidate_matrices()
        perm = Vector{Int}(undef, N)
        seen = falses(N)
        ok = true
        for i in 1:N
            p = c + R * (P[i] - c)
            j = get(lookup, key(p), 0)
            if j == 0 || seen[j]
                ok = false
                break
            end
            perm[i] = j
            seen[j] = true
        end
        ok && push!(ops, SymmetryOperation(R, perm))
    end
    return ops
end

function _resolve_center(center, P, per1, per2)
    center === :origin && return zero(SVector{2,Float64})
    center === :centroid && return sum(P) / length(P)
    if center === :auto
        return (per1 && per2) ? zero(SVector{2,Float64}) : sum(P) / length(P)
    end
    return SVector{2,Float64}(center[1], center[2])
end

# Canonical integer key of a position, reduced into the supercell on periodic
# axes so that torus-equivalent points collide.
@inline function _reduce_key(p, Minv, per1::Bool, per2::Bool, tol::Float64)
    f = Minv * p
    f1 = per1 ? mod(f[1], 1.0) : f[1]
    f2 = per2 ? mod(f[2], 1.0) : f[2]
    # snap values a hair below 1 back to 0 so the wrap seam does not split
    f1 = f1 > 1 - tol ? f1 - 1 : f1
    f2 = f2 > 1 - tol ? f2 - 1 : f2
    return (round(Int, f1 / tol), round(Int, f2 / tol))
end

"""
    symmetry_group_order(lat::Lattice; kwargs...) -> Int

Number of realized point-group operations of `lat`
(`length(symmetry_operations(lat; kwargs...))`).
"""
symmetry_group_order(lat::Lattice; kwargs...) = length(symmetry_operations(lat; kwargs...))

# ---- site orbits -----------------------------------------------------

"""
    symmetry_orbits(lat::Lattice; kwargs...) -> Vector{Vector{Int}}

Partition the sites of `lat` into **orbits** under its realized point group:
two sites are in the same orbit iff some symmetry maps one to the other. Each
orbit is returned sorted, and the orbits are ordered by their smallest site
index. By the orbit–stabilizer theorem every orbit size divides
[`symmetry_group_order`](@ref).

`kwargs` are forwarded to [`symmetry_operations`](@ref) (`center`, `tol`).

# Example
```julia
lat = square(6, 6)            # periodic C₄ᵥ about the origin
orbits = symmetry_orbits(lat) # sites grouped by C₄ᵥ equivalence
```
"""
function symmetry_orbits(lat::Lattice; kwargs...)
    N = num_sites(lat)
    ops = symmetry_operations(lat; kwargs...)
    seen = falses(N)
    orbits = Vector{Int}[]
    for i in 1:N
        seen[i] && continue
        orbit = sort!(unique(op.permutation[i] for op in ops))
        for j in orbit
            seen[j] = true
        end
        push!(orbits, orbit)
    end
    return orbits
end

"""
    site_orbit(lat::Lattice, i::Int; kwargs...) -> Vector{Int}

The orbit (sorted site indices) containing site `i` under the realized point
group of `lat`. `kwargs` are forwarded to [`symmetry_operations`](@ref).
"""
function site_orbit(lat::Lattice, i::Int; kwargs...)
    ops = symmetry_operations(lat; kwargs...)
    return sort!(unique(op.permutation[i] for op in ops))
end
