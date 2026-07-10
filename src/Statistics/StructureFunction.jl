"""
Isotropic second-order structure functions of 2D maps.

The structure function of a field `X`,

    SF(r) = ⟨ [X(x) − X(x + r)]² ⟩,

measured on rotation measure or polarization angle maps is a standard
diagnostic of magnetized turbulence: its slope and outer scale constrain the
turbulent cascade of the Faraday-rotating medium (e.g. Minter & Spangler
1996; Haverkorn et al. 2004). It complements the power-spectrum tools of
`PowerSpectrum.jl` and is more robust on masked/irregular fields since it
needs no gridding or apodization.

Pairs are drawn by Monte-Carlo sampling (`npairs` random pixel pairs), which
keeps the cost independent of map size and handles `NaN`-masked maps
naturally. Polarization angles are compared modulo their `π` ambiguity.
"""

"""
    StructureFunctionResult

Result of [`structure_function`](@ref).

# Fields
- `separation::Vector{Float64}`: bin centres (geometric mean of the bin
  edges), in units of `pixel_size`.
- `sf::Vector{Float64}`: structure function ⟨ΔX²⟩ per bin (`NaN` for empty
  bins).
- `counts::Vector{Int}`: number of sampled pairs per bin.
- `edges::Vector{Float64}`: the `nbins + 1` bin edges.
- `npairs::Int`: total number of pair draws requested.
- `angle::Bool`: whether differences were wrapped modulo `angle_period`.
"""
struct StructureFunctionResult
    separation::Vector{Float64}
    sf::Vector{Float64}
    counts::Vector{Int}
    edges::Vector{Float64}
    npairs::Int
    angle::Bool
end

"""
    structure_function(X; pixel_size=1.0, nbins=20, min_sep=pixel_size,
                       max_sep=nothing, npairs=1_000_000, angle=false,
                       angle_period=pi, rng=Random.default_rng())
        -> StructureFunctionResult

Estimate the isotropic second-order structure function of the 2D map `X` by
Monte-Carlo pair sampling.

# Keywords
- `pixel_size::Real`: linear pixel size; separations are reported in this
  unit.
- `nbins::Int`: number of logarithmic separation bins.
- `min_sep`, `max_sep`: separation range (defaults: one pixel to half the map
  diagonal).
- `npairs::Int`: number of random pixel pairs to draw. Pairs falling outside
  the separation range or hitting a non-finite pixel are discarded (they
  still count towards `npairs`).
- `angle::Bool`: set to `true` for orientation maps such as the polarization
  angle. Differences are then wrapped into `±angle_period/2` before squaring,
  so that e.g. `ψ = +89°` and `ψ = −89°` are treated as 2° apart
  (`angle_period = π` matches the `n·π` ambiguity of polarization angles).
- `rng`: random number generator (pass a seeded generator for reproducible
  estimates).

`NaN`/non-finite pixels (masked, HEALPix UNSEEN) are ignored. Use the
`counts` field to judge the reliability of each bin.

# Example
```julia
rm_map = ...                         # rotation measure map, rad/m²
sf = structure_function(rm_map; pixel_size = 2.0, npairs = 2_000_000)
loglog(sf.separation, sf.sf)         # slope → turbulence spectral index

psi = 0.5 .* atan.(U, Q)             # polarization angle map
sfa = structure_function(psi; angle = true)
```
"""
function structure_function(X::AbstractMatrix;
                            pixel_size::Real = 1.0,
                            nbins::Int = 20,
                            min_sep::Real = float(pixel_size),
                            max_sep::Union{Nothing,Real} = nothing,
                            npairs::Int = 1_000_000,
                            angle::Bool = false,
                            angle_period::Real = pi,
                            rng = Random.default_rng())

    pixel_size > 0 || throw(ArgumentError("pixel_size must be positive (got $pixel_size)."))
    nbins >= 1 || throw(ArgumentError("nbins must be at least 1 (got $nbins)."))
    npairs >= 1 || throw(ArgumentError("npairs must be at least 1 (got $npairs)."))
    angle_period > 0 || throw(ArgumentError("angle_period must be positive (got $angle_period)."))

    nx, ny = size(X)
    h = float(pixel_size)
    rmax = max_sep === nothing ? h * hypot(nx, ny) / 2 : float(max_sep)
    rmin = float(min_sep)
    0 < rmin < rmax ||
        throw(ArgumentError("Require 0 < min_sep < max_sep (got $rmin and $rmax)."))

    # Valid (finite) pixels only.
    valid = Tuple{Int,Int}[]
    sizehint!(valid, nx * ny)
    @inbounds for j in 1:ny, i in 1:nx
        isfinite(X[i, j]) && push!(valid, (i, j))
    end
    length(valid) >= 2 ||
        throw(ArgumentError("structure_function needs at least 2 finite pixels."))

    edges = collect(exp.(range(log(rmin), log(rmax); length = nbins + 1)))
    log_rmin = log(rmin)
    inv_dlog = nbins / (log(rmax) - log_rmin)

    sums = zeros(Float64, nbins)
    counts = zeros(Int, nbins)
    half = angle_period / 2

    @inbounds for _ in 1:npairs
        a = valid[rand(rng, 1:length(valid))]
        b = valid[rand(rng, 1:length(valid))]
        a == b && continue
        r = h * hypot(a[1] - b[1], a[2] - b[2])
        (rmin <= r < rmax) || continue
        bin = 1 + floor(Int, (log(r) - log_rmin) * inv_dlog)
        bin = clamp(bin, 1, nbins)
        d = X[a[1], a[2]] - X[b[1], b[2]]
        if angle
            d = mod(d + half, angle_period) - half
        end
        sums[bin] += d * d
        counts[bin] += 1
    end

    sf = [counts[k] > 0 ? sums[k] / counts[k] : NaN for k in 1:nbins]
    centers = [sqrt(edges[k] * edges[k+1]) for k in 1:nbins]

    return StructureFunctionResult(centers, sf, counts, edges, npairs, angle)
end
