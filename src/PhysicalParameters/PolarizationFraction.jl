"""
    PolarizationFraction(P::AbstractArray, I::AbstractArray) -> AbstractArray

Calculate the polarization fraction for given polarized intensity and total
intensity arrays.

# Arguments
- `P::AbstractArray`: An array representing the polarized intensity.
- `I::AbstractArray`: An array representing the total intensity.

# Returns
- `AbstractArray`: An array representing the polarization fraction.

# Description
This function calculates the polarization fraction using the formula:
PolarizationFraction = P ./ I where P is the polarized intensity and I is the
total intensity. Pixels with non-finite values or non-positive total intensity
are returned as `NaN`, because the fraction is undefined there.

# Example
```julia
P = rand(100, 100)  # Example polarized intensity data
I = rand(100, 100)  # Example total intensity data
polarization_fraction = PolarizationFraction(P, I)
println(polarization_fraction)
```
"""
function PolarizationFraction(P::AbstractArray, I::AbstractArray)
    axes(P) == axes(I) || error("P and I must have the same axes, got $(axes(P)) and $(axes(I)).")
    T = Base.promote_op(/, eltype(P), eltype(I))
    out = similar(P, T)
    nan = convert(T, NaN)

    @inbounds for idx in eachindex(P, I, out)
        p = P[idx]
        i = I[idx]
        out[idx] = isfinite(p) && isfinite(i) && i > zero(i) ? p / i : nan
    end

    return out
end
