"""
    PolarizationFraction(P::AbstractArray, I::AbstractArray) -> AbstractArray

Calculate the polarization fraction for given polarized intensity and total intensity arrays.

# Arguments
- `P::AbstractArray`: An array representing the polarized intensity.
- `I::AbstractArray`: An array representing the total intensity.

# Returns
- `AbstractArray`: An array representing the polarization fraction.

# Description
This function calculates the polarization fraction using the formula:
PolarizationFraction = P ./ I
where P is the polarized intensity and I is the total intensity.

# Example
```julia
# Example usage
P = rand(100, 100)  # Example polarized intensity data
I = rand(100, 100)  # Example total intensity data
polarization_fraction = PolarizationFraction(P, I)
println(polarization_fraction)
"""
PolarizationFraction(P::AbstractArray, I::AbstractArray) = P ./ I