"""
    PolarizationAngle(U::AbstractArray, Q::AbstractArray) -> AbstractArray

Calculate the polarization angle for given Stokes parameters U and Q.

# Arguments
- `U::AbstractArray`: An array representing the Stokes parameter U.
- `Q::AbstractArray`: An array representing the Stokes parameter Q.

# Returns
- `AbstractArray`: An array representing the polarization angle in degrees.

# Description
This function calculates the polarization angle using the formula:
PolarizationAngle = 1/2 .* atan.(U,Q)
The resulting angle is then converted from radians to degrees.

# Example
```julia
# Example usage
U = rand(100, 100)  # Example Stokes parameter U data
Q = rand(100, 100)  # Example Stokes parameter Q data
polarization_angle = PolarizationAngle(U, Q)
println(polarization_angle)

"""
PolarizationAngle(U::AbstractArray, Q::AbstractArray) = rad2deg.(1/2 .* atan.(U,Q))