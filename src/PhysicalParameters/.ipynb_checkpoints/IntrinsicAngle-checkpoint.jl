"""
    IntrinsicAngle(B2::AbstractArray, B1::AbstractArray) -> AbstractArray

Calculate the intrinsic angle of polarization for given magnetic field components.

# Arguments
- `B2::AbstractArray`: An array representing the second component of the magnetic field.
- `B1::AbstractArray`: An array representing the first component of the magnetic field.

# Returns
- `AbstractArray`: An array representing the intrinsic angle of polarization.

# Description
This function calculates the intrinsic angle of polarization using the formula:
angle = atan.(B2, B1) .+ π / 2
where B2 and B1 are the magnetic field components.

# Example
```julia
# Example usage
B1 = rand(100, 100)  # Example magnetic field component B1
B2 = rand(100, 100)  # Example magnetic field component B2
angle = IntrinsicAngle(B2, B1)
println(angle)
"""
IntrinsicAngle(B2::AbstractArray, B1::AbstractArray) = atan.(B2, B1)