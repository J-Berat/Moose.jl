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
# The π/2 offset is converted to the (promoted) element type of the inputs so
# that reduced-precision cubes (e.g. Float32 in `precision = "float32"` runs)
# are not silently promoted back to Float64. For Float64 inputs the result is
# bit-identical to the previous `.+ π / 2` form.
IntrinsicAngle(B2::AbstractArray, B1::AbstractArray) =
    atan.(B2, B1) .+ float(promote_type(eltype(B2), eltype(B1)))(π / 2)