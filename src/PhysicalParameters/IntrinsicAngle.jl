"""
    IntrinsicAngle(B1::AbstractArray, B2::AbstractArray) -> AbstractArray

Calculate the intrinsic angle of polarization for given magnetic field components.

# Arguments
- `B1::AbstractArray`: Magnetic-field component along the reference sky axis.
- `B2::AbstractArray`: Magnetic-field component along the second sky axis.

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
angle = IntrinsicAngle(B1, B2)
println(angle)
"""
# The π/2 offset is converted to the (promoted) element type of the inputs so
# that reduced-precision cubes (e.g. Float32 in `precision = "float32"` runs)
# are not silently promoted back to Float64. The public argument order is
# (reference-axis component, second-axis component), as used by all pipelines.
IntrinsicAngle(B1::AbstractArray, B2::AbstractArray) =
    atan.(B2, B1) .+ float(promote_type(eltype(B2), eltype(B1)))(π / 2)