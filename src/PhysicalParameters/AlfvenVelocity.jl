"""
    AlfvenVelocity(B::AbstractArray, density::AbstractArray) -> AbstractArray

Calculate the Alfvén velocity for given magnetic field strength and density arrays.

# Arguments
- `B::AbstractArray`: An array representing the magnetic field strength.
- `density::AbstractArray`: An array representing the mass density.

# Returns
- `AbstractArray`: An array representing the Alfvén velocity.

# Description
This function calculates the Alfvén velocity using the formula:
v_A = sqrt(B^2 / (4π * density))
where B is the magnetic field strength and density is the mass density.

# Example
```julia
# Example usage
B = rand(100, 100)  # Example magnetic field strength data
density = rand(100, 100)  # Example mass density data
alfven_velocity = AlfvenVelocity(B, density)
println(alfven_velocity)
"""
AlfvenVelocity(B,density) = @. sqrt(B^2 / (4π * density))