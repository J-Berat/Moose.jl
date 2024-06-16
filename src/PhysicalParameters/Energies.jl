"""
    MagE(B::AbstractArray) -> AbstractArray

Calculate the magnetic energy density for a given magnetic field array.

# Arguments
- `B::AbstractArray`: An array representing the magnetic field.

# Returns
- `AbstractArray`: An array representing the magnetic energy density.

# Description
This function calculates the magnetic energy density using the formula:
E = B.^2 / (8π)
where B is the magnetic field and pi is the mathematical constant pi.

# Example
```julia
# Example usage
B = rand(100, 100, 100)  # Example magnetic field data
energy_density = MagE(B)
println(energy_density)
"""
MagE(B::AbstractArray) = B.^2 / (8π)

"""
    KE(n::AbstractArray, v::AbstractArray) -> AbstractArray

Calculate the kinetic energy density for a given density and velocity arrays.

# Arguments
- `n::AbstractArray`: An array representing the mass density.
- `v::AbstractArray`: An array representing the velocity.

# Returns
- `AbstractArray`: An array representing the kinetic energy density.

# Description
This function calculates the kinetic energy density using the formula:
KE = 0.5 * n .* v.^2
where n is the mass density and v is the velocity.

# Example
```julia
# Example usage
n = rand(100, 100, 100)  # Example density data
v = rand(100, 100, 100)  # Example velocity data
kinetic_energy_density = KE(n, v)
println(kinetic_energy_density)
"""
KE(n::AbstractArray, v::AbstractArray) = 0.5 * n .* v.^2
