"""
    pressure(n::AbstractArray, T::AbstractArray) -> AbstractArray

Calculate the pressure for given density and temperature arrays.

# Arguments
- `n::AbstractArray`: An array representing the number density.
- `T::AbstractArray`: An array representing the temperature.

# Returns
- `AbstractArray`: An array representing the pressure.

# Description
This function calculates the pressure using the ideal gas law:
P = n * T 
where n is the number density and T is the temperature.

# Example
```julia
# Example usage
n = rand(100, 100)  # Example number density data
T = rand(100, 100)  # Example temperature data
pressure_values = pressure(n, T)
println(pressure_values)
"""
pressure(n::AbstractArray, T::AbstractArray) = n .* T