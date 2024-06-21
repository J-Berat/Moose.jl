h = 6.626e-34 # J.s\n,
kb = 1.380e-23 # J.K-1
c = 2.99e8 # m.s-1
"""
    Bnu(nu::AbstractArray, T::AbstractArray) -> AbstractArray

Calculate the spectral radiance (Planck function) for given frequency and temperature arrays.

# Arguments
- `nu::AbstractArray`: An array representing the frequency.
- `T::AbstractArray`: An array representing the temperature.

# Returns
- `AbstractArray`: An array representing the spectral radiance in Jy/sr.

# Description
This function calculates the spectral radiance (Bnu) using the Planck function:
nu = 2 * h * nu^3 / (c^2 * (exp((h * nu) / (kb * T)) - 1)) * Wm2Hz_to_Jy 
where:
- h is the Planck constant,
- nu is the frequency,
- c is the speed of light,
- k_B is the Boltzmann constant,
- T is the temperature,
- Wm2Hz_to_Jy is the conversion factor from W/m²/Hz to Jy.

# Example
```julia
# Example usage
nu = 1e9 .* rand(100)  # Example frequency data in Hz
T = 2.725 .* ones(100)  # Example temperature data in K (CMB temperature)
spectral_radiance = Bnu(nu, T)
println(spectral_radiance)
"""
Bnu(nu::AbstractArray, T::AbstractArray) = @. 2 * h * nu^3 / (c^2 * (exp((h * nu) / (kb * T)) - 1)) * Wm2Hz_to_Jy # Jy.sr-1