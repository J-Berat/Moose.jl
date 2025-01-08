"""
    BrightnessTemperature(nu_MHz::Float64, I::AbstractArray) -> AbstractArray

Calculate the brightness temperature for a given frequency and intensity in CGS units.

# Arguments
- `nu_MHz::Float64`: The frequency in megahertz (MHz).
- `I::AbstractArray`: The intensity array in CGS units (erg/s/cm²/Hz/sr).

# Returns
- `AbstractArray`: An array representing the brightness temperature in Kelvin (K).

# Description
The `BrightnessTemperature` function computes the brightness temperature using the formula:

Tb = (c^2 / (2 * k_B * (nu_MHz * 10^6)^2)) * I

where:
- Tb is the brightness temperature (K).
- c is the speed of light (2.99792 × 10^10 cm/s).
- k_B is the Boltzmann constant (1.38065 × 10^-16 erg/K).
- nu_MHz is the frequency in megahertz (MHz).
- I is the specific intensity in CGS units (erg/s/cm²/Hz/sr).

# Example
```julia
# Define the input frequency in MHz and example intensity array
nu_MHz = 1400.0  # Frequency in MHz
I = rand(100, 100)  # Example intensity data in CGS units (erg/s/cm²/Hz/sr)

# Compute the brightness temperature
brightness_temp = BrightnessTemperature(nu_MHz, I)
"""
BrightnessTemperature(nu_MHz, I) = (C^2 / (2 * K_B * ((nu_MHz * 1e6)^2))) .* I # K