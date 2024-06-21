"""
    BrightnessTemperature(nu_MHz::Float64, I::AbstractArray) -> AbstractArray

Calculate the brightness temperature for a given frequency and intensity.

# Arguments
- `nu_MHz::Float64`: The frequency in megahertz (MHz).
- `I::AbstractArray`: The intensity array.

# Returns
- `AbstractArray`: An array representing the brightness temperature in Kelvin (K).

# Description
This function calculates the brightness temperature for a given frequency and intensity using the formula:
BrightnessTemperature = (C^2 / (2 * K_B * ((nu_MHz * 1e6)^2))) * I 
where:
- T_b is the brightness temperature.
- C is the speed of light.
- K_B is the Boltzmann constant.
- nu_MHz is the frequency in megahertz.
- I is the intensity.

# Example
```julia
# Example usage
nu_MHz = 1400.0  # Frequency in MHz
I = rand(100, 100)  # Example intensity data
brightness_temp = BrightnessTemperature(nu_MHz, I)
println(brightness_temp)

"""
BrightnessTemperature(nu_MHz, I) = (C^2 / (2 * K_B * ((nu_MHz * 1e6)^2))) .* I # K