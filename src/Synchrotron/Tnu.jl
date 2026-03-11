"""
    Tnu(Bperp::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> AbstractArray

Calculate the brightness temperature T as a function of frequency.

# Arguments
- `Bperp::AbstractArray`: Array representing the perpendicular component of the magnetic field.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the brightness temperature.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular emissivity), and `e_para` (parallel emissivity).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `AbstractArray`: An array representing the brightness temperature T as a function of frequency.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperp = [1.0, 2.0]
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
T_nu = Tnu(Bperp, nuArray, df, PixelLength_cm)
"""

struct TemperatureInterpolator
    B::Vector{Float64}
    eps_interp::Spline2D
end

function TemperatureInterpolator(df::DataFrame)
    B = collect(unique(df.B))
    nu = collect(unique(df.nu))
    eps = reshape(df.e_para .+ df.e_perp, (length(nu), length(B)))
    eps_interp = Spline2D(B, nu, eps)
    return TemperatureInterpolator(B, eps_interp)
end

function emissivity_total_at_frequency!(buffer, B::Vector{Float64}, eps_interp::Spline2D, Bperp::AbstractArray, nui)
    eps_i = @. eps_interp(B, nui)
    eps_i_interp = linear_interpolation(B, eps_i, extrapolation_bc = Line())
    buffer .= @. eps_i_interp(Bperp)
    return buffer
end

function Tnu(Bperp, nuArray, df, PixelLength_cm; precomputed_interp = nothing)
    interpolator = precomputed_interp === nothing ? TemperatureInterpolator(df) : precomputed_interp
    Nfreq = length(nuArray)
    T_nu = zeros(Nfreq)
    eps_buffer = similar(Bperp, Float64)

    # Tnu computation
    for i in 1:Nfreq
        nui = nuArray[i]

        emissivity_total_at_frequency!(eps_buffer, interpolator.B, interpolator.eps_interp, Bperp, nui)

        Inui = sum(eps_buffer) .* PixelLength_cm
        T_nu[i] = BrightnessTemperature(nui, Inui)
    end
    return T_nu

end

"""
    Tnu3D(Bperpcube::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> AbstractArray

Calculate the brightness temperature T for a 3D cube as a function of frequency.

# Arguments
- `Bperpcube::AbstractArray`: 3D array representing the perpendicular component of the magnetic field for each pixel in the cube.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the brightness temperature.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular electric field), and `e_para` (parallel electric field).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `AbstractArray`: A 3D array representing the brightness temperature T as a function of frequency for each pixel in the cube.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperpcube = rand(10, 10, 2)  # 10x10 pixels, 2 depth slices
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
T_nu = Tnu3D(Bperpcube, nuArray, df, PixelLength_cm)
"""

function Tnu3D(Bperpcube, nuArray, df, PixelLength_cm)
    nx, ny = size(Bperpcube, 1), size(Bperpcube, 2)
    Nfreq = length(nuArray)
    T_nu = zeros(nx, ny, Nfreq)
    interpolator = TemperatureInterpolator(df)

    Threads.@threads for idx in CartesianIndices((1:nx, 1:ny))
        i, j = idx[1], idx[2]
        @views T_nu[i, j, :] = Tnu(Bperpcube[i, j, :], nuArray, df, PixelLength_cm; precomputed_interp = interpolator)
    end

    return T_nu
end
