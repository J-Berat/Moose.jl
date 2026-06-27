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
    B, nu, eps = emissivity_grid(df, df.e_para .+ df.e_perp)
    eps_interp = Spline2D(B, nu, eps)
    return TemperatureInterpolator(B, eps_interp)
end

function build_emissivity_frequency_cache(interpolator::TemperatureInterpolator, nuArray)
    Nfreq = length(nuArray)
    B = interpolator.B
    cache = Matrix{Float64}(undef, length(B), Nfreq)
    @inbounds for i in 1:Nfreq
        nui = nuArray[i]
        for j in eachindex(B)
            cache[j, i] = interpolator.eps_interp(B[j], nui)
        end
    end
    return cache
end

function emissivity_total_at_frequency!(buffer, B::Vector{Float64}, eps_interp::Spline2D, Bperp::AbstractArray, nui;
    eps_cache_col=nothing, eps_line_buffer=nothing)
    eps_i = eps_cache_col
    if eps_i === nothing
        eps_i = eps_line_buffer
        @inbounds for j in eachindex(B)
            eps_i[j] = eps_interp(B[j], nui)
        end
    end

    @inbounds for idx in eachindex(Bperp, buffer)
        buffer[idx] = linear_interp_extrapolated(B, eps_i, Float64(Bperp[idx]))
    end

    return buffer
end

function _Tnu!(T_nu, Bperp, nuArray, PixelLength_cm, interpolator; emissivity_cache=nothing)
    eps_buffer = similar(Bperp, Float64)
    eps_line_buffer = emissivity_cache === nothing ? similar(interpolator.B, Float64) : nothing

    for i in eachindex(nuArray)
        nui = nuArray[i]
        cache_col = emissivity_cache === nothing ? nothing : view(emissivity_cache, :, i)
        emissivity_total_at_frequency!(eps_buffer, interpolator.B, interpolator.eps_interp, Bperp, nui;
            eps_cache_col=cache_col, eps_line_buffer=eps_line_buffer)

        Inui = sum(eps_buffer) * PixelLength_cm
        T_nu[i] = BrightnessTemperature(nui, Inui)
    end

    return T_nu
end

function Tnu(Bperp, nuArray, df, PixelLength_cm; precomputed_interp = nothing, emissivity_cache=nothing)
    interpolator = precomputed_interp === nothing ? TemperatureInterpolator(df) : precomputed_interp
    Nfreq = length(nuArray)
    T_nu = zeros(Nfreq)
    return _Tnu!(T_nu, Bperp, nuArray, PixelLength_cm, interpolator; emissivity_cache=emissivity_cache)

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
    emissivity_cache = build_emissivity_frequency_cache(interpolator, nuArray)

    Threads.@threads for idx in CartesianIndices((1:nx, 1:ny))
        i, j = idx[1], idx[2]
        @views Bperp_vec = Bperpcube[i, j, :]
        @views tdest = T_nu[i, j, :]
        _Tnu!(tdest, Bperp_vec, nuArray, PixelLength_cm, interpolator; emissivity_cache=emissivity_cache)
    end

    return T_nu
end
