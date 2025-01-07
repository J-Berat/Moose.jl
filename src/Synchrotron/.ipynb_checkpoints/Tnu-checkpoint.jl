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

function Tnu(Bperp, nuArray, df, PixelLength_cm)
    
    B = unique(df.B)
    nu = unique(df.nu)
    eps = reshape(df.e_para .+ df.e_perp , (size(nu,1), size(B,1)))
   
    # 2D interpolation function 
    eps_interp = Spline2D(B, nu, eps)
    
    # create T_nu
    Nfreq = length(nuArray)
    T_nu = zeros(Nfreq)

    # Tnu computation       
    Threads.@threads for i in 1:Nfreq
        nui = nuArray[i]
        eps_i = @. eps_interp(B, nui)  # interpolation vector at nui frequency 
        eps_i_interp = linear_interpolation(B, eps_i, extrapolation_bc=Line()) # 1D interpolation function
        eps_i = @. eps_i_interp(Bperp)  # interpolate only on B over the full cube
        Inui = sum(eps_i) .* PixelLength_cm
        T_nu[i] = BrightnessTemperature(nui,Inui)
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
    
    Nfreq = length(nuArray)
    T_nu = zeros(size(Bperpcube,1), size(Bperpcube,2), Nfreq)
    
    Threads.@threads for i = 1:size(T_nu,1)
        Threads.@threads for j = 1:size(T_nu,2)
            T_nu[i,j,:] = Tnu(Bperpcube[i,j,:], nuArray, df, PixelLength_cm)
        end
    end
    return T_nu
end