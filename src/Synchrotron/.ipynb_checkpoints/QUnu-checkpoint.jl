"""
    QUnu(Bperp::AbstractArray, psi_src::AbstractArray, RM::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> Tuple{AbstractArray, AbstractArray}

Calculate the Stokes parameters Q and U as functions of frequency.

# Arguments
- `Bperp::AbstractArray`: Array representing the perpendicular component of the magnetic field.
- `psi_src::AbstractArray`: Array representing the source polarization angles.
- `RM::AbstractArray`: Array representing the Rotation Measure.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the Stokes parameters.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular emissivity), and `e_para` (parallel emissivity).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `Tuple{AbstractArray, AbstractArray}`: Two arrays, Qnu and Unu, representing the Stokes parameters Q and U as functions of frequency.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperp = [1.0, 2.0]
psi_src = [0.1, 0.2]
RM = [0.001, 0.002]
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
Qnu, Unu = QUnu(Bperp, psi_src, RM, nuArray, df, PixelLength_cm)

"""

function QUnu(Bperp, psi_src, RM, nuArray, df, PixelLength_cm)

    B = unique(df.B)
    nu = unique(df.nu)
    eps = reshape(df.e_perp .- df.e_para, (size(nu,1), size(B,1)))

    # 2D interpolation function 
    eps_interp = Spline2D(B, nu, eps)
    
    Nfreq = length(nuArray)
    Qnu = zeros(Nfreq)
    Unu = zeros(Nfreq)

    # QUnu computation       
    Threads.@threads for i = 1:Nfreq
        nui = nuArray[i]

        eps_i = @. eps_interp(B, nui)  #interpolation vector at frequency nui
        eps_i_interp = linear_interpolation(B, eps_i, extrapolation_bc=Line()) # 1D interpolation function
        eps_i = @. eps_i_interp(Bperp)  # interpolate only on B over the full cube

        FaradayAngle = RM .* (C_m / (nui * 1e6))^2
        argument = 2 .* (psi_src .+ FaradayAngle)
        
        integrande_U = eps_i .* sin.(argument)
        integrande_Q = eps_i .* cos.(argument)
        Unui =  sum(integrande_U) .* PixelLength_cm
        Qnui =  sum(integrande_Q) .* PixelLength_cm
    
        Unu[i] = @. BrightnessTemperature(nui,Unui)
        Qnu[i] = @. BrightnessTemperature(nui,Qnui)
    end
    return Qnu, Unu
    
end
"""
    QUnu3D(Bperpcube::AbstractArray, psi_src::AbstractArray, RM::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> Tuple{AbstractArray, AbstractArray}

Calculate the Stokes parameters Q and U for a 3D cube as functions of frequency.

# Arguments
- `Bperpcube::AbstractArray`: 3D array representing the perpendicular component of the magnetic field for each pixel in the cube.
- `psi_src::AbstractArray`: 3D array representing the source polarization angles for each pixel in the cube.
- `RM::AbstractArray`: 3D array representing the Rotation Measure for each pixel in the cube.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the Stokes parameters.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular emissivity), and `e_para` (parallel emissivity).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `Tuple{AbstractArray, AbstractArray}`: Two 3D arrays, Qnu and Unu, representing the Stokes parameters Q and U as functions of frequency for each pixel in the cube.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperpcube = rand(10, 10, 2)  # 10x10 pixels, 2 depth slices
psi_src = rand(10, 10, 2)
RM = rand(10, 10, 2)
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
Qnu, Unu = QUnu3D(Bperpcube, psi_src, RM, nuArray, df, PixelLength_cm)
"""
function QUnu3D(Bperpcube, psi_src, RM, nuArray, df, PixelLength_cm)
    
    Nfreq = length(nuArray)
    nx,ny = size(Bperpcube,1), size(Bperpcube,2)
    Qnu = zeros(nx, ny, Nfreq)
    Unu = zeros(nx, ny, Nfreq)
    
    for i = 1:nx
        for j = 1:ny
            Bperp_vec = Bperpcube[i,j,:]
            #computing RM
            RM_vec = RM[i,j,:]
            #computing Stokes Qnu, Unu
            psi_src_vec = psi_src[i,j,:] #instrinsic angle of polarization
            Qnu[i,j,:], Unu[i,j,:] = QUnu(Bperp_vec, psi_src_vec, RM_vec, nuArray, df, PixelLength_cm)
        end
    end

    return Qnu, Unu

end

"""
    QUnuNoFaraday(Bperp::AbstractArray, psi_src::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> Tuple{AbstractArray, AbstractArray}

Calculate the Stokes parameters Q and U as functions of frequency without considering Faraday rotation.

# Arguments
- `Bperp::AbstractArray`: Array representing the perpendicular component of the magnetic field.
- `psi_src::AbstractArray`: Array representing the source polarization angles.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the Stokes parameters.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular electric field), and `e_para` (parallel electric field).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `Tuple{AbstractArray, AbstractArray}`: Two arrays, Qnu and Unu, representing the Stokes parameters Q and U as functions of frequency.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperp = [1.0, 2.0]
psi_src = [0.1, 0.2]
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
Qnu, Unu = QUnuNoFaraday(Bperp, psi_src, nuArray, df, PixelLength_cm)
"""
function QUnuNoFaraday(Bperp, psi_src, nuArray, df, PixelLength_cm)
    
    B = unique(df.B)
    nu = unique(df.nu)
    eps = reshape(df.e_perp .- df.e_para, (size(nu,1), size(B,1)))

    # 2D interpolation function 
    eps_interp = Spline2D(B, nu, eps)
    
    Nfreq = length(nuArray)
    Qnu = zeros(Nfreq)
    Unu = zeros(Nfreq)

    # QUnu computation       
    Threads.@threads for i = 1:Nfreq
        nui = nuArray[i]

        eps_i = @. eps_interp(B, nui)  #interpolation vector at frequency nui
        eps_i_interp = linear_interpolation(B, eps_i, extrapolation_bc=Line()) # 1D interpolation function
        eps_i = @. eps_i_interp(Bperp)  # interpolate only on B over the full cube
        
        argument = 2 .* psi_src 
        
        integrande_U = eps_i .* sin.(argument)
        integrande_Q = eps_i .* cos.(argument)
        Unui =  sum(integrande_U) .* PixelLength_cm
        Qnui =  sum(integrande_Q) .* PixelLength_cm
    
        Unu[i] = @. BrightnessTemperature(nui,Unui)
        Qnu[i] = @. BrightnessTemperature(nui,Qnui)
    end
    return Qnu, Unu
    
end

"""
    QUnuNoFaraday3D(Bperpcube::AbstractArray, psi_src::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> Tuple{AbstractArray, AbstractArray}

Calculate the Stokes parameters Q and U for a 3D cube as functions of frequency without considering Faraday rotation.

# Arguments
- `Bperpcube::AbstractArray`: 3D array representing the perpendicular component of the magnetic field for each pixel in the cube.
- `psi_src::AbstractArray`: 3D array representing the source polarization angles for each pixel in the cube.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the Stokes parameters.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular emissivity), and `e_para` (parallel emissivity).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `Tuple{AbstractArray, AbstractArray}`: Two 3D arrays, Qnu and Unu, representing the Stokes parameters Q and U as functions of frequency for each pixel in the cube.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperpcube = rand(10, 10, 2)  # 10x10 pixels, 2 depth slices
psi_src = rand(10, 10, 2)
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
Qnu, Unu = QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, PixelLength_cm)
"""
function QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, PixelLength_cm)
    
    Nfreq = length(nuArray)
    Qnu = zeros(size(Bperpcube,1), size(Bperpcube,2), Nfreq)
    Unu = zeros(size(Bperpcube,1), size(Bperpcube,2), Nfreq)
    
    Threads.@threads for i = 1:size(Bperpcube,1)
        Threads.@threads for j = 1:size(Bperpcube,2)
            Bperp_vec = Bperpcube[i,j,:]
            #computing Stokes Qnu, Unu
            psi_src_vec = psi_src[i,j,:] #instrinsic angle of polarization
            Qnu[i,j,:], Unu[i,j,:] = QUnuNoFaraday(Bperp_vec, psi_src_vec, nuArray, df, PixelLength_cm)
        end
    end

    return Qnu, Unu

end