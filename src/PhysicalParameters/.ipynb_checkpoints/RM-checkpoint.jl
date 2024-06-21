"""
    deltaRM(BLOS::AbstractArray, ne::AbstractArray, PixelLength_pc::Float64) -> AbstractArray

Calculate the differential rotation measure (delta RM) for a given magnetic field, electron density, and pixel length.

# Arguments
- `BLOS::AbstractArray`: Array of the magnetic field component along the line of sight (LOS).
- `ne::AbstractArray`: Array of electron densities.
- `PixelLength_pc::Float64`: The length of a pixel in parsecs.

# Returns
- `AbstractArray`: An array representing the differential rotation measure.

# Description
This function calculates the differential rotation measure (delta RM) using the provided magnetic field component along the line of sight (BLOS), electron density (ne), and pixel length in parsecs (PixelLength_pc). The calculation is based on the formula:
DeltaRM = RM_PREFACTOR * ne * BLOS * PixelLength_pc 
where `RM_PREFACTOR` is a predefined constant.

# Example
```julia
# Example usage
BLOS = randn(100, 100, 100)  # Example BLOS data
ne = rand(100, 100, 100)  # Example electron density data
PixelLength_pc = 0.1  # Example pixel length in parsecs
delta_rm = deltaRM(BLOS, ne, PixelLength_pc)
println(delta_rm)
"""
deltaRM(BLOS::AbstractArray, ne::AbstractArray, PixelLength_pc::Float64) = RM_PREFACTOR .* ne .* BLOS .* PixelLength_pc 

"""
    RM(deltaRM::Array{Float64, 1}) -> Array{Float64, 1}

Calculate the rotation measure (RM) by taking the cumulative sum of the differential rotation measure (delta RM).

# Arguments
- `deltaRM::Array{Float64, 1}`: A 1D array representing the differential rotation measure.

# Returns
- `Array{Float64, 1}`: A 1D array representing the rotation measure.

# Description
This function calculates the rotation measure (RM) by computing the cumulative sum of the differential rotation measure (delta RM). The cumulative sum operation integrates the delta RM values to produce the RM values along the specified axis.

# Example
```julia
# Example usage
delta_rm = randn(100)  # Example differential rotation measure data
rm = RM(delta_rm)
println(rm)
"""
RM(deltaRM::Array{Float64, 1}) = cumsum(deltaRM)

"""
    RM(deltaRM::Array{Float64, 3}) -> Array{Float64, 3}

Calculate the rotation measure (RM) by taking the cumulative sum of the differential rotation measure (delta RM) along the third dimension.

# Arguments
- `deltaRM::Array{Float64, 3}`: A 3D array representing the differential rotation measure.

# Returns
- `Array{Float64, 3}`: A 3D array representing the rotation measure.

# Description
This function calculates the rotation measure (RM) by computing the cumulative sum of the differential rotation measure (delta RM) along the third dimension. The cumulative sum operation integrates the delta RM values along the specified dimension to produce the RM values.

# Example
```julia
# Example usage
delta_rm = randn(100, 100, 100)  # Example differential rotation measure data
rm = RM(delta_rm)
println(rm)
"""
RM(deltaRM::Array{Float64, 3}) = cumsum(deltaRM,dims=3)