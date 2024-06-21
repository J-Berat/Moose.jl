"""
    ReadSimulation(simu::String, LOS::String) -> Tuple{AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, Union{AbstractArray, Nothing}, Union{AbstractArray, Nothing}}

Read and process simulation data from FITS files for a specified line of sight (LOS).

# Arguments
- `simu`: The directory containing the simulation FITS files.
- `LOS::String`: The line of sight direction, either "x", "y", or "z".

# Returns
- `Tuple`: A tuple containing the following elements:
  - `B1::AbstractArray`: Magnetic field component perpendicular to the LOS.
  - `B2::AbstractArray`: Another magnetic field component perpendicular to the LOS.
  - `BLOS::AbstractArray`: Magnetic field component along the LOS.
  - `V1::AbstractArray`: Velocity component perpendicular to the LOS.
  - `V2::AbstractArray`: Another velocity component perpendicular to the LOS.
  - `VLOS::AbstractArray`: Velocity component along the LOS.
  - `T::AbstractArray`: Temperature array.
  - `n::AbstractArray`: Density array.
  - `nH2::Union{AbstractArray, Nothing}`: Molecular hydrogen density array (or `nothing` if the file is not present).
  - `nHp::Union{AbstractArray, Nothing}`: Ionized hydrogen density array (or `nothing` if the file is not present).

# Description
This function reads various physical quantities from FITS files located in the specified simulation directory. The files include magnetic field components, velocity components, temperature, density, and optionally, molecular and ionized hydrogen densities. The function processes these quantities based on the specified line of sight (LOS) direction and returns them in a tuple.

# Example
```julia
# Example usage
simu_dir = "path/to/simulation"
LOS = "z"
(B1, B2, BLOS, V1, V2, VLOS, T, n, nH2, nHp) = ReadSimulation(simu_dir, LOS)
"""

function ReadSimulation(simu, LOS, conversionn, conversionT, conversionV, conversionB)
    fileBx = "$simu/Bx.fits"
    fileBy = "$simu/By.fits"
    fileBz = "$simu/Bz.fits"
    filen = "$simu/density.fits"
    fileT = "$simu/temperature.fits"
    fileVx = "$simu/Vx.fits"
    fileVy = "$simu/Vy.fits"
    fileVz = "$simu/Vz.fits"
    filenH2 = "$simu/densityH2.fits"
    filenHp = "$simu/densityHp.fits"
    
    # Reading files that are always present
    T = read(FITS(fileT)[1]) .* conversionT # K
    n = read(FITS(filen)[1]) .* conversionn # cm^-3
    
    # Reading files based on LOS and their presence
    if LOS == "z"
        B1 = read(FITS(fileBx)[1]) .* conversionB # microG
        B2 = read(FITS(fileBy)[1]) .* conversionB # microG
        BLOS = read(FITS(fileBz)[1]) .* conversionB # microG
        V1 = read(FITS(fileVx)[1]) .* conversionV # km/s
        V2 = read(FITS(fileVy)[1]) .* conversionV # km/s
        VLOS = read(FITS(fileVz)[1]) .* conversionV # km/s
        if isfile(filenH2) && isfile(filenHp)
            nH2 = read(FITS(filenH2)[1]) .* conversionn # cm^-3
            nHp = read(FITS(filenHp)[1]) .* conversionn # cm^-3
        else
            nH2 = nothing
            nHp = nothing
        end
    elseif LOS == "y"
        B1 = read(FITS(fileBx)[1]) .* conversionB # microG
        B2 = read(FITS(fileBz)[1]) .* conversionB # microG
        BLOS = read(FITS(fileBy)[1]) .* conversionB # microG
        V1 = read(FITS(fileVx)[1]) .* conversionV # km/s
        V2 = read(FITS(fileVz)[1]) .* conversionV # km/s
        VLOS = read(FITS(fileVy)[1]) .* conversionV # km/s
        if isfile(filenH2) && isfile(filenHp)
            nH2 = read(FITS(filenH2)[1]) .* conversionn # cm^-3
            nHp = read(FITS(filenHp)[1]) .* conversionn # cm^-3
        else
            nH2 = nothing
            nHp = nothing
        end
        # Permuting dimensions for y LOS
        B1 = permutedims(B1, [1, 3, 2])
        B2 = permutedims(B2, [1, 3, 2])
        BLOS = permutedims(BLOS, [1, 3, 2])
        V1 = permutedims(V1, [1, 3, 2])
        V2 = permutedims(V2, [1, 3, 2])
        VLOS = permutedims(VLOS, [1, 3, 2])
        T = permutedims(T, [1, 3, 2])
        n = permutedims(n, [1, 3, 2])
    else
        B1 = read(FITS(fileBy)[1]) .* conversionB # microG
        B2 = read(FITS(fileBz)[1]) .* conversionB # microG
        BLOS = read(FITS(fileBx)[1]) .* conversionB # microG
        V1 = read(FITS(fileVy)[1]) .* conversionV # km/s
        V2 = read(FITS(fileVz)[1]) .* conversionV # km/s
        VLOS = read(FITS(fileVx)[1]) .* conversionV # km/s
        if isfile(filenH2) && isfile(filenHp)
            nH2 = read(FITS(filenH2)[1]) .* conversionn # cm^-3
            nHp = read(FITS(filenHp)[1]) .* conversionn # cm^-3
        else
            nH2 = nothing
            nHp = nothing
        end
        # Permuting dimensions for x LOS
        B1 = permutedims(B1, [3, 2, 1])
        B2 = permutedims(B2, [3, 2, 1])
        BLOS = permutedims(BLOS, [3, 2, 1])
        V1 = permutedims(V1, [3, 2, 1])
        V2 = permutedims(V2, [3, 2, 1])
        VLOS = permutedims(VLOS, [3, 2, 1])
        T = permutedims(T, [3, 2, 1])
        n = permutedims(n, [3, 2, 1])
    end
    return (B1, B2, BLOS, V1, V2, VLOS, T, n, nH2, nHp)
end