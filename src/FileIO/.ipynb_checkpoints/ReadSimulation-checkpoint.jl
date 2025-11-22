"""
    ReadSimulation(simu::String, LOS::String, conversionn::Number, conversionT::Number, conversionV::Number, conversionB::Number) 
        -> Tuple{AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, Union{AbstractArray, Nothing}, Union{AbstractArray, Nothing}}

Reads and processes simulation data from FITS files for a specified line of sight (LOS).

# Arguments
- `simu::String`: The directory containing the simulation FITS files.
- `LOS::String`: The line of sight direction, either "x", "y", or "z".
- `conversionn::Number`: Conversion factor for density.
- `conversionT::Number`: Conversion factor for temperature.
- `conversionV::Number`: Conversion factor for velocity.
- `conversionB::Number`: Conversion factor for magnetic field.

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
    
    T = read_file(fileT, conversionT)
    n = read_file(filen, conversionn)
    nH2 = read_optional_file(filenH2, conversionn, LOS)
    nHp = read_optional_file(filenHp, conversionn, LOS)

    B1, B2, BLOS = LOS == "z" ? (read_file(fileBx, conversionB), read_file(fileBy, conversionB), read_file(fileBz, conversionB)) :
                    LOS == "y" ? (read_file(fileBz, conversionB), read_file(fileBx, conversionB), read_file(fileBy, conversionB)) :
                                 (read_file(fileBy, conversionB), read_file(fileBz, conversionB), read_file(fileBx, conversionB))

    V1, V2, VLOS = LOS == "z" ? (read_file(fileVx, conversionV), read_file(fileVy, conversionV), read_file(fileVz, conversionV)) :
                    LOS == "y" ? (read_file(fileVz, conversionV), read_file(fileVx, conversionV), read_file(fileVy, conversionV)) :
                                 (read_file(fileVy, conversionV), read_file(fileVz, conversionV), read_file(fileVx, conversionV))

    B1, B2, BLOS = permute_dims(B1, LOS), permute_dims(B2, LOS), permute_dims(BLOS, LOS)
    V1, V2, VLOS = permute_dims(V1, LOS), permute_dims(V2, LOS), permute_dims(VLOS, LOS)
    T, n = permute_dims(T, LOS), permute_dims(n, LOS)

    return (B1, B2, BLOS, V1, V2, VLOS, T, n, nH2, nHp)
end

function ReadSimulation(simu, LOS, conversionn, conversionT, conversionB)
    fileBx = "$simu/Bx.fits"
    fileBy = "$simu/By.fits"
    fileBz = "$simu/Bz.fits"
    filen = "$simu/density.fits"
    fileT = "$simu/temperature.fits"
    filenH2 = "$simu/densityH2.fits"
    filenHp = "$simu/densityHp.fits"
    
    T = read_file(fileT, conversionT)
    n = read_file(filen, conversionn)
    nH2 = read_optional_file(filenH2, conversionn, LOS)
    nHp = read_optional_file(filenHp, conversionn, LOS)

    B1, B2, BLOS = LOS == "z" ? (read_file(fileBx, conversionB), read_file(fileBy, conversionB), read_file(fileBz, conversionB)) :
                    LOS == "y" ? (read_file(fileBx, conversionB), read_file(fileBz, conversionB), read_file(fileBy, conversionB)) :
                                 (read_file(fileBy, conversionB), read_file(fileBz, conversionB), read_file(fileBx, conversionB))

    B1, B2, BLOS = permute_dims(B1, LOS), permute_dims(B2, LOS), permute_dims(BLOS, LOS)
    T, n = permute_dims(T, LOS), permute_dims(n, LOS)

    return (B1, B2, BLOS, T, n, nH2, nHp)
end
