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

function validate_required_fits(file)
    validation_error = ensure_readable_file(file; expected_exts=[".fits"])
    validation_error === nothing || error(validation_error)
end

function read_required_cube(file, conversion)
    validate_required_fits(file)

    try
        return read_file(file, conversion)
    catch err
        error(user_error_message(:corrupted_file, file; reason=string(err)))
    end
end

function read_optional_cube(file, conversion, LOS)
    if !ispath(file) || !isfile(file)
        return nothing
    end

    validation_error = ensure_readable_file(file; expected_exts=[".fits"])
    validation_error === nothing || error(validation_error)

    try
        return permute_dims(read_file(file, conversion), LOS)
    catch err
        error(user_error_message(:corrupted_file, file; reason=string(err)))
    end
end

function ReadSimulation(simu, LOS, conversionn, conversionT, conversionV, conversionB)
    validation_error = ensure_directory_access(simu)
    validation_error === nothing || error(validation_error)

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

    required_files = [fileBx, fileBy, fileBz, filen, fileT, fileVx, fileVy, fileVz]
    foreach(validate_required_fits, required_files)

    T = read_required_cube(fileT, conversionT)
    n = read_required_cube(filen, conversionn)
    nH2 = read_optional_cube(filenH2, conversionn, LOS)
    nHp = read_optional_cube(filenHp, conversionn, LOS)

    B1, B2, BLOS = LOS == "z" ? (read_required_cube(fileBx, conversionB), read_required_cube(fileBy, conversionB), read_required_cube(fileBz, conversionB)) :
                    LOS == "y" ? (read_required_cube(fileBz, conversionB), read_required_cube(fileBx, conversionB), read_required_cube(fileBy, conversionB)) :
                                 (read_required_cube(fileBy, conversionB), read_required_cube(fileBz, conversionB), read_required_cube(fileBx, conversionB))

    V1, V2, VLOS = LOS == "z" ? (read_required_cube(fileVx, conversionV), read_required_cube(fileVy, conversionV), read_required_cube(fileVz, conversionV)) :
                    LOS == "y" ? (read_required_cube(fileVz, conversionV), read_required_cube(fileVx, conversionV), read_required_cube(fileVy, conversionV)) :
                                 (read_required_cube(fileVy, conversionV), read_required_cube(fileVz, conversionV), read_required_cube(fileVx, conversionV))

    B1, B2, BLOS = permute_dims(B1, LOS), permute_dims(B2, LOS), permute_dims(BLOS, LOS)
    V1, V2, VLOS = permute_dims(V1, LOS), permute_dims(V2, LOS), permute_dims(VLOS, LOS)
    T, n = permute_dims(T, LOS), permute_dims(n, LOS)

    return (B1, B2, BLOS, V1, V2, VLOS, T, n, nH2, nHp)
end

function ReadSimulation(simu, LOS, conversionn, conversionT, conversionB)
    validation_error = ensure_directory_access(simu)
    validation_error === nothing || error(validation_error)

    fileBx = "$simu/Bx.fits"
    fileBy = "$simu/By.fits"
    fileBz = "$simu/Bz.fits"
    filen = "$simu/density.fits"
    fileT = "$simu/temperature.fits"
    filenH2 = "$simu/densityH2.fits"
    filenHp = "$simu/densityHp.fits"

    required_files = [fileBx, fileBy, fileBz, filen, fileT]
    foreach(validate_required_fits, required_files)

    T = read_required_cube(fileT, conversionT)
    n = read_required_cube(filen, conversionn)
    nH2 = read_optional_cube(filenH2, conversionn, LOS)
    nHp = read_optional_cube(filenHp, conversionn, LOS)

    B1, B2, BLOS = LOS == "z" ? (read_required_cube(fileBx, conversionB), read_required_cube(fileBy, conversionB), read_required_cube(fileBz, conversionB)) :
                    LOS == "y" ? (read_required_cube(fileBx, conversionB), read_required_cube(fileBz, conversionB), read_required_cube(fileBy, conversionB)) :
                                 (read_required_cube(fileBy, conversionB), read_required_cube(fileBz, conversionB), read_required_cube(fileBx, conversionB))

    B1, B2, BLOS = permute_dims(B1, LOS), permute_dims(B2, LOS), permute_dims(BLOS, LOS)
    T, n = permute_dims(T, LOS), permute_dims(n, LOS)

    return (B1, B2, BLOS, T, n, nH2, nHp)
end

