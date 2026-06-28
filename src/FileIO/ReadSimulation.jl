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

"""
    los_basis(Ax, Ay, Az, LOS) -> (A1, A2, ALOS)

Map cartesian vector components onto the right-handed (plane-of-sky ⊕ LOS)
basis using cyclic permutations (determinant +1, chirality preserved):
LOS = "z" → (Ax, Ay, Az) ; LOS = "x" → (Ay, Az, Ax) ; LOS = "y" → (Az, Ax, Ay).

This is the single source of truth for the LOS frame convention: the intrinsic
polarization angle ψ_src = atan(A2, A1) + π/2 is only consistent across the
three lines of sight if the same orientation-preserving mapping is used
everywhere.
"""
function los_basis(Ax, Ay, Az, LOS::AbstractString)
    LOS == "z" && return (Ax, Ay, Az)
    LOS == "x" && return (Ay, Az, Ax)
    LOS == "y" && return (Az, Ax, Ay)
    error("Unknown LOS: $LOS (expected \"x\", \"y\" or \"z\")")
end

function validate_required_fits(file)
    validation_error = ensure_readable_file(file; expected_exts=[".fits"])
    validation_error === nothing || error(validation_error)
end

const REQUIRED_SIMULATION_FIELDS = ("Bx", "By", "Bz", "density", "temperature")

function _fits_file_candidates(dir::AbstractString)
    isdir(dir) || return String[]
    files = String[]
    for name in sort(readdir(dir))
        ext = lowercase(splitext(name)[2])
        ext in (".fits", ".fit", ".fts") && push!(files, joinpath(dir, name))
    end
    return files
end

function simulation_field_source(simu::AbstractString, field::AbstractString)
    file_path = joinpath(simu, "$(field).fits")
    isfile(file_path) && return file_path

    dir_path = joinpath(simu, field)
    dir_files = _fits_file_candidates(dir_path)
    !isempty(dir_files) && return dir_files

    prefix = field * "_"
    prefixed_files = String[]
    for name in sort(readdir(simu))
        ext = lowercase(splitext(name)[2])
        if startswith(name, prefix) && ext in (".fits", ".fit", ".fts")
            push!(prefixed_files, joinpath(simu, name))
        end
    end
    !isempty(prefixed_files) && return prefixed_files

    return file_path
end

function _validate_fits_source(source)
    if source isa AbstractString
        validate_required_fits(source)
    else
        isempty(source) && error("Cannot read an empty FITS source list.")
        foreach(validate_required_fits, source)
    end
end

function _source_grid_kind(source)
    if source isa AbstractString
        return detect_fits_grid(source)
    end

    kinds = detect_fits_grid.(source)
    length(unique(kinds)) == 1 || error("Cannot mix HEALPix FITS tables and regular FITS images for one simulation field.")
    return first(kinds)
end

function simulation_grid_kind(simu::AbstractString)
    validation_error = ensure_directory_access(simu)
    validation_error === nothing || error(validation_error)

    sources = [simulation_field_source(simu, field) for field in REQUIRED_SIMULATION_FIELDS]
    foreach(_validate_fits_source, sources)
    kinds = _source_grid_kind.(sources)
    length(unique(kinds)) == 1 ||
        error("Simulation $(simu) mixes regular FITS images and HEALPix FITS tables. Use one grid type per run.")

    return first(kinds)
end

function healpix_simulation_metadata(simu::AbstractString)
    simulation_grid_kind(simu) == :healpix ||
        error("Simulation $(simu) is not backed by HEALPix FITS inputs.")
    stack = read_fits_grid_stack(simulation_field_source(simu, "Bx"); column=:all)
    stack isa HealpixStack || error("Expected HEALPix stack for $(joinpath(simu, "Bx")).")
    return (; nside = stack.nside, order = stack.order)
end

function read_required_cube(file, conversion)
    validate_required_fits(file)

    try
        return read_file(file, conversion; expected_ndims=3)
    catch err
        error(user_error_message(:corrupted_file, file; reason=string(err)))
    end
end

function _grid_to_processing_cube(grid)
    if grid isa HealpixStack
        return reshape(grid.pixels, size(grid.pixels, 1), 1, size(grid.pixels, 2))
    elseif grid isa AbstractVector && hasproperty(grid, :resolution)
        return reshape(collect(grid), length(grid), 1, 1)
    elseif grid isa AbstractArray
        return grid
    end

    error("Unsupported simulation grid input $(typeof(grid)).")
end

function read_required_grid(source, conversion)
    _validate_fits_source(source)

    try
        grid = read_fits_grid_stack(source, conversion; column=:all, expected_ndims = nothing)
        if grid isa AbstractArray
            _validate_fits_array(grid, source isa AbstractString ? source : first(source); expected_ndims = 3)
        end
        return _grid_to_processing_cube(grid)
    catch err
        label = source isa AbstractString ? source : join(source, ", ")
        error(user_error_message(:corrupted_file, label; reason=string(err)))
    end
end

function read_optional_cube(file, conversion, LOS)
    if !ispath(file) || !isfile(file)
        return nothing
    end

    validation_error = ensure_readable_file(file; expected_exts=[".fits"])
    validation_error === nothing || error(validation_error)

    try
        return permute_dims(read_file(file, conversion; expected_ndims=3), LOS)
    catch err
        error(user_error_message(:corrupted_file, file; reason=string(err)))
    end
end

function read_optional_grid(source, conversion, LOS, grid_kind)
    exists = source isa AbstractString ? isfile(source) : !isempty(source)
    exists || return nothing

    if grid_kind == :healpix
        return read_required_grid(source, conversion)
    end

    source isa AbstractString || error("Regular cube input for optional fields must be a single FITS image.")
    return read_optional_cube(source, conversion, LOS)
end

function ReadSimulation(simu, LOS, conversionn, conversionT, conversionV, conversionB)
    validation_error = ensure_directory_access(simu)
    validation_error === nothing || error(validation_error)

    fileBx = simulation_field_source(simu, "Bx")
    fileBy = simulation_field_source(simu, "By")
    fileBz = simulation_field_source(simu, "Bz")
    filen = simulation_field_source(simu, "density")
    fileT = simulation_field_source(simu, "temperature")
    fileVx = simulation_field_source(simu, "Vx")
    fileVy = simulation_field_source(simu, "Vy")
    fileVz = simulation_field_source(simu, "Vz")
    filenH2 = simulation_field_source(simu, "densityH2")
    filenHp = simulation_field_source(simu, "densityHp")
    grid_kind = simulation_grid_kind(simu)

    required_files = [fileBx, fileBy, fileBz, filen, fileT, fileVx, fileVy, fileVz]
    foreach(_validate_fits_source, required_files)

    T = read_required_grid(fileT, conversionT)
    n = read_required_grid(filen, conversionn)
    nH2 = read_optional_grid(filenH2, conversionn, LOS, grid_kind)
    nHp = read_optional_grid(filenHp, conversionn, LOS, grid_kind)

    B1, B2, BLOS = los_basis(read_required_grid(fileBx, conversionB),
                             read_required_grid(fileBy, conversionB),
                             read_required_grid(fileBz, conversionB), LOS)

    V1, V2, VLOS = los_basis(read_required_grid(fileVx, conversionV),
                             read_required_grid(fileVy, conversionV),
                             read_required_grid(fileVz, conversionV), LOS)

    if grid_kind == :image
        B1, B2, BLOS = permute_dims(B1, LOS), permute_dims(B2, LOS), permute_dims(BLOS, LOS)
        V1, V2, VLOS = permute_dims(V1, LOS), permute_dims(V2, LOS), permute_dims(VLOS, LOS)
        T, n = permute_dims(T, LOS), permute_dims(n, LOS)
    end

    return (B1, B2, BLOS, V1, V2, VLOS, T, n, nH2, nHp)
end

function ReadSimulation(simu, LOS, conversionn, conversionT, conversionB)
    validation_error = ensure_directory_access(simu)
    validation_error === nothing || error(validation_error)

    fileBx = simulation_field_source(simu, "Bx")
    fileBy = simulation_field_source(simu, "By")
    fileBz = simulation_field_source(simu, "Bz")
    filen = simulation_field_source(simu, "density")
    fileT = simulation_field_source(simu, "temperature")
    filenH2 = simulation_field_source(simu, "densityH2")
    filenHp = simulation_field_source(simu, "densityHp")
    grid_kind = simulation_grid_kind(simu)

    required_files = [fileBx, fileBy, fileBz, filen, fileT]
    foreach(_validate_fits_source, required_files)

    T = read_required_grid(fileT, conversionT)
    n = read_required_grid(filen, conversionn)
    nH2 = read_optional_grid(filenH2, conversionn, LOS, grid_kind)
    nHp = read_optional_grid(filenHp, conversionn, LOS, grid_kind)

    B1, B2, BLOS = los_basis(read_required_grid(fileBx, conversionB),
                             read_required_grid(fileBy, conversionB),
                             read_required_grid(fileBz, conversionB), LOS)

    if grid_kind == :image
        B1, B2, BLOS = permute_dims(B1, LOS), permute_dims(B2, LOS), permute_dims(BLOS, LOS)
        T, n = permute_dims(T, LOS), permute_dims(n, LOS)
    end

    return (B1, B2, BLOS, T, n, nH2, nHp)
end
