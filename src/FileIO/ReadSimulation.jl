"""
    ReadSimulation(simu::String, LOS::String, conversionn::Number, conversionT::Number, conversionV::Number, conversionB::Number)
        -> Tuple{AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, Union{AbstractArray, Nothing}, Union{AbstractArray, Nothing}}

Reads and processes simulation data from FITS or HDF5 files for a specified line of sight (LOS).

# Arguments
- `simu::String`: The directory containing the simulation files.
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

const FITS_EXTS = (".fits", ".fit", ".fts")
const SIMULATION_FIELD_EXTS = (FITS_EXTS..., HDF5_EXTS...)

function _is_fits_path(path::AbstractString)
    return lowercase(splitext(path)[2]) in FITS_EXTS
end

function validate_required_fits(file)
    validation_error = ensure_readable_file(file; expected_exts=collect(FITS_EXTS))
    validation_error === nothing || error(validation_error)
end

const REQUIRED_SIMULATION_FIELDS = ("Bx", "By", "Bz", "density", "temperature")

function _fits_file_candidates(dir::AbstractString)
    isdir(dir) || return String[]
    files = String[]
    for name in sort(readdir(dir))
        ext = lowercase(splitext(name)[2])
        ext in FITS_EXTS && push!(files, joinpath(dir, name))
    end
    return files
end

function _field_file_candidates(dir::AbstractString, field::AbstractString)
    isdir(dir) || return String[]
    files = String[]
    for name in sort(readdir(dir))
        stem, ext = splitext(name)
        lowercase(ext) in SIMULATION_FIELD_EXTS && stem == field && push!(files, joinpath(dir, name))
    end
    return files
end

function _shared_hdf5_field_source(simu::AbstractString, field::AbstractString)
    for name in sort(readdir(simu))
        file = joinpath(simu, name)
        isfile(file) && is_hdf5_path(file) || continue
        dataset = find_hdf5_dataset(file, field)
        dataset === nothing || return HDF5DatasetSource(file, dataset)
    end

    return nothing
end

function simulation_field_source(simu::AbstractString, field::AbstractString)
    exact_files = _field_file_candidates(simu, field)
    !isempty(exact_files) && return first(exact_files)

    shared_source = _shared_hdf5_field_source(simu, field)
    shared_source === nothing || return shared_source

    dir_path = joinpath(simu, field)
    dir_files = _fits_file_candidates(dir_path)
    !isempty(dir_files) && return dir_files

    prefix = field * "_"
    prefixed_files = String[]
    for name in sort(readdir(simu))
        ext = lowercase(splitext(name)[2])
        if startswith(name, prefix) && ext in FITS_EXTS
            push!(prefixed_files, joinpath(simu, name))
        end
    end
    !isempty(prefixed_files) && return prefixed_files

    return joinpath(simu, "$(field).fits")
end

function _validate_simulation_source(source)
    if source isa HDF5DatasetSource
        validation_error = ensure_readable_file(source.file; expected_exts=collect(HDF5_EXTS))
        validation_error === nothing || error(validation_error)
        find_hdf5_dataset(source.file, source.dataset) !== nothing ||
            error("Dataset $(source.dataset) was not found in HDF5 file $(source.file).")
    elseif source isa AbstractString
        if is_hdf5_path(source)
            validation_error = ensure_readable_file(source; expected_exts=collect(HDF5_EXTS))
            validation_error === nothing || error(validation_error)
        else
            validate_required_fits(source)
        end
    else
        isempty(source) && error("Cannot read an empty simulation source list.")
        foreach(_validate_simulation_source, source)
    end
end

function _validate_fits_source(source)
    if source isa AbstractString
        validate_required_fits(source)
    else
        isempty(source) && error("Cannot read an empty FITS source list.")
        foreach(validate_required_fits, source)
    end
end

function _single_source_grid_kind(source)
    source isa HDF5DatasetSource && return :image
    if source isa AbstractString
        is_hdf5_path(source) && return :image
        return detect_fits_grid(source)
    end
    error("Unsupported simulation field source $(typeof(source)).")
end

function _source_grid_kind(source)
    if source isa AbstractString || source isa HDF5DatasetSource
        return _single_source_grid_kind(source)
    end

    kinds = _single_source_grid_kind.(source)
    length(unique(kinds)) == 1 || error("Cannot mix HEALPix FITS tables and regular FITS images for one simulation field.")
    return first(kinds)
end

function simulation_grid_kind(simu::AbstractString)
    validation_error = ensure_directory_access(simu)
    validation_error === nothing || error(validation_error)

    sources = [simulation_field_source(simu, field) for field in REQUIRED_SIMULATION_FIELDS]
    foreach(_validate_simulation_source, sources)
    kinds = _source_grid_kind.(sources)
    length(unique(kinds)) == 1 ||
        error("Simulation $(simu) mixes regular image cubes and HEALPix FITS tables. Use one grid type per run.")

    return first(kinds)
end

"""
    _validate_healpix_los(LOS)

HEALPix simulations are radial by construction: each pixel's line of sight is
the local radial direction `e_r`, so the cartesian LOS permutations `x`/`y`
are meaningless on the sphere and are rejected.

Required field convention for HEALPix inputs (documented, not verifiable from
the data): `Bx`/`Vx` must hold the component along the local colatitude
direction `e_θ`, `By`/`Vy` the component along the local azimuth direction
`e_φ`, and `Bz`/`Vz` the radial (LOS) component along `e_r`, for every pixel.
With `LOS = "z"` these map unchanged onto `(B1, B2, BLOS)`, and the intrinsic
polarization angle `ψ = atan(B2, B1) + π/2` is measured in the local
`(e_θ, e_φ)` basis — the same formula as the cartesian path. Global cartesian
components would make the polarization angle wrong everywhere except at the
pole; project them onto the per-pixel tangent basis before writing the FITS
inputs.
"""
function _validate_healpix_los(LOS)
    String(LOS) == "z" || throw_config_error(
        "HEALPix simulations only support LOS=\"z\" (got \"$(LOS)\"): the line of sight is the local radial direction, " *
        "and Bx/By/Bz must already be the per-pixel (e_theta, e_phi, e_r) tangent components. " *
        "Cartesian LOS permutations \"x\"/\"y\" are meaningless on the sphere.";
        code=:unsupported_grid_operation)
    return nothing
end

# NSIDE/ordering/coordsys of the reference field (Bx), read from the FITS
# headers only (no pixel data is loaded).
function _healpix_target_from_source(source)
    file = source_path(source isa AbstractVector ? first(source) : source)
    info = _healpix_table_hdu_info(file)
    return (; nside = info.nside, order = info.order, coordsys = info.coordsys)
end

function healpix_simulation_metadata(simu::AbstractString)
    simulation_grid_kind(simu) == :healpix ||
        error("Simulation $(simu) is not backed by HEALPix FITS inputs.")
    return _healpix_target_from_source(simulation_field_source(simu, "Bx"))
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

function read_grid_stack(files::AbstractVector, conversion::Real=1.0;
    column=1,
    T::Type=Float64,
    expected_ndims=nothing,
    allow_nonfinite::Bool=false,
    unseen_to_nan::Bool=true)

    isempty(files) && error("Cannot read an empty simulation stack.")
    kinds = _single_source_grid_kind.(files)
    length(unique(kinds)) == 1 ||
        error("Cannot mix HEALPix FITS tables and regular image files in one stack.")

    if first(kinds) == :healpix
        all(source -> source isa AbstractString && _is_fits_path(source), files) ||
            error("HEALPix stacks must be FITS files.")
        return read_fits_grid_stack(files, conversion; column=column, T=T, expected_ndims=expected_ndims,
            allow_nonfinite=allow_nonfinite, unseen_to_nan=unseen_to_nan)
    end

    if length(files) == 1
        return read_file(first(files), conversion; expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite)
    end

    planes = [read_file(file, conversion; expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite) for file in files]
    first_size = size(first(planes))
    all(size(plane) == first_size for plane in planes) ||
        error("All regular image planes in a stack must have the same shape.")

    return cat(planes...; dims=ndims(first(planes)) + 1)
end

function read_grid_stack(source, conversion::Real=1.0; kwargs...)
    source isa AbstractVector && return read_grid_stack(source, conversion; kwargs...)
    return read_grid_stack([source], conversion; kwargs...)
end

function read_required_grid(source, conversion; hp_target=nothing)
    _validate_simulation_source(source)

    try
        grid = read_grid_stack(source, conversion; column=:all, expected_ndims = nothing)
        if grid isa AbstractArray
            _validate_fits_array(grid, source_label(source isa AbstractVector ? first(source) : source); expected_ndims = 3)
        end
        if hp_target !== nothing && grid isa HealpixStack
            label = source_label(source isa AbstractVector ? first(source) : source)
            grid = _conform_healpix_stack(grid, hp_target.nside, hp_target.order; label = label)
        end
        return _grid_to_processing_cube(grid)
    catch err
        label = source isa AbstractVector ? join(source_label.(source), ", ") : source_label(source)
        error(user_error_message(:corrupted_file, label; reason=string(err)))
    end
end

function read_optional_cube(file, conversion, LOS)
    if !ispath(file) || !isfile(file)
        return nothing
    end

    validation_error = ensure_readable_file(file; expected_exts=collect(FITS_EXTS))
    validation_error === nothing || error(validation_error)

    try
        return permute_dims(read_file(file, conversion; expected_ndims=3), LOS)
    catch err
        error(user_error_message(:corrupted_file, file; reason=string(err)))
    end
end

function read_optional_grid(source, conversion, LOS, grid_kind; hp_target=nothing)
    exists = source isa HDF5DatasetSource ? isfile(source.file) : source isa AbstractString ? isfile(source) : !isempty(source)
    exists || return nothing

    if grid_kind == :healpix
        return read_required_grid(source, conversion; hp_target = hp_target)
    end

    if source isa AbstractString && _is_fits_path(source)
        return read_optional_cube(source, conversion, LOS)
    end

    return permute_dims(read_required_grid(source, conversion), LOS)
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
    grid_kind == :healpix && _validate_healpix_los(LOS)
    hp_target = grid_kind == :healpix ? _healpix_target_from_source(fileBx) : nothing

    required_files = [fileBx, fileBy, fileBz, filen, fileT, fileVx, fileVy, fileVz]
    foreach(_validate_simulation_source, required_files)

    T = read_required_grid(fileT, conversionT; hp_target = hp_target)
    n = read_required_grid(filen, conversionn; hp_target = hp_target)
    nH2 = read_optional_grid(filenH2, conversionn, LOS, grid_kind; hp_target = hp_target)
    nHp = read_optional_grid(filenHp, conversionn, LOS, grid_kind; hp_target = hp_target)

    B1, B2, BLOS = los_basis(read_required_grid(fileBx, conversionB; hp_target = hp_target),
                             read_required_grid(fileBy, conversionB; hp_target = hp_target),
                             read_required_grid(fileBz, conversionB; hp_target = hp_target), LOS)

    V1, V2, VLOS = los_basis(read_required_grid(fileVx, conversionV; hp_target = hp_target),
                             read_required_grid(fileVy, conversionV; hp_target = hp_target),
                             read_required_grid(fileVz, conversionV; hp_target = hp_target), LOS)

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
    grid_kind == :healpix && _validate_healpix_los(LOS)
    hp_target = grid_kind == :healpix ? _healpix_target_from_source(fileBx) : nothing

    required_files = [fileBx, fileBy, fileBz, filen, fileT]
    foreach(_validate_simulation_source, required_files)

    T = read_required_grid(fileT, conversionT; hp_target = hp_target)
    n = read_required_grid(filen, conversionn; hp_target = hp_target)
    nH2 = read_optional_grid(filenH2, conversionn, LOS, grid_kind; hp_target = hp_target)
    nHp = read_optional_grid(filenHp, conversionn, LOS, grid_kind; hp_target = hp_target)

    B1, B2, BLOS = los_basis(read_required_grid(fileBx, conversionB; hp_target = hp_target),
                             read_required_grid(fileBy, conversionB; hp_target = hp_target),
                             read_required_grid(fileBz, conversionB; hp_target = hp_target), LOS)

    if grid_kind == :image
        B1, B2, BLOS = permute_dims(B1, LOS), permute_dims(B2, LOS), permute_dims(BLOS, LOS)
        T, n = permute_dims(T, LOS), permute_dims(n, LOS)
    end

    return (B1, B2, BLOS, T, n, nH2, nHp)
end
