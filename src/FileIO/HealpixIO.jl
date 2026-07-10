"""
Helpers for using MOOSE algorithms with HEALPix maps.

The core MOOSE RM-synthesis routine works on arrays whose last dimension is
frequency.  These helpers bridge Healpix.jl maps to that representation while
preserving NSIDE, ordering and coordinate-system (`COORDSYS`) metadata for
writing results back to HEALPix FITS.

Masked pixels follow the HEALPix "UNSEEN" convention (sentinel value
`-1.6375e30`): they are converted to `NaN` when reading stacks so they
propagate correctly through MOOSE algorithms, and converted back to the
sentinel when writing FITS output.

Vector-field convention for HEALPix *simulation inputs* (synchrotron
pipeline): the line of sight is the local radial direction `e_r` of each
pixel, so only `LOS = "z"` is accepted, and the `Bx`/`By`/`Bz` (and
`Vx`/`Vy`/`Vz`) files must hold the per-pixel tangent-basis components
`(B·e_θ, B·e_φ, B·e_r)` — NOT global cartesian components. See
`_validate_healpix_los` in `ReadSimulation.jl`.
"""

const HealpixOrderName = Union{Symbol, AbstractString}
const FITSGridKind = Symbol

"""
Sentinel value used by the HEALPix convention to mark masked ("bad") pixels.
"""
const HEALPIX_UNSEEN = -1.6375e30

# Threshold comparison instead of equality: the sentinel loses precision when
# round-tripped through Float32 FITS columns.
_is_unseen(x::Real) = isfinite(x) && x <= -1.63e30
_is_unseen(x) = false

function _mask_unseen!(pixels::AbstractArray{T}) where {T}
    T <: AbstractFloat || return pixels
    @inbounds for i in eachindex(pixels)
        _is_unseen(pixels[i]) && (pixels[i] = T(NaN))
    end
    return pixels
end

function _restore_unseen(pixels::AbstractArray{T}) where {T}
    T <: AbstractFloat || return pixels
    any(isnan, pixels) || return pixels
    out = copy(pixels)
    @inbounds for i in eachindex(out)
        isnan(out[i]) && (out[i] = T(HEALPIX_UNSEEN))
    end
    return out
end

_as_matrix(x::Matrix) = x
_as_matrix(x::AbstractMatrix) = Matrix(x)

function _normalize_coordsys(coordsys)
    coordsys === nothing && return nothing
    value = uppercase(strip(String(coordsys)))
    return isempty(value) ? nothing : value
end

struct HealpixStack{T}
    pixels::Matrix{T}
    nside::Int
    order::Symbol
    coordsys::Union{Nothing, String}
end

HealpixStack(pixels::Matrix, nside::Integer, order::Symbol) =
    HealpixStack(pixels, Int(nside), order, nothing)

struct HealpixRMResult{T}
    fdf::Matrix{T}
    realFDF::Matrix{T}
    imagFDF::Matrix{T}
    phi::Vector{Float64}
    nside::Int
    order::Symbol
    coordsys::Union{Nothing, String}
end

HealpixRMResult(fdf::Matrix, realFDF::Matrix, imagFDF::Matrix,
    phi::Vector{Float64}, nside::Integer, order::Symbol) =
    HealpixRMResult(fdf, realFDF, imagFDF, phi, Int(nside), order, nothing)

function _fits_header_value(header, key::AbstractString, default=nothing)
    haskey(header, key) || return default
    return header[key]
end

function _parse_header_int(value)
    value === nothing && return nothing
    value isa Integer && return Int(value)
    return tryparse(Int, strip(string(value)))
end

"""
    detect_fits_grid(filename) -> Symbol

Detect whether a FITS file stores a HEALPix map (`:healpix`) or a regular FITS
image (`:image`) by inspecting HDU headers. HEALPix FITS files are recognised
from the standard `PIXTYPE=HEALPIX` metadata, with `ORDERING`/`NSIDE` as a
fallback for older files.
"""
function detect_fits_grid(filename::AbstractString)::FITSGridKind
    validation_error = ensure_readable_file(filename; expected_exts=[".fits", ".fit", ".fts"])
    validation_error === nothing || error(validation_error)

    detected = FITS(filename) do fits
        saw_image = false
        for hdu in fits
            header = read_header(hdu)
            pixtype = uppercase(strip(String(_fits_header_value(header, "PIXTYPE", ""))))
            ordering = _fits_header_value(header, "ORDERING", nothing)
            nside = _fits_header_value(header, "NSIDE", nothing)
            xtension = uppercase(strip(String(_fits_header_value(header, "XTENSION", ""))))
            naxis = _fits_header_value(header, "NAXIS", 0)

            if pixtype == "HEALPIX" || (ordering !== nothing && nside !== nothing && xtension == "BINTABLE")
                return :healpix
            end

            parsed_naxis = _parse_header_int(naxis)
            if parsed_naxis !== nothing && parsed_naxis > 0 && xtension != "BINTABLE"
                saw_image = true
            end
        end

        return saw_image ? :image : nothing
    end

    detected !== nothing && return detected
    error("Could not detect whether $(filename) is a HEALPix FITS table or a regular FITS image.")
end

is_healpix_fits(filename::AbstractString) = detect_fits_grid(filename) == :healpix
is_image_fits(filename::AbstractString) = detect_fits_grid(filename) == :image

Base.size(stack::HealpixStack) = size(stack.pixels)
Base.size(stack::HealpixStack, dim::Integer) = size(stack.pixels, dim)
Base.getindex(stack::HealpixStack, args...) = getindex(stack.pixels, args...)

function normalize_healpix_order(order::HealpixOrderName)
    normalized = lowercase(strip(String(order)))
    if normalized in ("ring", "ringorder")
        return :ring
    elseif normalized in ("nested", "nest", "nestedorder")
        return :nested
    end

    error("Invalid HEALPix ordering $(order). Use :ring or :nested.")
end

healpix_order_type(order::HealpixOrderName) =
    normalize_healpix_order(order) == :nested ? Healpix.NestedOrder : Healpix.RingOrder

function healpix_nside_from_npix(npix::Integer)
    return Healpix.npix2nside(npix)
end

function healpix_map(pixels::AbstractVector; nside::Union{Nothing, Integer}=nothing, order::HealpixOrderName=:ring)
    inferred_nside = healpix_nside_from_npix(length(pixels))
    if nside !== nothing && Int(nside) != inferred_nside
        error("Pixel vector length $(length(pixels)) corresponds to NSIDE=$(inferred_nside), not NSIDE=$(nside).")
    end

    order_type = healpix_order_type(order)
    return Healpix.HealpixMap{eltype(pixels), order_type}(collect(pixels))
end

function healpix_map(nside::Integer; T::Type=Float64, order::HealpixOrderName=:ring)
    order_type = healpix_order_type(order)
    return Healpix.HealpixMap{T, order_type}(nside)
end

healpix_order(::Healpix.HealpixMap{T, Healpix.RingOrder}) where {T} = :ring
healpix_order(::Healpix.HealpixMap{T, Healpix.NestedOrder}) where {T} = :nested

# Fallback for duck-typed map-like objects that are not Healpix.HealpixMap.
function healpix_order(map)
    type_name = string(typeof(map))
    occursin("NestedOrder", type_name) && return :nested
    occursin("RingOrder", type_name) && return :ring
    error("Could not infer HEALPix ordering from map type $(typeof(map)).")
end

healpix_nside(map) = Int(map.resolution.nside)

"""
    read_healpix_map(filename; column=1, T=Float64) -> Healpix.HealpixMap

Thin wrapper around `Healpix.readMapFromFITS`. Note that UNSEEN sentinel
values are kept as-is here (the Healpix.jl ecosystem understands them);
[`read_healpix_stack`](@ref) converts them to `NaN` for MOOSE processing.
"""
function read_healpix_map(filename::AbstractString; column=1, T::Type=Float64)
    validation_error = ensure_readable_file(filename; expected_exts=[".fits", ".fit", ".fts"])
    validation_error === nothing || error(validation_error)

    return Healpix.readMapFromFITS(String(filename), column, T)
end

function _healpix_table_hdu_info(fits::FITS)
    for idx in 1:length(fits)
        hdu = fits[idx]
        header = read_header(hdu)
        pixtype = uppercase(strip(String(_fits_header_value(header, "PIXTYPE", ""))))
        ordering = _fits_header_value(header, "ORDERING", nothing)
        nside = _fits_header_value(header, "NSIDE", nothing)
        xtension = uppercase(strip(String(_fits_header_value(header, "XTENSION", ""))))

        if pixtype == "HEALPIX" || (ordering !== nothing && nside !== nothing && xtension == "BINTABLE")
            ordering === nothing && error("HEALPix HDU $(idx) lacks the ORDERING keyword.")
            nside_value = _parse_header_int(nside)
            nside_value === nothing && error("HEALPix HDU $(idx) has an invalid NSIDE header value $(repr(nside)).")

            return (;
                hdu_index = idx,
                nside = nside_value,
                order = normalize_healpix_order(String(ordering)),
                coordsys = _normalize_coordsys(_fits_header_value(header, "COORDSYS", nothing)),
                colnames = FITSIO.colnames(hdu),
                nrows = something(_parse_header_int(_fits_header_value(header, "NAXIS2", 0)), 0),
            )
        end
    end

    return nothing
end

function _healpix_table_hdu_info(filename::AbstractString)
    validation_error = ensure_readable_file(filename; expected_exts=[".fits", ".fit", ".fts"])
    validation_error === nothing || error(validation_error)

    info = FITS(filename) do fits
        _healpix_table_hdu_info(fits)
    end

    info === nothing && error("No HEALPix binary table HDU found in $(filename).")
    return info
end

function _normalize_healpix_columns(column)
    column === :all && return :all
    if column isa Integer || column isa AbstractString
        return [column]
    elseif column isa AbstractVector || column isa Tuple
        return collect(column)
    end

    error("Invalid HEALPix column selector $(column). Use an integer, name, list, tuple, or :all.")
end

function _healpix_matrix_from_column_data(data::AbstractArray, info)
    npix = Healpix.nside2npix(info.nside)
    if ndims(data) == 1
        length(data) == npix || error("HEALPix column has $(length(data)) values, expected $(npix) for NSIDE=$(info.nside).")
        return reshape(collect(data), npix, 1)
    elseif ndims(data) == 2
        if size(data, 1) == npix
            return _as_matrix(data)
        elseif size(data, 2) == npix
            return permutedims(data)
        elseif length(data) % npix == 0
            @warn "Ambiguous HEALPix column shape $(size(data)) for NSIDE=$(info.nside); assuming column-major Npix x Nslice layout." npix
            return reshape(vec(data), npix, length(data) ÷ npix)
        end
    end

    error("Cannot interpret HEALPix table column with size $(size(data)) as an NSIDE=$(info.nside) map or cube.")
end

"""
    read_healpix_stack(filename; column=:all, T=Float64, unseen_to_nan=true) -> HealpixStack

Read one or more HEALPix table columns from `filename` into an `Npix × Nslice`
[`HealpixStack`](@ref). UNSEEN sentinel pixels are converted to `NaN` unless
`unseen_to_nan=false`.
"""
function read_healpix_stack(filename::AbstractString; column=:all, T::Type=Float64, unseen_to_nan::Bool=true)
    validation_error = ensure_readable_file(filename; expected_exts=[".fits", ".fit", ".fts"])
    validation_error === nothing || error(validation_error)

    stack = FITS(filename) do fits
        info = _healpix_table_hdu_info(fits)
        info === nothing && error("No HEALPix binary table HDU found in $(filename).")

        columns = _normalize_healpix_columns(column)
        columns = columns === :all ? collect(eachindex(info.colnames)) : columns
        isempty(columns) && error("Cannot read an empty HEALPix column selection.")

        hdu = fits[info.hdu_index]
        chunks = Matrix{T}[]
        for col in columns
            name = col isa Integer ? info.colnames[Int(col)] : String(col)
            data = T.(read(hdu, name; case_sensitive=false))
            push!(chunks, _healpix_matrix_from_column_data(data, info))
        end

        pixels = length(chunks) == 1 ? chunks[1] : reduce(hcat, chunks)
        HealpixStack(pixels; nside=info.nside, order=info.order, coordsys=info.coordsys)
    end

    unseen_to_nan && _mask_unseen!(stack.pixels)
    return stack
end

function read_healpix_stack(files::AbstractVector{<:AbstractString}; column=1, T::Type=Float64, unseen_to_nan::Bool=true)
    isempty(files) && error("Cannot read an empty HEALPix stack.")
    stack = HealpixStack([read_healpix_map(file; column=column, T=T) for file in files])
    unseen_to_nan && _mask_unseen!(stack.pixels)
    return stack
end

function read_fits_grid(filename::AbstractString, conversion::Real=1.0;
    column=1,
    T::Type=Float64,
    expected_ndims=nothing,
    allow_nonfinite::Bool=false,
    unseen_to_nan::Bool=true)

    kind = detect_fits_grid(filename)
    if kind == :healpix
        expected_ndims === nothing ||
            error("$(filename) is a HEALPix FITS table, not a $(expected_ndims)D image cube.")
        stack = read_healpix_stack(filename; column=column, T=T, unseen_to_nan=unseen_to_nan)
        if size(stack.pixels, 2) == 1
            map = healpix_map(view(stack.pixels, :, 1); nside=stack.nside, order=stack.order)
            return conversion == 1 ? map : healpix_map(collect(map) .* conversion; nside=stack.nside, order=stack.order)
        end

        return conversion == 1 ? stack :
            HealpixStack(stack.pixels .* conversion; nside=stack.nside, order=stack.order, coordsys=stack.coordsys)
    end

    return read_file(filename, conversion; expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite)
end

function read_fits_grid_stack(files::AbstractVector{<:AbstractString}, conversion::Real=1.0;
    column=1,
    T::Type=Float64,
    expected_ndims=nothing,
    allow_nonfinite::Bool=false,
    unseen_to_nan::Bool=true)

    isempty(files) && error("Cannot read an empty FITS stack.")
    kinds = detect_fits_grid.(files)
    length(unique(kinds)) == 1 ||
        error("Cannot mix HEALPix FITS tables and regular FITS images in one stack.")

    if first(kinds) == :healpix
        expected_ndims === nothing ||
            error("HEALPix FITS stacks do not have a regular image dimensionality.")
        stack = length(files) == 1 ?
            read_healpix_stack(first(files); column=column, T=T, unseen_to_nan=unseen_to_nan) :
            read_healpix_stack(files; column=column, T=T, unseen_to_nan=unseen_to_nan)
        return conversion == 1 ? stack :
            HealpixStack(stack.pixels .* conversion; nside=stack.nside, order=stack.order, coordsys=stack.coordsys)
    end

    if length(files) == 1
        return read_file(first(files), conversion; expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite)
    end

    planes = [read_file(file, conversion; expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite) for file in files]
    first_size = size(first(planes))
    all(size(plane) == first_size for plane in planes) ||
        error("All regular FITS image planes in a stack must have the same shape.")

    return cat(planes...; dims=ndims(first(planes)) + 1)
end

read_fits_grid_stack(file::AbstractString, conversion::Real=1.0; kwargs...) =
    read_fits_grid_stack([file], conversion; kwargs...)

function HealpixStack(maps::AbstractVector)
    isempty(maps) && error("Cannot build a HEALPix stack from an empty map list.")

    first_map = first(maps)
    nside = healpix_nside(first_map)
    order = healpix_order(first_map)
    npix = length(first_map)
    T = promote_type(map(eltype, maps)...)
    pixels = Matrix{T}(undef, npix, length(maps))

    for (idx, map) in enumerate(maps)
        healpix_nside(map) == nside || error("All HEALPix maps must have the same NSIDE.")
        healpix_order(map) == order || error("All HEALPix maps must use the same ordering.")
        length(map) == npix || error("All HEALPix maps must have the same number of pixels.")
        pixels[:, idx] .= collect(map)
    end

    return HealpixStack(pixels, nside, order)
end

function HealpixStack(pixels::AbstractMatrix; nside::Union{Nothing, Integer}=nothing, order::HealpixOrderName=:ring, coordsys=nothing)
    inferred_nside = healpix_nside_from_npix(size(pixels, 1))
    if nside !== nothing && Int(nside) != inferred_nside
        error("Matrix has $(size(pixels, 1)) rows, corresponding to NSIDE=$(inferred_nside), not NSIDE=$(nside).")
    end

    return HealpixStack(_as_matrix(pixels), inferred_nside, normalize_healpix_order(order), _normalize_coordsys(coordsys))
end

function _check_healpix_stack_consistency(q_stack::HealpixStack, u_stack::HealpixStack)
    size(q_stack.pixels) == size(u_stack.pixels) ||
        error("Q and U HEALPix stacks must have the same Npix x Nfreq shape.")
    q_stack.order == u_stack.order ||
        error("Q and U HEALPix stacks must use the same ordering (Q is $(q_stack.order), U is $(u_stack.order)).")
    q_stack.nside == u_stack.nside ||
        error("Q and U HEALPix stacks must have the same NSIDE (Q is $(q_stack.nside), U is $(u_stack.nside)).")
    if q_stack.coordsys !== nothing && u_stack.coordsys !== nothing && q_stack.coordsys != u_stack.coordsys
        error("Q and U HEALPix stacks must use the same coordinate system (Q is $(q_stack.coordsys), U is $(u_stack.coordsys)).")
    end
    return nothing
end

function _healpix_stack_from_input(input; nside=nothing, order::HealpixOrderName=:ring)
    input isa HealpixStack && return input
    input isa AbstractVector && hasproperty(input, :resolution) && return HealpixStack([input])
    input isa AbstractMatrix && return HealpixStack(input; nside=nside, order=order)
    input isa AbstractVector && !isempty(input) && hasproperty(first(input), :resolution) && return HealpixStack(input)

    error("Expected a HealpixStack, a matrix of size Npix x Nfreq, or a vector of Healpix.jl maps.")
end

_is_path_or_path_stack(input) =
    input isa AbstractString ||
    (input isa AbstractVector && !isempty(input) && all(item -> item isa AbstractString, input))

function _read_auto_grid_input(input; conversion::Real=1.0, column=1, T::Type=Float64, expected_ndims=nothing, allow_nonfinite::Bool=false)
    if input isa AbstractString
        return read_fits_grid_stack(input, conversion; column=column, T=T, expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite)
    elseif input isa AbstractVector && !isempty(input) && all(item -> item isa AbstractString, input)
        return read_fits_grid_stack(input, conversion; column=column, T=T, expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite)
    end

    return input
end

_is_healpix_input(input; nside=nothing) =
    input isa HealpixStack ||
    (input isa AbstractVector && hasproperty(input, :resolution)) ||
    (input isa AbstractVector && !isempty(input) && hasproperty(first(input), :resolution)) ||
    (nside !== nothing && input isa AbstractMatrix)

function healpix_maps_from_stack(stack::HealpixStack)
    return [healpix_map(view(stack.pixels, :, i); nside=stack.nside, order=stack.order) for i in axes(stack.pixels, 2)]
end

function healpix_maps_from_stack(pixels::AbstractMatrix; nside::Union{Nothing, Integer}=nothing, order::HealpixOrderName=:ring)
    return healpix_maps_from_stack(HealpixStack(pixels; nside=nside, order=order))
end

function RMSynthesisHealpix(Q, U, nuArray::AbstractArray, PhiArray::AbstractArray;
    nside::Union{Nothing, Integer}=nothing,
    order::HealpixOrderName=:ring,
    log_progress::Bool=false)

    q_stack = _healpix_stack_from_input(Q; nside=nside, order=order)
    u_stack = _healpix_stack_from_input(U; nside=q_stack.nside, order=q_stack.order)

    _check_healpix_stack_consistency(q_stack, u_stack)
    length(nuArray) == size(q_stack.pixels, 2) ||
        error("nuArray length ($(length(nuArray))) must match the number of HEALPix frequency maps ($(size(q_stack.pixels, 2))).")

    fdf, realFDF, imagFDF = RMSynthesis(q_stack.pixels, u_stack.pixels, nuArray, PhiArray; log_progress=log_progress)

    return HealpixRMResult(
        _as_matrix(fdf),
        _as_matrix(realFDF),
        _as_matrix(imagFDF),
        Float64.(collect(PhiArray)),
        q_stack.nside,
        q_stack.order,
        q_stack.coordsys,
    )
end

"""
    RMSynthesisAuto(Q, U, nuArray, PhiArray; kwargs...)

Run RM synthesis while automatically dispatching regular FITS images/cubes to
[`RMSynthesis`](@ref) and HEALPix FITS maps/stacks to
[`RMSynthesisHealpix`](@ref). `Q` and `U` may be arrays, `HealpixStack`s,
vectors of Healpix maps, single FITS image paths, or vectors of HEALPix FITS
paths ordered like `nuArray`.
"""
function RMSynthesisAuto(Q, U, nuArray::AbstractArray, PhiArray::AbstractArray;
    conversion::Real=1.0,
    column=1,
    T::Type=Float64,
    nside::Union{Nothing, Integer}=nothing,
    order::HealpixOrderName=:ring,
    log_progress::Bool=false,
    allow_nonfinite::Bool=false)

    q_input = _read_auto_grid_input(Q; conversion=conversion, column=column, T=T, allow_nonfinite=allow_nonfinite)
    u_input = _read_auto_grid_input(U; conversion=conversion, column=column, T=T, allow_nonfinite=allow_nonfinite)
    q_is_healpix = _is_healpix_input(q_input; nside=nside)
    u_is_healpix = _is_healpix_input(u_input; nside=nside)

    q_is_healpix == u_is_healpix ||
        error("Q and U must both be HEALPix grids or both be regular image/cube grids.")

    if q_is_healpix
        return RMSynthesisHealpix(q_input, u_input, nuArray, PhiArray;
            nside=nside, order=order, log_progress=log_progress)
    end

    return RMSynthesis(q_input, u_input, nuArray, PhiArray; log_progress=log_progress)
end

function _write_coordsys_keyword(filename::AbstractString, coordsys::AbstractString)
    fits = FITS(String(filename), "r+")
    try
        FITSIO.write_key(fits[2], "COORDSYS", String(coordsys), "Pixelisation coordinate system")
    finally
        close(fits)
    end
    return nothing
end

"""
    write_healpix_map(filename, pixels; nside, order, typechar, unit, extname,
                      coordsys, nan_to_unseen=true, overwrite=true)

Write a single HEALPix map to a FITS file. `NaN` pixels are stored as the
UNSEEN sentinel unless `nan_to_unseen=false`. When `coordsys` is given (e.g.
`"G"`, `"E"`, `"C"`), it is stored in the `COORDSYS` header keyword.
"""
function write_healpix_map(filename::AbstractString, pixels::AbstractVector;
    nside::Union{Nothing, Integer}=nothing,
    order::HealpixOrderName=:ring,
    typechar::AbstractString="D",
    unit::AbstractString="",
    extname::AbstractString="MAP",
    coordsys=nothing,
    nan_to_unseen::Bool=true,
    overwrite::Bool=true)

    values = nan_to_unseen ? _restore_unseen(pixels) : pixels
    map = healpix_map(values; nside=nside, order=order)
    fits_path = overwrite ? "!" * String(filename) : String(filename)
    Healpix.saveToFITS(map, fits_path; typechar=typechar, unit=unit, extname=extname)

    coordsys_value = _normalize_coordsys(coordsys)
    coordsys_value === nothing || _write_coordsys_keyword(filename, coordsys_value)
    return filename
end

function _healpix_slice_label(value)
    rounded = round(Float64(value); digits=6)
    prefix = rounded < 0 ? "m" : "p"
    safe_value = replace(string(abs(rounded)), "." => "p")
    return prefix * safe_value
end

"""
    write_healpix_cube(filename, pixels, coordinates; nside, order, unit,
                       extname, coordname, coordunit, coordsys,
                       nan_to_unseen=true, overwrite=true)

Write an `Npix × Nslice` HEALPix stack as a **single** FITS file: one HEALPix
binary-table HDU holding all slices (one vector cell per pixel row) plus a
`COORDS` extension listing the slice coordinates (Faraday depths,
frequencies, ...). The result is readable with [`read_healpix_stack`](@ref)
and [`read_healpix_cube`](@ref).
"""
function write_healpix_cube(filename::AbstractString, pixels::AbstractMatrix, coordinates::AbstractVector;
    nside::Union{Nothing, Integer}=nothing,
    order::HealpixOrderName=:ring,
    unit::AbstractString="",
    extname::AbstractString="CUBE",
    coordname::AbstractString="COORD",
    coordunit::AbstractString="",
    coordsys=nothing,
    nan_to_unseen::Bool=true,
    overwrite::Bool=true)

    size(pixels, 2) == length(coordinates) ||
        error("The number of HEALPix slices ($(size(pixels, 2))) must match the number of coordinates ($(length(coordinates))).")
    overwrite || !isfile(filename) ||
        error("File $(filename) already exists. Pass overwrite=true to replace it.")

    stack = HealpixStack(pixels; nside=nside, order=order, coordsys=coordsys)
    values = nan_to_unseen ? _restore_unseen(stack.pixels) : stack.pixels
    data = permutedims(values)  # Nslice x Npix: one vector cell per pixel row.

    header_keys = String["PIXTYPE", "ORDERING", "NSIDE", "FIRSTPIX", "LASTPIX", "INDXSCHM"]
    header_values = Any["HEALPIX", stack.order == :nested ? "NESTED" : "RING",
                        stack.nside, 0, size(stack.pixels, 1) - 1, "IMPLICIT"]
    header_comments = String["", "", "", "", "", ""]
    if !isempty(unit)
        push!(header_keys, "TUNIT1"); push!(header_values, String(unit)); push!(header_comments, "Pixel unit")
    end
    if stack.coordsys !== nothing
        push!(header_keys, "COORDSYS"); push!(header_values, stack.coordsys); push!(header_comments, "Pixelisation coordinate system")
    end
    header = FITSHeader(header_keys, header_values, header_comments)

    coord_keys = String["COORDNAM"]
    coord_values = Any[String(coordname)]
    coord_comments = String["Physical meaning of the slice coordinate"]
    if !isempty(coordunit)
        push!(coord_keys, "TUNIT1"); push!(coord_values, String(coordunit)); push!(coord_comments, "Coordinate unit")
    end
    coord_header = FITSHeader(coord_keys, coord_values, coord_comments)

    atomic_write_path(String(filename)) do tmp_path
        FITS(tmp_path, "w") do fits
            write(fits, ["PIXELS"], [data]; header=header, name=String(extname))
            write(fits, [String(coordname)], [Float64.(collect(coordinates))]; header=coord_header, name="COORDS")
        end
    end

    return filename
end

"""
    read_healpix_cube(filename; T=Float64, unseen_to_nan=true) -> (stack, coordinates)

Read a single-file HEALPix cube written by [`write_healpix_cube`](@ref).
Returns the [`HealpixStack`](@ref) and the slice coordinate vector (or
`nothing` if the file has no `COORDS` extension).
"""
function read_healpix_cube(filename::AbstractString; T::Type=Float64, unseen_to_nan::Bool=true)
    stack = read_healpix_stack(filename; column=:all, T=T, unseen_to_nan=unseen_to_nan)

    coordinates = FITS(filename) do fits
        for idx in 1:length(fits)
            hdu = fits[idx]
            header = read_header(hdu)
            hdu_extname = uppercase(strip(String(_fits_header_value(header, "EXTNAME", ""))))
            if hdu_extname == "COORDS"
                names = FITSIO.colnames(hdu)
                isempty(names) && return nothing
                return Float64.(read(hdu, first(names)))
            end
        end
        return nothing
    end

    return stack, coordinates
end

"""
    write_healpix_stack(output_dir, pixels, basename, coordinates; format=:files, kwargs...)

Write an `Npix × Nslice` HEALPix stack. With `format=:files` (default), one
FITS file per slice is written, labelled by its coordinate. With
`format=:cube`, a single `basename.fits` file containing all slices plus a
`COORDS` extension is written instead (see [`write_healpix_cube`](@ref)).
Returns the list of written paths.
"""
function write_healpix_stack(output_dir::AbstractString, pixels::AbstractMatrix, basename::AbstractString, coordinates::AbstractVector;
    nside::Union{Nothing, Integer}=nothing,
    order::HealpixOrderName=:ring,
    unit::AbstractString="",
    extname::AbstractString="MAP",
    coordsys=nothing,
    format::Symbol=:files,
    nan_to_unseen::Bool=true,
    overwrite::Bool=true)

    size(pixels, 2) == length(coordinates) ||
        error("The number of HEALPix slices ($(size(pixels, 2))) must match the number of coordinates ($(length(coordinates))).")

    mkpath(output_dir)
    stack = HealpixStack(pixels; nside=nside, order=order, coordsys=coordsys)

    if format == :cube
        path = joinpath(output_dir, "$(basename).fits")
        write_healpix_cube(path, stack.pixels, collect(coordinates);
            nside=stack.nside, order=stack.order, unit=unit, extname=extname,
            coordsys=stack.coordsys, nan_to_unseen=nan_to_unseen, overwrite=overwrite)
        return String[path]
    end
    format == :files || error("Unknown HEALPix stack format $(repr(format)). Use :files or :cube.")

    paths = String[]
    for idx in axes(stack.pixels, 2)
        label = _healpix_slice_label(coordinates[idx])
        path = joinpath(output_dir, "$(basename)_$(lpad(idx, 4, "0"))_$(label).fits")
        write_healpix_map(path, view(stack.pixels, :, idx); nside=stack.nside, order=stack.order,
            unit=unit, extname=extname, coordsys=stack.coordsys,
            nan_to_unseen=nan_to_unseen, overwrite=overwrite)
        push!(paths, path)
    end

    return paths
end

function write_healpix_rm_result(output_dir::AbstractString, result::HealpixRMResult;
    prefix::AbstractString="",
    fdf_unit::AbstractString="",
    format::Symbol=:files,
    overwrite::Bool=true)

    mkpath(output_dir)
    stem(name) = isempty(prefix) ? name : "$(prefix)_$(name)"

    fdf_paths = write_healpix_stack(output_dir, result.fdf, stem("FDF"), result.phi;
        nside=result.nside, order=result.order, coordsys=result.coordsys,
        unit=fdf_unit, extname="FDF", format=format, overwrite=overwrite)
    real_paths = write_healpix_stack(output_dir, result.realFDF, stem("realFDF"), result.phi;
        nside=result.nside, order=result.order, coordsys=result.coordsys,
        unit=fdf_unit, extname="REALFDF", format=format, overwrite=overwrite)
    imag_paths = write_healpix_stack(output_dir, result.imagFDF, stem("imagFDF"), result.phi;
        nside=result.nside, order=result.order, coordsys=result.coordsys,
        unit=fdf_unit, extname="IMAGFDF", format=format, overwrite=overwrite)

    return (; fdf=fdf_paths, realFDF=real_paths, imagFDF=imag_paths)
end

# ---------------------------------------------------------------------------
# Resolution / ordering / smoothing helpers
# ---------------------------------------------------------------------------

function _require_pow2_nside(nside::Integer, context::AbstractString)
    ispow2(Int(nside)) || error("$(context) requires a power-of-two NSIDE, got $(nside).")
    return Int(nside)
end

# Permutation such that out[p] = in[perm[p]] converts `from` ordering into
# `to` ordering (1-based HEALPix pixel indices).
function _healpix_order_permutation(nside::Int, from::Symbol, to::Symbol)
    res = Healpix.Resolution(nside)
    npix = Healpix.nside2npix(nside)
    perm = Vector{Int}(undef, npix)
    if from == :ring && to == :nested
        @inbounds for p in 1:npix
            perm[p] = Healpix.nest2ring(res, p)
        end
    elseif from == :nested && to == :ring
        @inbounds for p in 1:npix
            perm[p] = Healpix.ring2nest(res, p)
        end
    else
        error("Unsupported HEALPix reordering $(from) -> $(to).")
    end
    return perm
end

"""
    healpix_reorder(stack, order) -> HealpixStack

Convert a [`HealpixStack`](@ref) between `:ring` and `:nested` pixel ordering.
Returns the input unchanged when it already uses the requested ordering.
"""
function healpix_reorder(stack::HealpixStack, order::HealpixOrderName)
    target = normalize_healpix_order(order)
    stack.order == target && return stack
    _require_pow2_nside(stack.nside, "HEALPix reordering")

    perm = _healpix_order_permutation(stack.nside, stack.order, target)
    return HealpixStack(stack.pixels[perm, :], stack.nside, target, stack.coordsys)
end

# NaN-aware up/degrade of a single map vector using NESTED hierarchy
# arithmetic (children of nested pixel p are p*r .. p*r + r - 1, 0-based).
function _healpix_udgrade_vector(values::AbstractVector{T}, nside_in::Int, nside_out::Int, order::Symbol) where {T}
    nside_out == nside_in && return collect(float(T), values)

    res_in = Healpix.Resolution(nside_in)
    res_out = Healpix.Resolution(nside_out)
    npix_out = Healpix.nside2npix(nside_out)
    F = float(T)
    out = Vector{F}(undef, npix_out)

    if nside_out > nside_in
        # Upgrade: each output pixel inherits its parent's value.
        ratio = (nside_out ÷ nside_in)^2
        @inbounds for p_out in 1:npix_out
            nest_out = order == :ring ? Healpix.ring2nest(res_out, p_out) - 1 : p_out - 1
            nest_in = nest_out ÷ ratio
            p_in = order == :ring ? Healpix.nest2ring(res_in, nest_in + 1) : nest_in + 1
            out[p_out] = values[p_in]
        end
    else
        # Degrade: mean of the children, ignoring NaN (masked) pixels.
        ratio = (nside_in ÷ nside_out)^2
        @inbounds for p_out in 1:npix_out
            nest_out = order == :ring ? Healpix.ring2nest(res_out, p_out) - 1 : p_out - 1
            acc = zero(F)
            cnt = 0
            for child in (nest_out * ratio):(nest_out * ratio + ratio - 1)
                p_in = order == :ring ? Healpix.nest2ring(res_in, child + 1) : child + 1
                v = F(values[p_in])
                if !isnan(v)
                    acc += v
                    cnt += 1
                end
            end
            out[p_out] = cnt == 0 ? F(NaN) : acc / cnt
        end
    end

    return out
end

"""
    healpix_udgrade(stack, nside_out) -> HealpixStack
    healpix_udgrade(pixels::AbstractVector, nside_out; order=:ring) -> Vector

Up- or degrade a HEALPix stack (or a single map vector) to `nside_out`.
Upgrading replicates each parent pixel into its children; degrading averages
the children, ignoring `NaN` (masked) pixels — an output pixel whose children
are all masked stays `NaN`. Both NSIDE values must be powers of two.
"""
function healpix_udgrade(stack::HealpixStack, nside_out::Integer)
    nside_out = Int(nside_out)
    nside_out == stack.nside && return stack
    _require_pow2_nside(stack.nside, "HEALPix up/degrading")
    _require_pow2_nside(nside_out, "HEALPix up/degrading")

    npix_out = Healpix.nside2npix(nside_out)
    F = float(eltype(stack.pixels))
    out = Matrix{F}(undef, npix_out, size(stack.pixels, 2))
    for j in axes(stack.pixels, 2)
        out[:, j] = _healpix_udgrade_vector(view(stack.pixels, :, j), stack.nside, nside_out, stack.order)
    end

    return HealpixStack(out, nside_out, stack.order, stack.coordsys)
end

function healpix_udgrade(pixels::AbstractVector, nside_out::Integer; order::HealpixOrderName=:ring)
    nside_in = healpix_nside_from_npix(length(pixels))
    nside_out = Int(nside_out)
    nside_out == nside_in && return collect(float(eltype(pixels)), pixels)
    _require_pow2_nside(nside_in, "HEALPix up/degrading")
    _require_pow2_nside(nside_out, "HEALPix up/degrading")

    return _healpix_udgrade_vector(pixels, nside_in, nside_out, normalize_healpix_order(order))
end

# Conform a stack to a target NSIDE/ordering, logging what is being changed.
function _conform_healpix_stack(stack::HealpixStack, nside::Integer, order::Symbol; label::AbstractString="HEALPix input")
    if stack.order != order
        @warn "Reordering $(label) from $(stack.order) to $(order) to match the reference field."
        stack = healpix_reorder(stack, order)
    end
    if stack.nside != Int(nside)
        @warn "Resampling $(label) from NSIDE=$(stack.nside) to NSIDE=$(nside) to match the reference field."
        stack = healpix_udgrade(stack, nside)
    end
    return stack
end

function _gaussian_fwhm_rad(; fwhm_rad=nothing, fwhm_deg=nothing, fwhm_arcmin=nothing)
    given = count(!isnothing, (fwhm_rad, fwhm_deg, fwhm_arcmin))
    given == 1 || error("Specify exactly one of `fwhm_rad`, `fwhm_deg`, or `fwhm_arcmin`.")
    value = fwhm_rad !== nothing ? Float64(fwhm_rad) :
            fwhm_deg !== nothing ? deg2rad(Float64(fwhm_deg)) :
            deg2rad(Float64(fwhm_arcmin) / 60.0)
    value > 0 || error("The smoothing FWHM must be positive, got $(value) rad.")
    return value
end

# Smooth one ring-ordered Float64 map vector: map -> alm, multiply by the
# Gaussian beam window, alm -> map.
function _healpix_smooth_ring_vector(values::Vector{Float64}, nside::Int, fwhm_rad::Float64, lmax::Int)
    map_in = Healpix.HealpixMap{Float64, Healpix.RingOrder}(values)
    alm = Healpix.map2alm(map_in; lmax=lmax, mmax=lmax)
    beam = Healpix.gaussbeam(fwhm_rad, lmax)
    @inbounds for l in 0:lmax, m in 0:l
        alm.alm[Healpix.almIndex(alm, l, m)] *= beam[l + 1]
    end
    return collect(Healpix.alm2map(alm, nside))
end

"""
    healpix_smooth(stack; fwhm_arcmin | fwhm_deg | fwhm_rad, lmax=3*nside-1) -> HealpixStack
    healpix_smooth(pixels::AbstractVector; kwargs..., order=:ring, nside=nothing) -> Vector

Convolve HEALPix maps with a Gaussian beam in spherical-harmonic space
(`map2alm` → multiply by `gaussbeam` → `alm2map`). Give the beam width with
exactly one of `fwhm_arcmin`, `fwhm_deg`, or `fwhm_rad`.

Masked (`NaN`) pixels are handled by smoothing the map with masked pixels set
to zero, smoothing the binary validity mask the same way, and dividing the two
(masked smoothing). Pixels that were masked on input remain `NaN` on output.

Nested-ordered inputs are converted to ring for the transform and converted
back, so the output ordering always matches the input.
"""
function healpix_smooth(stack::HealpixStack;
    fwhm_rad=nothing, fwhm_deg=nothing, fwhm_arcmin=nothing,
    lmax::Union{Nothing, Integer}=nothing)

    fwhm = _gaussian_fwhm_rad(; fwhm_rad=fwhm_rad, fwhm_deg=fwhm_deg, fwhm_arcmin=fwhm_arcmin)
    lmax_value = lmax === nothing ? 3 * stack.nside - 1 : Int(lmax)
    lmax_value >= 0 || error("`lmax` must be non-negative, got $(lmax_value).")

    original_order = stack.order
    ring_stack = healpix_reorder(stack, :ring)

    npix = size(ring_stack.pixels, 1)
    out = Matrix{Float64}(undef, npix, size(ring_stack.pixels, 2))

    for j in axes(ring_stack.pixels, 2)
        column = Float64.(view(ring_stack.pixels, :, j))
        invalid = isnan.(column)
        filled = ifelse.(invalid, 0.0, column)
        weights = Float64.(.!invalid)
        smoothed = _healpix_smooth_ring_vector(filled, ring_stack.nside, fwhm, lmax_value)
        weight_sm = _healpix_smooth_ring_vector(weights, ring_stack.nside, fwhm, lmax_value)
        @inbounds for p in 1:npix
            out[p, j] = (invalid[p] || weight_sm[p] <= 1e-6) ? NaN : smoothed[p] / weight_sm[p]
        end
    end

    result = HealpixStack(out, ring_stack.nside, :ring, ring_stack.coordsys)
    return original_order == :ring ? result : healpix_reorder(result, original_order)
end

function healpix_smooth(pixels::AbstractVector;
    fwhm_rad=nothing, fwhm_deg=nothing, fwhm_arcmin=nothing,
    lmax::Union{Nothing, Integer}=nothing,
    order::HealpixOrderName=:ring,
    nside::Union{Nothing, Integer}=nothing)

    stack = HealpixStack(reshape(collect(pixels), :, 1); nside=nside, order=order)
    smoothed = healpix_smooth(stack; fwhm_rad=fwhm_rad, fwhm_deg=fwhm_deg, fwhm_arcmin=fwhm_arcmin, lmax=lmax)
    return vec(smoothed.pixels)
end
