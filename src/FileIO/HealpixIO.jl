"""
Helpers for using MOOSE algorithms with HEALPix maps.

The core MOOSE RM-synthesis routine works on arrays whose last dimension is
frequency.  These helpers bridge Healpix.jl maps to that representation while
preserving NSIDE and ordering metadata for writing results back to HEALPix FITS.
"""

const HealpixOrderName = Union{Symbol, AbstractString}
const FITSGridKind = Symbol

struct HealpixStack{T}
    pixels::Matrix{T}
    nside::Int
    order::Symbol
end

struct HealpixRMResult{T}
    fdf::Matrix{T}
    realFDF::Matrix{T}
    imagFDF::Matrix{T}
    phi::Vector{Float64}
    nside::Int
    order::Symbol
end

function _fits_header_value(header, key::AbstractString, default=nothing)
    try
        return header[key]
    catch
        return default
    end
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
            pixtype = _fits_header_value(header, "PIXTYPE", "")
            ordering = _fits_header_value(header, "ORDERING", nothing)
            nside = _fits_header_value(header, "NSIDE", nothing)
            xtension = uppercase(String(_fits_header_value(header, "XTENSION", "")))
            naxis = _fits_header_value(header, "NAXIS", 0)

            if uppercase(String(pixtype)) == "HEALPIX" || (ordering !== nothing && nside !== nothing && xtension == "BINTABLE")
                return :healpix
            end

            parsed_naxis = tryparse(Int, string(naxis))
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
    normalized = lowercase(String(order))
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

function healpix_order(map)
    type_name = string(typeof(map))
    occursin("NestedOrder", type_name) && return :nested
    occursin("RingOrder", type_name) && return :ring
    error("Could not infer HEALPix ordering from map type $(typeof(map)).")
end

healpix_nside(map) = Int(map.resolution.nside)

function read_healpix_map(filename::AbstractString; column=1, T::Type=Float64)
    validation_error = ensure_readable_file(filename; expected_exts=[".fits", ".fit", ".fts"])
    validation_error === nothing || error(validation_error)

    return Healpix.readMapFromFITS(String(filename), column, T)
end

function _healpix_table_hdu_info(filename::AbstractString)
    validation_error = ensure_readable_file(filename; expected_exts=[".fits", ".fit", ".fts"])
    validation_error === nothing || error(validation_error)

    info = FITS(filename) do fits
        for idx in 1:length(fits)
            hdu = fits[idx]
            header = read_header(hdu)
            pixtype = uppercase(String(_fits_header_value(header, "PIXTYPE", "")))
            ordering = _fits_header_value(header, "ORDERING", nothing)
            nside = _fits_header_value(header, "NSIDE", nothing)
            xtension = uppercase(String(_fits_header_value(header, "XTENSION", "")))
            if pixtype == "HEALPIX" || (ordering !== nothing && nside !== nothing && xtension == "BINTABLE")
                return (;
                    hdu_index = idx,
                    nside = Int(nside),
                    order = normalize_healpix_order(ordering),
                    colnames = FITSIO.colnames(hdu),
                    nrows = Int(_fits_header_value(header, "NAXIS2", 0)),
                )
            end
        end

        return nothing
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

function _read_healpix_column_data(filename::AbstractString, column, info, ::Type{T}) where {T}
    FITS(filename) do fits
        hdu = fits[info.hdu_index]
        name = column isa Integer ? info.colnames[Int(column)] : String(column)
        data = read(hdu, name; case_sensitive=false)
        return T.(data)
    end
end

function _healpix_matrix_from_column_data(data::AbstractArray, info)
    npix = Healpix.nside2npix(info.nside)
    if ndims(data) == 1
        length(data) == npix || error("HEALPix column has $(length(data)) values, expected $(npix) for NSIDE=$(info.nside).")
        return reshape(collect(data), npix, 1)
    elseif ndims(data) == 2
        if size(data, 1) == npix
            return Matrix(data)
        elseif size(data, 2) == npix
            return permutedims(data)
        elseif length(data) % npix == 0
            return reshape(vec(data), npix, length(data) ÷ npix)
        end
    end

    error("Cannot interpret HEALPix table column with size $(size(data)) as an NSIDE=$(info.nside) map or cube.")
end

function read_healpix_stack(filename::AbstractString; column=:all, T::Type=Float64)
    info = _healpix_table_hdu_info(filename)
    columns = _normalize_healpix_columns(column)
    columns = columns === :all ? collect(eachindex(info.colnames)) : columns
    isempty(columns) && error("Cannot read an empty HEALPix column selection.")

    chunks = Matrix{T}[]
    for col in columns
        data = _read_healpix_column_data(filename, col, info, T)
        push!(chunks, _healpix_matrix_from_column_data(data, info))
    end

    pixels = reduce(hcat, chunks)
    return HealpixStack(pixels; nside=info.nside, order=info.order)
end

function read_healpix_stack(files::AbstractVector{<:AbstractString}; column=1, T::Type=Float64)
    isempty(files) && error("Cannot read an empty HEALPix stack.")
    return HealpixStack([read_healpix_map(file; column=column, T=T) for file in files])
end

function read_fits_grid(filename::AbstractString, conversion::Real=1.0;
    column=1,
    T::Type=Float64,
    expected_ndims=nothing,
    allow_nonfinite::Bool=false)

    kind = detect_fits_grid(filename)
    if kind == :healpix
        expected_ndims === nothing ||
            error("$(filename) is a HEALPix FITS table, not a $(expected_ndims)D image cube.")
        stack = read_healpix_stack(filename; column=column, T=T)
        if size(stack.pixels, 2) == 1
            map = healpix_map(view(stack.pixels, :, 1); nside=stack.nside, order=stack.order)
            return conversion == 1 ? map : healpix_map(collect(map) .* conversion; nside=stack.nside, order=stack.order)
        end

        return conversion == 1 ? stack : HealpixStack(stack.pixels .* conversion; nside=stack.nside, order=stack.order)
    end

    return read_file(filename, conversion; expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite)
end

function read_fits_grid_stack(files::AbstractVector{<:AbstractString}, conversion::Real=1.0;
    column=1,
    T::Type=Float64,
    expected_ndims=nothing,
    allow_nonfinite::Bool=false)

    isempty(files) && error("Cannot read an empty FITS stack.")
    kinds = detect_fits_grid.(files)
    length(unique(kinds)) == 1 ||
        error("Cannot mix HEALPix FITS tables and regular FITS images in one stack.")

    if first(kinds) == :healpix
        expected_ndims === nothing ||
            error("HEALPix FITS stacks do not have a regular image dimensionality.")
        stack = length(files) == 1 ? read_healpix_stack(first(files); column=column, T=T) : read_healpix_stack(files; column=column, T=T)
        return conversion == 1 ? stack : HealpixStack(stack.pixels .* conversion; nside=stack.nside, order=stack.order)
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

function HealpixStack(pixels::AbstractMatrix; nside::Union{Nothing, Integer}=nothing, order::HealpixOrderName=:ring)
    inferred_nside = healpix_nside_from_npix(size(pixels, 1))
    if nside !== nothing && Int(nside) != inferred_nside
        error("Matrix has $(size(pixels, 1)) rows, corresponding to NSIDE=$(inferred_nside), not NSIDE=$(nside).")
    end

    return HealpixStack(Matrix(pixels), inferred_nside, normalize_healpix_order(order))
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

    size(q_stack.pixels) == size(u_stack.pixels) ||
        error("Q and U HEALPix stacks must have the same Npix x Nfreq shape.")
    length(nuArray) == size(q_stack.pixels, 2) ||
        error("nuArray length ($(length(nuArray))) must match the number of HEALPix frequency maps ($(size(q_stack.pixels, 2))).")

    fdf, realFDF, imagFDF = RMSynthesis(q_stack.pixels, u_stack.pixels, nuArray, PhiArray; log_progress=log_progress)

    return HealpixRMResult(
        Matrix(fdf),
        Matrix(realFDF),
        Matrix(imagFDF),
        Float64.(collect(PhiArray)),
        q_stack.nside,
        q_stack.order,
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

function write_healpix_map(filename::AbstractString, pixels::AbstractVector;
    nside::Union{Nothing, Integer}=nothing,
    order::HealpixOrderName=:ring,
    typechar::AbstractString="D",
    unit::AbstractString="",
    extname::AbstractString="MAP",
    overwrite::Bool=true)

    map = healpix_map(pixels; nside=nside, order=order)
    fits_path = overwrite ? "!" * String(filename) : String(filename)
    Healpix.saveToFITS(map, fits_path; typechar=typechar, unit=unit, extname=extname)
    return filename
end

function _healpix_slice_label(value)
    rounded = round(Float64(value); digits=6)
    prefix = rounded < 0 ? "m" : "p"
    safe_value = replace(string(abs(rounded)), "." => "p")
    return prefix * safe_value
end

function write_healpix_stack(output_dir::AbstractString, pixels::AbstractMatrix, basename::AbstractString, coordinates::AbstractVector;
    nside::Union{Nothing, Integer}=nothing,
    order::HealpixOrderName=:ring,
    unit::AbstractString="",
    extname::AbstractString="MAP",
    overwrite::Bool=true)

    size(pixels, 2) == length(coordinates) ||
        error("The number of HEALPix slices ($(size(pixels, 2))) must match the number of coordinates ($(length(coordinates))).")

    mkpath(output_dir)
    stack = HealpixStack(pixels; nside=nside, order=order)
    paths = String[]

    for idx in axes(stack.pixels, 2)
        label = _healpix_slice_label(coordinates[idx])
        path = joinpath(output_dir, "$(basename)_$(lpad(idx, 4, "0"))_$(label).fits")
        write_healpix_map(path, view(stack.pixels, :, idx); nside=stack.nside, order=stack.order,
            unit=unit, extname=extname, overwrite=overwrite)
        push!(paths, path)
    end

    return paths
end

function write_healpix_rm_result(output_dir::AbstractString, result::HealpixRMResult;
    prefix::AbstractString="",
    fdf_unit::AbstractString="",
    overwrite::Bool=true)

    mkpath(output_dir)
    stem(name) = isempty(prefix) ? name : "$(prefix)_$(name)"

    fdf_paths = write_healpix_stack(output_dir, result.fdf, stem("FDF"), result.phi;
        nside=result.nside, order=result.order, unit=fdf_unit, extname="FDF", overwrite=overwrite)
    real_paths = write_healpix_stack(output_dir, result.realFDF, stem("realFDF"), result.phi;
        nside=result.nside, order=result.order, unit=fdf_unit, extname="REALFDF", overwrite=overwrite)
    imag_paths = write_healpix_stack(output_dir, result.imagFDF, stem("imagFDF"), result.phi;
        nside=result.nside, order=result.order, unit=fdf_unit, extname="IMAGFDF", overwrite=overwrite)

    return (; fdf=fdf_paths, realFDF=real_paths, imagFDF=imag_paths)
end
