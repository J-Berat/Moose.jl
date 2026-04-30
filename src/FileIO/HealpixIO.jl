"""
Helpers for using MOOSE algorithms with HEALPix maps.

The core MOOSE RM-synthesis routine works on arrays whose last dimension is
frequency.  These helpers bridge Healpix.jl maps to that representation while
preserving NSIDE and ordering metadata for writing results back to HEALPix FITS.
"""

const HealpixOrderName = Union{Symbol, AbstractString}

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

function read_healpix_stack(files::AbstractVector{<:AbstractString}; column=1, T::Type=Float64)
    isempty(files) && error("Cannot read an empty HEALPix stack.")
    return HealpixStack([read_healpix_map(file; column=column, T=T) for file in files])
end

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
    input isa AbstractMatrix && return HealpixStack(input; nside=nside, order=order)
    input isa AbstractVector && !isempty(input) && hasproperty(first(input), :resolution) && return HealpixStack(input)

    error("Expected a HealpixStack, a matrix of size Npix x Nfreq, or a vector of Healpix.jl maps.")
end

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
