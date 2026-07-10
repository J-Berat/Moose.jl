"""
Helpers for reading FITS files and aligning data cubes with the chosen line of sight.
"""

function _validate_fits_array(data, file; expected_ndims=nothing, allow_nonfinite::Bool=false)
    data isa AbstractArray || error("No image data found in $(file).")
    isempty(data) && error("The first readable image HDU in $(file) is empty.")
    eltype(data) <: Number || error("FITS image $(file) must contain numeric data, got $(eltype(data)).")

    if expected_ndims !== nothing && ndims(data) != expected_ndims
        error("FITS image $(file) has $(ndims(data)) dimensions; expected $(expected_ndims)D data.")
    end

    if !allow_nonfinite
        bad_index = findfirst(x -> !isfinite(x), data)
        bad_index === nothing || error("FITS image $(file) contains a non-finite value at index $(bad_index).")
    end

    return data
end

function read_FITS_file(file; expected_ndims=nothing, allow_nonfinite::Bool=false)
    last_error = nothing

    data = FITS(file) do fits
        for hdu in fits
            candidate = try
                read(hdu)
            catch err
                last_error = err
                continue
            end

            candidate isa AbstractArray || continue
            isempty(candidate) && continue
            return _validate_fits_array(candidate, file; expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite)
        end

        return nothing
    end

    if data === nothing
        detail = last_error === nothing ? "" : " Last read error: $(last_error)."
        error("No readable non-empty image HDU found in $(file)." * detail)
    end

    return data
end

function read_file(file::AbstractString, conversion; expected_ndims=nothing, allow_nonfinite::Bool=false)
    if isdefined(@__MODULE__, :is_hdf5_path) && is_hdf5_path(file)
        return read_HDF5_file(file; expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite) .* conversion
    end

    return read_FITS_file(file; expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite) .* conversion
end

function permute_dims(array, LOS)
    if LOS == "x"
        permutedims(array, [2, 3, 1])
    elseif LOS == "y"
        permutedims(array, [3, 1, 2])
    else
        array
    end
end

read_optional_file(file, conversion, LOS) = isfile(file) ? permute_dims(read_file(file, conversion; expected_ndims=3), LOS) : nothing
