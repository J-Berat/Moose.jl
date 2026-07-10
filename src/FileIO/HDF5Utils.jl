"""
Helpers for reading regular simulation cubes from HDF5 files.
"""

const HDF5_EXTS = (".h5", ".hdf5")

struct HDF5DatasetSource
    file::String
    dataset::String
end

is_hdf5_path(path::AbstractString) = lowercase(splitext(path)[2]) in HDF5_EXTS
source_path(source::AbstractString) = String(source)
source_path(source::HDF5DatasetSource) = source.file
source_label(source::AbstractString) = String(source)
source_label(source::HDF5DatasetSource) = "$(source.file):$(source.dataset)"

function _hdf5_dataset_paths(group, prefix::AbstractString="")
    paths = String[]
    for key in keys(group)
        name = String(key)
        path = isempty(prefix) ? name : string(prefix, "/", name)
        child = group[name]
        if child isa HDF5.Dataset
            push!(paths, path)
        elseif child isa HDF5.Group || child isa HDF5.File
            append!(paths, _hdf5_dataset_paths(child, path))
        end
    end
    return paths
end

function hdf5_dataset_paths(file::AbstractString)
    validation_error = ensure_readable_file(file; expected_exts=collect(HDF5_EXTS))
    validation_error === nothing || error(validation_error)

    return h5open(file, "r") do h5
        _hdf5_dataset_paths(h5)
    end
end

function _matching_hdf5_dataset(paths::Vector{String}, dataset::AbstractString)
    requested = strip(String(dataset), '/')
    exact = findfirst(==(requested), paths)
    exact !== nothing && return paths[exact]

    matches = [path for path in paths if basename(path) == requested]
    if length(matches) == 1
        return only(matches)
    elseif length(matches) > 1
        error("Multiple HDF5 datasets named $(dataset): $(join(matches, ", ")). Use unique dataset names.")
    end

    return nothing
end

function find_hdf5_dataset(file::AbstractString, dataset::AbstractString)
    is_hdf5_path(file) || return nothing
    paths = hdf5_dataset_paths(file)
    return _matching_hdf5_dataset(paths, dataset)
end

function _resolve_hdf5_dataset(file::AbstractString, dataset)
    paths = hdf5_dataset_paths(file)
    isempty(paths) && error("No datasets found in HDF5 file $(file).")

    if dataset !== nothing
        match = _matching_hdf5_dataset(paths, String(dataset))
        match !== nothing && return match
        error("Dataset $(dataset) was not found in HDF5 file $(file). Available datasets: $(join(paths, ", ")).")
    end

    length(paths) == 1 && return only(paths)

    inferred = _matching_hdf5_dataset(paths, splitext(basename(file))[1])
    inferred !== nothing && return inferred

    error("HDF5 file $(file) contains multiple datasets ($(join(paths, ", "))). Name one dataset $(splitext(basename(file))[1]) or use a shared simulation file with datasets named Bx, By, Bz, density, and temperature.")
end

function _validate_hdf5_array(data, file; dataset=nothing, expected_ndims=nothing, allow_nonfinite::Bool=false)
    label = dataset === nothing ? file : "$(file):$(dataset)"
    data isa AbstractArray || error("HDF5 dataset $(label) must contain array data.")
    isempty(data) && error("HDF5 dataset $(label) is empty.")
    eltype(data) <: Number || error("HDF5 dataset $(label) must contain numeric data, got $(eltype(data)).")

    if expected_ndims !== nothing && ndims(data) != expected_ndims
        error("HDF5 dataset $(label) has $(ndims(data)) dimensions; expected $(expected_ndims)D data.")
    end

    if !allow_nonfinite
        bad_index = findfirst(x -> !isfinite(x), data)
        bad_index === nothing || error("HDF5 dataset $(label) contains a non-finite value at index $(bad_index).")
    end

    return data
end

function read_HDF5_file(file::AbstractString; dataset=nothing, expected_ndims=nothing, allow_nonfinite::Bool=false)
    resolved_dataset = _resolve_hdf5_dataset(file, dataset)
    data = h5open(file, "r") do h5
        read(h5[resolved_dataset])
    end
    return _validate_hdf5_array(data, file; dataset=resolved_dataset, expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite)
end

function read_file(source::HDF5DatasetSource, conversion; expected_ndims=nothing, allow_nonfinite::Bool=false)
    return read_HDF5_file(source.file; dataset=source.dataset, expected_ndims=expected_ndims, allow_nonfinite=allow_nonfinite) .* conversion
end
