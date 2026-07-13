"""
Utilities for rasterizing leaf-cell adaptive meshes onto the regular cartesian
grid consumed by the synthetic-observation pipeline.

The AMR input is deliberately format-neutral: geometry and fields live in
HDF5 datasets, while `field_sources.amr` describes their names.  This avoids
tying MOOSE to one simulation code's native binary format.
"""

struct AMRGeometry
    centers::Matrix{Float64} # Ncell × 3
    widths::Matrix{Float64}  # Ncell × 3
    bounds::NTuple{3, Tuple{Float64, Float64}}
    shape::NTuple{3, Int}
    strict::Bool
    tolerance::Float64
end

function _amr_get(config::AbstractDict, key::AbstractString, default=nothing)
    for candidate in (key, Symbol(key))
        try
            haskey(config, candidate) && return config[candidate]
        catch
        end
    end
    return default
end

function amr_config(field_sources)
    field_sources isa AbstractDict || return nothing
    value = _amr_get(field_sources, "amr", nothing)
    value === nothing && return nothing
    value isa AbstractDict || throw_config_error(
        "`field_sources.amr` must be an object describing the AMR cell geometry.";
        code=:invalid_field_sources)
    return value
end

function _amr_shape(config::AbstractDict)
    raw = _amr_get(config, "shape", nothing)
    raw === nothing && throw_config_error(
        "`field_sources.amr.shape` is required and must contain [nx, ny, nz].";
        code=:invalid_field_sources)
    values = if raw isa AbstractVector
        length(raw) == 3 || throw_config_error("`field_sources.amr.shape` must have three elements."; code=:invalid_field_sources)
        Tuple(Int.(raw))
    elseif raw isa AbstractDict
        (Int(_amr_get(raw, "x", 0)), Int(_amr_get(raw, "y", 0)), Int(_amr_get(raw, "z", 0)))
    else
        throw_config_error("`field_sources.amr.shape` must be [nx, ny, nz] or an object with x/y/z."; code=:invalid_field_sources)
    end
    all(>(0), values) || throw_config_error("Every AMR output dimension must be positive; got $(values)."; code=:invalid_field_sources)
    return values
end

function _amr_bounds(config::AbstractDict)
    raw = _amr_get(config, "bounds", [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
    axes = if raw isa AbstractVector && length(raw) == 6
        ((raw[1], raw[2]), (raw[3], raw[4]), (raw[5], raw[6]))
    elseif raw isa AbstractVector && length(raw) == 3
        Tuple(begin
            pair isa AbstractVector && length(pair) == 2 || throw_config_error(
                "Each entry in `field_sources.amr.bounds` must be [min, max]."; code=:invalid_field_sources)
            (pair[1], pair[2])
        end for pair in raw)
    elseif raw isa AbstractDict
        Tuple(begin
            pair = _amr_get(raw, axis, nothing)
            pair isa AbstractVector && length(pair) == 2 || throw_config_error(
                "`field_sources.amr.bounds.$(axis)` must be [min, max]."; code=:invalid_field_sources)
            (pair[1], pair[2])
        end for axis in ("x", "y", "z"))
    else
        throw_config_error("`field_sources.amr.bounds` must describe x/y/z min/max pairs."; code=:invalid_field_sources)
    end
    result = ntuple(i -> (Float64(axes[i][1]), Float64(axes[i][2])), 3)
    all(isfinite(pair[1]) && isfinite(pair[2]) && pair[2] > pair[1] for pair in result) ||
        throw_config_error("AMR bounds must be finite with max > min."; code=:invalid_field_sources)
    return result
end

function _amr_geometry_file(simu::AbstractString, config::AbstractDict, fallback_file::AbstractString)
    raw = _amr_get(config, "file", fallback_file)
    path = expanduser(String(raw))
    return isabspath(path) ? path : joinpath(simu, path)
end

function _read_amr_vector(file::AbstractString, dataset, label::AbstractString)
    dataset === nothing && throw_config_error("`field_sources.amr.$(label)` must name an HDF5 dataset."; code=:invalid_field_sources)
    data = vec(read_HDF5_file(file; dataset=String(dataset), expected_ndims=nothing))
    all(isfinite, data) || throw_config_error("AMR dataset $(dataset) contains non-finite values."; code=:invalid_field_sources)
    return Float64.(data)
end

function _amr_widths(file::AbstractString, config::AbstractDict, ncell::Int, bounds)
    size_dataset = _amr_get(config, "size", nothing)
    level_dataset = _amr_get(config, "level", nothing)
    (size_dataset === nothing) == (level_dataset === nothing) && throw_config_error(
        "`field_sources.amr` must define exactly one of `size` or `level`.";
        code=:invalid_field_sources)

    if level_dataset !== nothing
        levels = _read_amr_vector(file, level_dataset, "level")
        length(levels) == ncell || throw_config_error("AMR level and coordinate datasets have different lengths."; code=:invalid_field_sources)
        offset = Int(_amr_get(config, "level_offset", 0))
        exponents = levels .- offset
        all(level -> isinteger(level) && level >= 0, exponents) || throw_config_error(
            "AMR levels minus `level_offset` must be non-negative integers."; code=:invalid_field_sources)
        widths = Matrix{Float64}(undef, ncell, 3)
        for axis in 1:3
            extent = bounds[axis][2] - bounds[axis][1]
            widths[:, axis] .= extent ./ exp2.(exponents)
        end
        return widths
    end

    raw = read_HDF5_file(file; dataset=String(size_dataset), expected_ndims=nothing)
    values = Float64.(raw)
    if length(values) == ncell
        return repeat(reshape(vec(values), ncell, 1), 1, 3)
    elseif ndims(values) == 2 && size(values) == (ncell, 3)
        return Matrix(values)
    elseif ndims(values) == 2 && size(values) == (3, ncell)
        return permutedims(values)
    end
    throw_config_error("AMR size dataset must have N, N×3, or 3×N values (N=$(ncell)); got size $(size(values))."; code=:invalid_field_sources)
end

function load_amr_geometry(simu::AbstractString, config::AbstractDict, fallback_file::AbstractString)
    file = _amr_geometry_file(simu, config, fallback_file)
    validation_error = ensure_readable_file(file; expected_exts=collect(HDF5_EXTS))
    validation_error === nothing || throw_config_error(validation_error; code=:invalid_field_sources)

    x = _read_amr_vector(file, _amr_get(config, "x", nothing), "x")
    y = _read_amr_vector(file, _amr_get(config, "y", nothing), "y")
    z = _read_amr_vector(file, _amr_get(config, "z", nothing), "z")
    length(x) == length(y) == length(z) || throw_config_error("AMR x/y/z datasets have different lengths."; code=:invalid_field_sources)
    isempty(x) && throw_config_error("AMR geometry contains no leaf cells."; code=:invalid_field_sources)

    bounds = _amr_bounds(config)
    widths = _amr_widths(file, config, length(x), bounds)
    all(isfinite.(widths)) && all(widths .> 0) || throw_config_error("AMR cell sizes must be finite and positive."; code=:invalid_field_sources)
    centers = hcat(x, y, z)
    tolerance = Float64(_amr_get(config, "tolerance", 1e-8))
    isfinite(tolerance) && tolerance >= 0 || throw_config_error("AMR tolerance must be finite and non-negative."; code=:invalid_field_sources)
    strict = Bool(_amr_get(config, "strict", true))
    return AMRGeometry(centers, widths, bounds, _amr_shape(config), strict, tolerance)
end

"""Rasterize an intensive leaf-cell field by exact volume-overlap averaging."""
function rasterize_amr_field(values::AbstractArray, geometry::AMRGeometry; label::AbstractString="field")
    field = vec(values)
    ncell = size(geometry.centers, 1)
    length(field) == ncell || throw_config_error(
        "AMR $(label) has $(length(field)) values but the geometry has $(ncell) cells.";
        code=:cube_shape_mismatch)
    all(isfinite, field) || throw_config_error("AMR $(label) contains non-finite values."; code=:invalid_field_sources)

    shape = geometry.shape
    voxel_width = ntuple(axis -> (geometry.bounds[axis][2] - geometry.bounds[axis][1]) / shape[axis], 3)
    voxel_volume = prod(voxel_width)
    weighted = zeros(Float64, shape)
    covered = zeros(Float64, shape)

    for cell in 1:ncell
        lower = ntuple(axis -> geometry.centers[cell, axis] - geometry.widths[cell, axis] / 2, 3)
        upper = ntuple(axis -> geometry.centers[cell, axis] + geometry.widths[cell, axis] / 2, 3)
        ranges = ntuple(3) do axis
            lo, hi = geometry.bounds[axis]
            first_index = clamp(floor(Int, (lower[axis] - lo) / voxel_width[axis]) + 1, 1, shape[axis])
            last_index = clamp(ceil(Int, (upper[axis] - lo) / voxel_width[axis]), 1, shape[axis])
            first_index:last_index
        end

        for k in ranges[3], j in ranges[2], i in ranges[1]
            indices = (i, j, k)
            overlap = 1.0
            for axis in 1:3
                voxel_lo = geometry.bounds[axis][1] + (indices[axis] - 1) * voxel_width[axis]
                voxel_hi = voxel_lo + voxel_width[axis]
                overlap *= max(0.0, min(upper[axis], voxel_hi) - max(lower[axis], voxel_lo))
            end
            overlap == 0 && continue
            weighted[i, j, k] += Float64(field[cell]) * overlap
            covered[i, j, k] += overlap
        end
    end

    tolerance_volume = geometry.tolerance * voxel_volume
    uncovered = findfirst(<(voxel_volume - tolerance_volume), covered)
    overlap = findfirst(>(voxel_volume + tolerance_volume), covered)
    if geometry.strict && uncovered !== nothing
        throw_config_error("AMR leaf cells do not cover output voxel $(Tuple(uncovered)) completely."; code=:cube_shape_mismatch)
    end
    if geometry.strict && overlap !== nothing
        throw_config_error("AMR leaf cells overlap in output voxel $(Tuple(overlap)); provide leaf cells only."; code=:cube_shape_mismatch)
    end

    result = fill(NaN, shape)
    valid = covered .> tolerance_volume
    result[valid] .= weighted[valid] ./ covered[valid]
    return result
end

function read_amr_field(source, conversion::Real, geometry::AMRGeometry)
    source isa HDF5DatasetSource || throw_config_error(
        "AMR fields must use HDF5 dataset sources (`path` plus `dataset`).";
        code=:invalid_field_sources)
    values = read_file(source, conversion; expected_ndims=nothing)
    return rasterize_amr_field(values, geometry; label=source_label(source))
end
