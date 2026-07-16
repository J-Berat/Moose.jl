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

"""Precomputed, exact assignment of every regular output voxel to one AMR leaf cell."""
struct AMRRasterPlan
    geometry::AMRGeometry
    cell_index::Array{Int, 3} # zero marks an uncovered voxel when strict=false
end

const _AMR_PLAN_CACHE = Dict{Any, AMRRasterPlan}()
const _AMR_PLAN_CACHE_LOCK = ReentrantLock()
# One plan is enough to reuse geometry across the sequential x/y/z workflow,
# while keeping memory bounded for very large target grids.
const _AMR_PLAN_CACHE_LIMIT = 1

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

function _amr_bool(value, label::AbstractString)
    value isa Bool && return value
    normalized = lowercase(strip(String(value)))
    normalized in ("true", "yes", "y", "1") && return true
    normalized in ("false", "no", "n", "0") && return false
    throw_config_error("`$(label)` must be a boolean; got $(value)."; code=:invalid_field_sources)
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
    strict = _amr_bool(_amr_get(config, "strict", true), "field_sources.amr.strict")
    return AMRGeometry(centers, widths, bounds, _amr_shape(config), strict, tolerance)
end

function build_amr_raster_plan(geometry::AMRGeometry)
    shape = geometry.shape
    voxel_width = ntuple(axis -> (geometry.bounds[axis][2] - geometry.bounds[axis][1]) / shape[axis], 3)
    cell_index = zeros(Int, shape)

    for cell in axes(geometry.centers, 1)
        ranges = ntuple(3) do axis
            bound_lo, bound_hi = geometry.bounds[axis]
            cell_lo = geometry.centers[cell, axis] - geometry.widths[cell, axis] / 2
            cell_hi = geometry.centers[cell, axis] + geometry.widths[cell, axis] / 2
            scale = voxel_width[axis]
            # Alignment is checked in voxel units. Account for floating-point
            # coordinate noise without letting a domain-relative tolerance
            # grow to a sizeable fraction of a voxel on extremely fine grids.
            atol = max(geometry.tolerance, 64 * eps(Float64) * shape[axis])
            qlo = (cell_lo - bound_lo) / scale
            qhi = (cell_hi - bound_lo) / scale

            cell_lo >= bound_lo - geometry.tolerance * (bound_hi - bound_lo) &&
                cell_hi <= bound_hi + geometry.tolerance * (bound_hi - bound_lo) ||
                throw_config_error(
                    "AMR cell $(cell) extends outside configured bounds on axis $(axis).";
                    code=:cube_shape_mismatch)
            geometry.widths[cell, axis] + geometry.tolerance * scale >= scale ||
                throw_config_error(
                    "The AMR target grid is coarser than cell $(cell) on axis $(axis). " *
                    "Increase `field_sources.amr.shape` to at least the finest AMR resolution.";
                    code=:amr_resolution_too_coarse)
            isapprox(qlo, round(qlo); atol=atol, rtol=0) &&
                isapprox(qhi, round(qhi); atol=atol, rtol=0) ||
                throw_config_error(
                    "AMR cell $(cell) is not aligned with the target grid on axis $(axis). " *
                    "Choose a shape matching the finest AMR level to avoid averaging physical states before emissivity calculation.";
                    code=:amr_grid_misaligned)

            first_index = round(Int, qlo) + 1
            last_index = round(Int, qhi)
            first_index <= last_index || throw_config_error(
                "The AMR target grid is coarser than cell $(cell) on axis $(axis). " *
                "Increase `field_sources.amr.shape` to at least the finest AMR resolution.";
                code=:amr_resolution_too_coarse)
            clamp(first_index, 1, shape[axis]):clamp(last_index, 1, shape[axis])
        end

        for k in ranges[3], j in ranges[2], i in ranges[1]
            previous = cell_index[i, j, k]
            previous == 0 || throw_config_error(
                "AMR leaf cells $(previous) and $(cell) overlap in output voxel $((i, j, k)); provide leaf cells only.";
                code=:cube_shape_mismatch)
            cell_index[i, j, k] = cell
        end
    end

    if geometry.strict
        uncovered = findfirst(==(0), cell_index)
        uncovered === nothing || throw_config_error(
            "AMR leaf cells do not cover output voxel $(Tuple(uncovered)) completely.";
            code=:cube_shape_mismatch)
    end
    return AMRRasterPlan(geometry, cell_index)
end

function _amr_plan_cache_key(simu::AbstractString, config::AbstractDict, fallback_file::AbstractString)
    file = abspath(_amr_geometry_file(simu, config, fallback_file))
    validation_error = ensure_readable_file(file; expected_exts=collect(HDF5_EXTS))
    validation_error === nothing || throw_config_error(validation_error; code=:invalid_field_sources)
    info = stat(file)
    geometry_keys = ("x", "y", "z", "size", "level", "level_offset", "strict", "tolerance")
    settings = Tuple(string(_amr_get(config, key, nothing)) for key in geometry_keys)
    return (file, info.size, info.mtime, _amr_shape(config), _amr_bounds(config), settings)
end

"""Load and cache the geometry-to-voxel assignment for reuse by every field and LOS."""
function load_amr_raster_plan(simu::AbstractString, config::AbstractDict, fallback_file::AbstractString)
    key = _amr_plan_cache_key(simu, config, fallback_file)
    cached = lock(_AMR_PLAN_CACHE_LOCK) do
        get(_AMR_PLAN_CACHE, key, nothing)
    end
    cached === nothing || return cached

    plan = build_amr_raster_plan(load_amr_geometry(simu, config, fallback_file))
    return lock(_AMR_PLAN_CACHE_LOCK) do
        if length(_AMR_PLAN_CACHE) >= _AMR_PLAN_CACHE_LIMIT
            delete!(_AMR_PLAN_CACHE, first(keys(_AMR_PLAN_CACHE)))
        end
        get!(_AMR_PLAN_CACHE, key, plan)
    end
end

"""Rasterize a field without mixing distinct AMR leaf-cell states."""
function rasterize_amr_field(values::AbstractArray, plan::AMRRasterPlan; label::AbstractString="field")
    field = vec(values)
    ncell = size(plan.geometry.centers, 1)
    length(field) == ncell || throw_config_error(
        "AMR $(label) has $(length(field)) values but the geometry has $(ncell) cells.";
        code=:cube_shape_mismatch)
    all(isfinite, field) || throw_config_error("AMR $(label) contains non-finite values."; code=:invalid_field_sources)

    result = fill(NaN, plan.geometry.shape)
    @inbounds for index in eachindex(plan.cell_index)
        cell = plan.cell_index[index]
        cell == 0 || (result[index] = Float64(field[cell]))
    end
    return result
end

function rasterize_amr_field(values::AbstractArray, geometry::AMRGeometry; kwargs...)
    return rasterize_amr_field(values, build_amr_raster_plan(geometry); kwargs...)
end

function read_amr_field(source, conversion::Real, plan::AMRRasterPlan)
    source isa HDF5DatasetSource || throw_config_error(
        "AMR fields must use HDF5 dataset sources (`path` plus `dataset`).";
        code=:invalid_field_sources)
    values = read_file(source, conversion; expected_ndims=nothing)
    return rasterize_amr_field(values, plan; label=source_label(source))
end
