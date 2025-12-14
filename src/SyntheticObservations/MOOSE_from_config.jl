module MOOSEFromConfig

using CSV
using Crayons
using DataFrames
using Dates
using JSON

using ..MOOSE: PARSEC_TO_CM, RunConfig, ValidationResult, ensure_directory_access,
               run_moose_processing, throw_config_error, validation_failure, validation_success

function normalize_base_dir(cfg, config_path)
    raw_dir = get(cfg, "base_dir") do
        return validation_failure(String, "Config file $(config_path) must define `base_dir`.")
    end

    raw_dir isa ValidationResult && return raw_dir

    normalized_dir = abspath(expanduser(String(raw_dir)))
    validation_error = ensure_directory_access(normalized_dir)
    validation_error === nothing || return validation_failure(String, validation_error)

    resolved_dir = try
        realpath(normalized_dir)
    catch
        normalized_dir
    end

    return validation_success(resolved_dir)
end

function collect_simulations(cfg, base_dir, config_path)
    sims = haskey(cfg, "simulations") ? cfg["simulations"] : get(cfg, "chosen_simu") do
        return validation_failure(Vector{String}, "Config file $(config_path) must define either `simulations` or `chosen_simu` (array of simulations).")
    end
    sims isa ValidationResult && return sims
    sims isa AbstractVector || return validation_failure(Vector{String}, "`simulations` must be an array of paths or folder names.")

    paths = String[]
    errors = String[]

    for simu in sims
        candidate = isabspath(simu) ? String(simu) : joinpath(base_dir, String(simu))
        expanded = abspath(expanduser(candidate))

        validation_error = ensure_directory_access(expanded)
        if validation_error !== nothing
            push!(errors, validation_error)
            continue
        end

        resolved = try
            realpath(expanded)
        catch
            expanded
        end

        push!(paths, resolved)
    end

    if !isempty(errors)
        return validation_failure(Vector{String}, join(errors, "\n"))
    end

    return validation_success(paths)
end

function collect_los(cfg)
    los_values = get(cfg, "chosen_LOS", ["x", "y", "z"])
    los_values isa AbstractVector ||
        throw_config_error("`chosen_LOS` must be an array of LOS identifiers (x, y, z)."; code=:invalid_los)

    normalized = unique(lowercase.(String.(los_values)))
    invalid = filter(los -> !(los in ("x", "y", "z")), normalized)
    !isempty(invalid) &&
        throw_config_error("Invalid line(s) of sight: $(join(invalid, ", ")). Allowed values are x, y, or z."; code=:invalid_los)

    return normalized
end

function build_frequency_array(cfg)
    freq_cfg = get(cfg, "freq", nothing)
    if freq_cfg isa AbstractDict
        start_val = get(freq_cfg, "start", 120.0)
        end_val = get(freq_cfg, "end", 167.0)
        step_val = get(freq_cfg, "step", 0.098)
    else
        start_val = get(cfg, "nustart", 120.0)
        end_val = get(cfg, "nuend", 167.0)
        step_val = get(cfg, "dnu", 0.098)
    end

    step_val <= 0 && throw_config_error("Frequency step must be positive (received $(step_val))."; code=:invalid_frequency)
    end_val <= start_val && throw_config_error(
        "The end frequency ($(end_val)) must be greater than the start frequency ($(start_val)).";
        code=:invalid_frequency,
    )

    return Float64(start_val), Float64(end_val), Float64(step_val)
end

function build_faraday(cfg)
    faraday_cfg = get(cfg, "faraday", nothing)
    if faraday_cfg isa AbstractDict
        enabled = get(faraday_cfg, "enabled", false)
        phimin = get(faraday_cfg, "phimin", -20.0)
        phimax = get(faraday_cfg, "phimax", 20.0)
        dphi = get(faraday_cfg, "dphi", 0.1)
        return enabled ? "Y" : "N", Float64(phimin), Float64(phimax), Float64(dphi)
    else
        rotation_flag = uppercase(get(cfg, "FaradayRotation", "N"))
        phimin = get(cfg, "phimin", -20.0)
        phimax = get(cfg, "phimax", 20.0)
        dphi = get(cfg, "dphi", 0.1)
        return rotation_flag, Float64(phimin), Float64(phimax), Float64(dphi)
    end

    dphi <= 0 && throw_config_error("The Faraday step dphi must be positive (received $(dphi))."; code=:invalid_faraday_range)
    phimax <= phimin && throw_config_error(
        "Faraday rotation range is invalid: phimax ($(phimax)) must be greater than phimin ($(phimin)).";
        code=:invalid_faraday_range,
    )
end

function normalize_box_lengths(box_length)
    if box_length isa AbstractVector
        length(box_length) == 3 || throw_config_error("Box length array must have three elements (x, y, z).")
        return (; x = Float64(box_length[1]), y = Float64(box_length[2]), z = Float64(box_length[3]))
    elseif box_length isa AbstractDict
        return (; x = Float64(get(box_length, "x", get(box_length, "X", get(box_length, "size_pc", 50.0)))),
                 y = Float64(get(box_length, "y", get(box_length, "Y", get(box_length, "size_pc", 50.0)))),
                 z = Float64(get(box_length, "z", get(box_length, "Z", get(box_length, "size_pc", 50.0)))))
    else
        return (; x = Float64(box_length), y = Float64(box_length), z = Float64(box_length))
    end
end

function normalize_box_pixels(box_length_pix)
    if box_length_pix isa AbstractVector
        length(box_length_pix) == 3 || throw_config_error("Box pixel array must have three elements (x, y, z).")
        return (; x = Int(box_length_pix[1]), y = Int(box_length_pix[2]), z = Int(box_length_pix[3]))
    elseif box_length_pix isa AbstractDict
        return (; x = Int(get(box_length_pix, "x", get(box_length_pix, "X", get(box_length_pix, "npix", 256)))),
                 y = Int(get(box_length_pix, "y", get(box_length_pix, "Y", get(box_length_pix, "npix", 256)))),
                 z = Int(get(box_length_pix, "z", get(box_length_pix, "Z", get(box_length_pix, "npix", 256)))))
    else
        return (; x = Int(box_length_pix), y = Int(box_length_pix), z = Int(box_length_pix))
    end
end

function build_distance_parameters(cfg)
    box_cfg = get(cfg, "box", nothing)
    raw_box_length_pc = box_cfg isa AbstractDict ? get(box_cfg, "size_pc", get(cfg, "BoxLength_pc", 50.0)) : get(cfg, "BoxLength_pc", 50.0)
    raw_box_length_pix = box_cfg isa AbstractDict ? get(box_cfg, "npix", get(cfg, "BoxLength_pix", 256)) : get(cfg, "BoxLength_pix", 256)

    box_length_pc = normalize_box_lengths(raw_box_length_pc)
    box_length_pix = normalize_box_pixels(raw_box_length_pix)

    return Float64(box_length_pc.x), Int(box_length_pix.x)
end

"""
    MOOSE_from_config(config_path::AbstractString; quiet::Bool = false)

Run the MOOSE pipeline non-interactively using parameters defined in a JSON configuration
file. The configuration keys mirror those written by the interactive `MOOSE` function
(e.g. `base_dir`, `chosen_simu`, `chosen_LOS`, `conversionB`, `FaradayRotation`, etc.),
and also supports a streamlined schema with nested keys as produced by the Streamlit
frontend example (`freq`, `box`, `faraday`, `ne`, `emissivity`).

# Arguments
- `config_path`: path to the JSON configuration file.
- `quiet`: when `true`, skips the startup logo from `MOOSE`.

The function will raise an error if required keys are missing from the configuration.
"""
function build_config(cfg, config_path)
    base_dir_result = normalize_base_dir(cfg, config_path)
    base_dir_result.error === nothing || throw_config_error(base_dir_result.error; code=:missing_base_dir)
    base_dir = base_dir_result.value

    simu_paths_result = collect_simulations(cfg, base_dir, config_path)
    simu_paths_result.error === nothing || throw_config_error(simu_paths_result.error; code=:missing_simulation)
    simu_paths = simu_paths_result.value

    chosen_LOS = collect_los(cfg)
    conversionB = get(cfg, "conversionB", 1.0)
    conversionn = get(cfg, "conversionn", 1.0)
    conversionT = get(cfg, "conversionT", 1.0)
    log_progress = get(cfg, "log_progress", false)
    responseSynchrotron = uppercase(get(cfg, "responseSynchrotron", "N"))
    kernel_size_synchrotron = get(cfg, "kernel_size_synchrotron", nothing)
    add_noise = uppercase(get(cfg, "add_noise", "N"))
    SNR_nu = get(cfg, "SNR_nu", nothing)

    interpolation_file_path = get(cfg, "interpolation_file_path") do
        emiss_cfg = get(cfg, "emissivity", nothing)
        emiss_cfg isa AbstractDict ? get(emiss_cfg, "path") : nothing
    end
    interpolation_file_path === nothing && throw_config_error(
        "Config file $(config_path) must define `interpolation_file_path` or `emissivity.path`.";
        code=:missing_interpolation_path,
    )
    interpolation_file_path = isabspath(interpolation_file_path) ? interpolation_file_path : joinpath(base_dir, interpolation_file_path)

    ne_option = string(get(cfg, "ne_option", get(get(cfg, "ne", Dict()), "mode", 1)))
    IonizationFraction = get(cfg, "IonizationFraction", get(get(cfg, "ne", Dict()), "ion_fraction", 0.01))

    BoxLength_pc, BoxLength_pix = build_distance_parameters(cfg)
    nustart, nuend, dnu = build_frequency_array(cfg)
    FaradayRotation, phimin, phimax, dphi = build_faraday(cfg)

    return RunConfig(
        base_dir,
        simu_paths,
        chosen_LOS,
        conversionB,
        conversionn,
        conversionT,
        FaradayRotation,
        phimin,
        phimax,
        dphi,
        responseSynchrotron,
        kernel_size_synchrotron,
        add_noise,
        SNR_nu,
        interpolation_file_path,
        ne_option,
        IonizationFraction,
        nustart,
        nuend,
        dnu,
        BoxLength_pc,
        BoxLength_pix,
        config_path,
    ), simu_paths
end

function MOOSE_from_config(config_path::AbstractString; quiet::Bool = false)
    cfg = JSON.parsefile(config_path)

    run_config, _ = build_config(cfg, config_path)
    run_moose_processing(run_config; quiet = quiet, persisted_config = cfg)

    return nothing
end

end # module
