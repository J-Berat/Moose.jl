module MOOSEFromConfig

using CSV
using Crayons
using DataFrames
using Dates
using JSON

include("MOOSE.jl")
using ..MOOSE: WolfireConstants, ProcessSynchrotron, print_logo, save_config, write_summary_log

const PARSEC_TO_CM = 3.0857e18

function normalize_base_dir(cfg, config_path)
    base_dir = get(cfg, "base_dir") do
        error("Config file $(config_path) must define `base_dir`.")
    end
    isdir(base_dir) || error("Configured base_dir $(base_dir) does not exist.")
    return base_dir
end

function collect_simulations(cfg, base_dir, config_path)
    sims = haskey(cfg, "simulations") ? cfg["simulations"] : get(cfg, "chosen_simu") do
        error("Config file $(config_path) must define either `simulations` or `chosen_simu` (array of simulations).")
    end
    sims isa AbstractVector || error("`simulations` must be an array of paths or folder names.")
    return map(sims) do simu
        isabspath(simu) ? simu : joinpath(base_dir, simu)
    end
end

function collect_los(cfg)
    return get(cfg, "chosen_LOS", ["x", "y", "z"])
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

    return range(start = Float64(start_val), stop = Float64(end_val), step = Float64(step_val))
end

function build_faraday(cfg)
    faraday_cfg = get(cfg, "faraday", nothing)
    if faraday_cfg isa AbstractDict
        enabled = get(faraday_cfg, "enabled", false)
        phimin = get(faraday_cfg, "phimin", -20.0)
        phimax = get(faraday_cfg, "phimax", 20.0)
        dphi = get(faraday_cfg, "dphi", 0.1)
        return enabled ? "Y" : "N", range(start = Float64(phimin), stop = Float64(phimax), step = Float64(dphi))
    else
        rotation_flag = uppercase(get(cfg, "FaradayRotation", "N"))
        phimin = get(cfg, "phimin", -20.0)
        phimax = get(cfg, "phimax", 20.0)
        dphi = get(cfg, "dphi", 0.1)
        phi_array = rotation_flag == "Y" ? range(start = Float64(phimin), stop = Float64(phimax), step = Float64(dphi)) : nothing
        return rotation_flag, phi_array
    end
end

function build_distance_parameters(cfg)
    box_cfg = get(cfg, "box", nothing)
    if box_cfg isa AbstractDict
        box_length_pc = get(box_cfg, "size_pc", 50.0)
        box_length_pix = get(box_cfg, "npix", 256)
    else
        box_length_pc = get(cfg, "BoxLength_pc", 50.0)
        box_length_pix = get(cfg, "BoxLength_pix", 256)
    end

    pixel_length_pc = Float64(box_length_pc) / Float64(box_length_pix)
    pixel_length_cm = pixel_length_pc * PARSEC_TO_CM
    distance_array = range(start = 0.0, stop = Float64(box_length_pc), step = pixel_length_pc)

    return pixel_length_pc, pixel_length_cm, Float64(box_length_pc), distance_array
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
function MOOSE_from_config(config_path::AbstractString; quiet::Bool = false)
    cfg = JSON.parsefile(config_path)

    base_dir = normalize_base_dir(cfg, config_path)
    simu_paths = collect_simulations(cfg, base_dir, config_path)

    chosen_LOS = collect_los(cfg)
    conversionB = get(cfg, "conversionB", 1.0)
    conversionn = get(cfg, "conversionn", 1.0)
    conversionT = get(cfg, "conversionT", 1.0)
    responseSynchrotron = uppercase(get(cfg, "responseSynchrotron", "N"))
    kernel_size_synchrotron = get(cfg, "kernel_size_synchrotron", nothing)
    add_noise = uppercase(get(cfg, "add_noise", "N"))
    SNR_nu = get(cfg, "SNR_nu", nothing)

    interpolation_file_path = get(cfg, "interpolation_file_path") do
        emiss_cfg = get(cfg, "emissivity", nothing)
        emiss_cfg isa AbstractDict ? get(emiss_cfg, "path") : nothing
    end
    interpolation_file_path === nothing && error("Config file $(config_path) must define `interpolation_file_path` or `emissivity.path`.")
    interpolation_file_path = isabspath(interpolation_file_path) ? interpolation_file_path : joinpath(base_dir, interpolation_file_path)
    isfile(interpolation_file_path) || error("Interpolation file $(interpolation_file_path) not found.")

    ne_option = string(get(cfg, "ne_option", get(get(cfg, "ne", Dict()), "mode", 1)))
    IonizationFraction = get(cfg, "IonizationFraction", get(get(cfg, "ne", Dict()), "ion_fraction", 0.01))

    PixelLength_pc, PixelLength_cm, BoxLength_pc, DistanceArray = build_distance_parameters(cfg)
    nuArray = build_frequency_array(cfg)
    FaradayRotation, PhiArray = build_faraday(cfg)
    df = CSV.File(interpolation_file_path) |> DataFrame

    if !quiet
        print_logo()
    end

    start_time = now()

    if ne_option == "3"
        missing_cubes = [simu for simu in simu_paths if !isfile(joinpath(simu, "densityHp.fits"))]
        !isempty(missing_cubes) && error("Electron density cube 'densityHp.fits' is missing for: $(join(missing_cubes, ", ")).")
    end

    for (i, simu_path) in enumerate(simu_paths)
        simu_name = basename(simu_path)
        println("Processing Simulation: $(simu_name)")

        for LOS in chosen_LOS
            println(Crayon(foreground = :yellow, bold = true)("→ Processing LOS: $(LOS)"))

            if ne_option == "1"
                zeta, Geff, omegaPAH, XC = WolfireConstants()
                ProcessSynchrotron(simu_path, LOS, FaradayRotation, responseSynchrotron, df, add_noise, SNR_nu,
                    kernel_size_synchrotron, zeta, Geff, omegaPAH, XC, nuArray, PhiArray, PixelLength_pc, PixelLength_cm,
                    BoxLength_pc, DistanceArray, conversionn, conversionT, conversionB)
            elseif ne_option == "2"
                ProcessSynchrotron(simu_path, LOS, FaradayRotation, responseSynchrotron, df, add_noise, SNR_nu,
                    kernel_size_synchrotron, IonizationFraction, nuArray, PhiArray, PixelLength_pc, PixelLength_cm,
                    BoxLength_pc, DistanceArray, conversionn, conversionT, conversionB)
            else
                ProcessSynchrotron(simu_path, LOS, FaradayRotation, responseSynchrotron, df, add_noise, SNR_nu, kernel_size_synchrotron,
                    nuArray, PhiArray, PixelLength_pc, PixelLength_cm, BoxLength_pc, DistanceArray,
                    conversionn, conversionT, conversionB)
            end
        end

        if length(simu_paths) > 1
            println("Finished processing all chosen LOS for simulation: $simu_name")
            print_progress(i, length(simu_paths))
        end
    end

    println("Finished processing all simulations.")

    elapsed = now() - start_time
    println(Crayon(foreground = :green, bold = true)("Summary:"))
    println(Crayon(foreground = :green)("Simulations processed: $(join(map(basename, simu_paths), ", "))"))
    println(Crayon(foreground = :green)("Lines of sight: $(join(chosen_LOS, ", "))"))
    println(Crayon(foreground = :green)("Output directory: $base_dir"))
    println(Crayon(foreground = :green)("Total execution time: $(format_duration(elapsed))"))

    save_config(cfg, config_path)
    write_summary_log(base_dir, map(basename, simu_paths), chosen_LOS, elapsed; config_path=config_path, faraday=FaradayRotation, responseSynchrotron=responseSynchrotron, add_noise=add_noise, interpolation_file_path=interpolation_file_path, conversionB=conversionB, conversionn=conversionn, conversionT=conversionT, ne_option=ne_option)

    return nothing
end

end # module
