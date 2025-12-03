using JSON
using MOOSE

function parse_numeric(name, value)
    parsed = tryparse(Float64, value)
    parsed === nothing && error("$(name) expects a numeric value, got '$(value)'.")
    return parsed
end

function parse_cli_args(args)
    config_path = nothing
    quiet = false
    overrides = Dict{String, Any}()
    simulations = String[]
    los_values = String[]
    interpolation_file = nothing

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            i += 1
            i > length(args) && error("--config expects a path")
            config_path = args[i]
        elseif arg == "--base-dir"
            i += 1
            i > length(args) && error("--base-dir expects a path")
            overrides["base_dir"] = args[i]
        elseif arg == "--simu"
            i += 1
            i > length(args) && error("--simu expects a simulation path")
            push!(simulations, args[i])
        elseif arg == "--los"
            i += 1
            i > length(args) && error("--los expects x, y, or z")
            raw_los = lowercase(args[i])
            if raw_los == "all"
                los_values = ["x", "y", "z"]
            else
                for candidate in split(raw_los, ",")
                    los = strip(candidate)
                    isempty(los) && continue
                    los in ("x", "y", "z") || error("Invalid LOS '$(los)'. Use x, y, z, or 'all'.")
                    push!(los_values, los)
                end
            end
        elseif arg == "--interpolation"
            i += 1
            i > length(args) && error("--interpolation expects a file path")
            interpolation_file = args[i]
        elseif arg == "--conversionB"
            i += 1
            i > length(args) && error("--conversionB expects a value")
            overrides["conversionB"] = parse_numeric("--conversionB", args[i])
        elseif arg == "--conversionn"
            i += 1
            i > length(args) && error("--conversionn expects a value")
            overrides["conversionn"] = parse_numeric("--conversionn", args[i])
        elseif arg == "--conversionT"
            i += 1
            i > length(args) && error("--conversionT expects a value")
            overrides["conversionT"] = parse_numeric("--conversionT", args[i])
        elseif arg == "--faraday"
            i += 1
            i > length(args) && error("--faraday expects Y or N")
            overrides["FaradayRotation"] = uppercase(args[i])
        elseif arg == "--phimin"
            i += 1
            i > length(args) && error("--phimin expects a value")
            overrides["phimin"] = parse_numeric("--phimin", args[i])
        elseif arg == "--phimax"
            i += 1
            i > length(args) && error("--phimax expects a value")
            overrides["phimax"] = parse_numeric("--phimax", args[i])
        elseif arg == "--dphi"
            i += 1
            i > length(args) && error("--dphi expects a value")
            overrides["dphi"] = parse_numeric("--dphi", args[i])
        elseif arg == "--filtering"
            i += 1
            i > length(args) && error("--filtering expects Y or N")
            overrides["responseSynchrotron"] = uppercase(args[i])
        elseif arg == "--kernel-size"
            i += 1
            i > length(args) && error("--kernel-size expects a value")
            overrides["kernel_size_synchrotron"] = parse_numeric("--kernel-size", args[i])
        elseif arg == "--noise"
            i += 1
            i > length(args) && error("--noise expects Y or N")
            overrides["add_noise"] = uppercase(args[i])
        elseif arg == "--snr"
            i += 1
            i > length(args) && error("--snr expects a value")
            overrides["SNR_nu"] = parse_numeric("--snr", args[i])
        elseif arg == "--ne-option"
            i += 1
            i > length(args) && error("--ne-option expects 1, 2, or 3")
            overrides["ne_option"] = string(args[i])
        elseif arg == "--quiet"
            quiet = true
        else
            if !startswith(arg, "--") && config_path === nothing
                config_path = arg
            else
                error("Unknown argument: $(arg)")
            end
        end
        i += 1
    end

    isempty(simulations) || (overrides["simulations"] = simulations)
    isempty(los_values) || (overrides["chosen_LOS"] = los_values)
    interpolation_file === nothing || (overrides["interpolation_file_path"] = interpolation_file)

    return config_path, quiet, overrides
end

function load_base_config(config_path)
    if config_path !== nothing && isfile(config_path)
        return JSON.parsefile(config_path)
    else
        return Dict{String, Any}()
    end
end

function run_with_config(config_path, quiet, overrides)
    cfg = merge(load_base_config(config_path), overrides)

    if config_path === nothing
        mktemp(; cleanup = false) do path, io
            try
                write(io, JSON.json(cfg))
                close(io)
                MOOSE.MOOSE_from_config(path; quiet = quiet)
            finally
                isopen(io) && close(io)
                rm(path; force = true)
            end
        end
    else
        open(config_path, "w") do io
            write(io, JSON.json(cfg))
        end
        MOOSE.MOOSE_from_config(config_path; quiet = quiet)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    config_path, quiet, overrides = parse_cli_args(ARGS)

    run_with_config(config_path, quiet, overrides)
end
