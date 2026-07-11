using JSON
using Moose: MooseError, cli_error, throw_cli_error, preflight_plan
using Moose.MooseFromConfig: MOOSE_from_config_dict, build_config

function parse_numeric(name, value)
    parsed = tryparse(Float64, value)
    parsed === nothing && throw_cli_error("$(name) expects a numeric value, got '$(value)'.")
    return parsed
end

function parse_flag(name, value; allowed=("Y", "N"))
    normalized = uppercase(strip(String(value)))
    if normalized in ("Y", "YES", "TRUE", "1")
        flag = "Y"
    elseif normalized in ("N", "NO", "FALSE", "0")
        flag = "N"
    else
        throw_cli_error("$(name) expects $(join(allowed, " or ")), got '$(value)'.")
    end
    flag in allowed || throw_cli_error("$(name) expects $(join(allowed, " or ")), got '$(value)'.")
    return flag
end

function parse_ne_option(value)
    normalized = string(value)
    normalized in ("1", "2", "3") || throw_cli_error("--ne-option expects 1, 2, or 3, got '$(value)'.")
    return normalized
end

function parse_density_kind(value)
    normalized = lowercase(strip(String(value)))
    normalized in ("number_density", "mass_density") ||
        throw_cli_error("--density-kind expects number_density or mass_density, got '$(value)'.")
    return normalized
end

function parse_integer(name, value)
    parsed = tryparse(Int, value)
    parsed === nothing && throw_cli_error("$(name) expects an integer value, got '$(value)'.")
    return parsed
end

function parse_cli_args(args)
    config_path = nothing
    quiet = false
    write_back = false
    overrides = Dict{String, Any}()
    simulations = String[]
    los_values = String[]
    interpolation_file = nothing

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            i += 1
            i > length(args) && throw_cli_error("--config expects a path")
            config_path = args[i]
        elseif arg == "--base-dir"
            i += 1
            i > length(args) && throw_cli_error("--base-dir expects a path")
            overrides["base_dir"] = args[i]
        elseif arg == "--simu"
            i += 1
            i > length(args) && throw_cli_error("--simu expects a simulation path")
            push!(simulations, args[i])
        elseif arg == "--los"
            i += 1
            i > length(args) && throw_cli_error("--los expects x, y, or z")
            raw_los = lowercase(args[i])
            if raw_los == "all"
                los_values = ["x", "y", "z"]
            else
                for candidate in split(raw_los, ",")
                    los = strip(candidate)
                    isempty(los) && continue
                    los in ("x", "y", "z") || throw_cli_error("Invalid LOS '$(los)'. Use x, y, z, or 'all'.")
                    push!(los_values, los)
                end
            end
        elseif arg == "--interpolation"
            i += 1
            i > length(args) && throw_cli_error("--interpolation expects a file path")
            interpolation_file = args[i]
        elseif arg == "--conversionB"
            i += 1
            i > length(args) && throw_cli_error("--conversionB expects a value")
            overrides["conversionB"] = parse_numeric("--conversionB", args[i])
        elseif arg == "--conversionn"
            i += 1
            i > length(args) && throw_cli_error("--conversionn expects a value")
            overrides["conversionn"] = parse_numeric("--conversionn", args[i])
        elseif arg == "--conversionT"
            i += 1
            i > length(args) && throw_cli_error("--conversionT expects a value")
            overrides["conversionT"] = parse_numeric("--conversionT", args[i])
        elseif arg == "--density-kind"
            i += 1
            i > length(args) && throw_cli_error("--density-kind expects number_density or mass_density")
            overrides["density_kind"] = parse_density_kind(args[i])
        elseif arg == "--mean-molecular-weight"
            i += 1
            i > length(args) && throw_cli_error("--mean-molecular-weight expects a value")
            overrides["mean_molecular_weight"] = parse_numeric("--mean-molecular-weight", args[i])
        elseif arg == "--hydrogen-mass-g"
            i += 1
            i > length(args) && throw_cli_error("--hydrogen-mass-g expects a value")
            overrides["hydrogen_mass_g"] = parse_numeric("--hydrogen-mass-g", args[i])
        elseif arg == "--faraday"
            i += 1
            i > length(args) && throw_cli_error("--faraday expects Y or N")
            overrides["FaradayRotation"] = parse_flag("--faraday", args[i])
        elseif arg == "--phimin"
            i += 1
            i > length(args) && throw_cli_error("--phimin expects a value")
            overrides["phimin"] = parse_numeric("--phimin", args[i])
        elseif arg == "--phimax"
            i += 1
            i > length(args) && throw_cli_error("--phimax expects a value")
            overrides["phimax"] = parse_numeric("--phimax", args[i])
        elseif arg == "--dphi"
            i += 1
            i > length(args) && throw_cli_error("--dphi expects a value")
            overrides["dphi"] = parse_numeric("--dphi", args[i])
        elseif arg == "--filtering"
            i += 1
            i > length(args) && throw_cli_error("--filtering expects Y or N")
            overrides["responseSynchrotron"] = parse_flag("--filtering", args[i])
        elseif arg == "--kernel-size"
            i += 1
            i > length(args) && throw_cli_error("--kernel-size expects a value")
            overrides["kernel_size_synchrotron"] = parse_numeric("--kernel-size", args[i])
        elseif arg == "--noise"
            i += 1
            i > length(args) && throw_cli_error("--noise expects Y or N")
            overrides["add_noise"] = parse_flag("--noise", args[i])
        elseif arg == "--snr"
            i += 1
            i > length(args) && throw_cli_error("--snr expects a value")
            overrides["SNR_nu"] = parse_numeric("--snr", args[i])
        elseif arg == "--rng-seed"
            i += 1
            i > length(args) && throw_cli_error("--rng-seed expects an integer value")
            overrides["rng_seed"] = parse_integer("--rng-seed", args[i])
        elseif arg == "--ne-option"
            i += 1
            i > length(args) && throw_cli_error("--ne-option expects 1, 2, or 3")
            overrides["ne_option"] = parse_ne_option(args[i])
        elseif arg == "--zeta"
            i += 1
            i > length(args) && throw_cli_error("--zeta expects a value")
            overrides["zeta"] = parse_numeric("--zeta", args[i])
        elseif arg == "--Geff"
            i += 1
            i > length(args) && throw_cli_error("--Geff expects a value")
            overrides["Geff"] = parse_numeric("--Geff", args[i])
        elseif arg == "--phiPAH"
            i += 1
            i > length(args) && throw_cli_error("--phiPAH expects a value")
            overrides["phiPAH"] = parse_numeric("--phiPAH", args[i])
        elseif arg == "--XC"
            i += 1
            i > length(args) && throw_cli_error("--XC expects a value")
            overrides["XC"] = parse_numeric("--XC", args[i])
        elseif arg == "--precision"
            i += 1
            i > length(args) && throw_cli_error("--precision expects float64 or float32")
            precision_value = lowercase(strip(args[i]))
            precision_value in ("float64", "float32") || throw_cli_error("--precision expects float64 or float32, got '$(args[i])'.")
            overrides["precision"] = precision_value
        elseif arg == "--tile-size"
            i += 1
            i > length(args) && throw_cli_error("--tile-size expects a positive integer")
            tile_value = parse_integer("--tile-size", args[i])
            tile_value > 0 || throw_cli_error("--tile-size expects a positive integer, got '$(args[i])'.")
            overrides["tile_size"] = tile_value
        elseif arg == "--resume"
            i += 1
            i > length(args) && throw_cli_error("--resume expects off or safe")
            resume_value = lowercase(strip(args[i]))
            resume_value in ("off", "safe") || throw_cli_error("--resume expects off or safe, got '$(args[i])'.")
            overrides["resume"] = resume_value
        elseif arg == "--outputs"
            i += 1
            i > length(args) && throw_cli_error("--outputs expects a comma-separated list")
            output_values = [strip(lowercase(value)) for value in split(args[i], ",") if !isempty(strip(value))]
            isempty(output_values) && throw_cli_error("--outputs expects at least one output group")
            overrides["outputs"] = output_values
        elseif arg == "--quiet"
            quiet = true
        elseif arg == "--write-back"
            write_back = true
        elseif arg == "--plan"
            overrides["__plan_only"] = true
        else
            if !startswith(arg, "--") && config_path === nothing
                config_path = arg
            else
                throw_cli_error("Unknown argument: $(arg)")
            end
        end
        i += 1
    end

    isempty(simulations) || (overrides["simulations"] = simulations)
    isempty(los_values) || (overrides["chosen_LOS"] = los_values)
    interpolation_file === nothing || (overrides["interpolation_file_path"] = interpolation_file)

    write_back && config_path === nothing && throw_cli_error("--write-back requires a config path (positional or via --config).")
    config_path !== nothing && !isfile(config_path) && throw_cli_error("Config file not found: $(config_path)")

    return config_path, quiet, write_back, overrides
end

function load_base_config(config_path)
    if config_path !== nothing && isfile(config_path)
        return JSON.parsefile(config_path)
    else
        return Dict{String, Any}()
    end
end

function run_with_config(config_path, quiet, write_back, overrides)
    cfg = merge(load_base_config(config_path), overrides)
    plan_only = pop!(cfg, "__plan_only", false)

    if plan_only
        effective_path = config_path === nothing ? "<cli-overrides>" : "$(config_path) + CLI overrides"
        run_config, _ = build_config(cfg, effective_path)
        preflight_plan(run_config)
        return nothing
    end

    if write_back
        Moose.save_config(cfg, config_path)
        MOOSE_from_config_dict(
            cfg;
            config_path = config_path,
            quiet = quiet,
            source_config_path = config_path,
            saved_config_path = config_path,
            write_config_file = true,
        )
    else
        effective_path = config_path === nothing ? "<cli-overrides>" : "$(config_path) + CLI overrides"
        MOOSE_from_config_dict(
            cfg;
            config_path = effective_path,
            quiet = quiet,
            source_config_path = config_path,
            saved_config_path = nothing,
            write_config_file = false,
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        config_path, quiet, write_back, overrides = parse_cli_args(ARGS)
        run_with_config(config_path, quiet, write_back, overrides)
    catch err
        if err isa MooseError
            println(stderr, err.message)
            exit(err.exit_code)
        else
            rethrow()
        end
    end
end
