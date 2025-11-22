include(joinpath(@__DIR__, "SyntheticObservations", "MOOSE_from_config.jl"))
using .MOOSEFromConfig

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        error("Usage: julia --project=. src/MOOSE_cli.jl <config.json> [--quiet]")
    end

    config_path = ARGS[1]
    quiet = any(arg -> arg == "--quiet", ARGS[2:end])

    MOOSEFromConfig.MOOSE_from_config(config_path; quiet = quiet)
end
