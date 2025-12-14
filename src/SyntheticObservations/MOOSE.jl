"""
    run_moose()

Interactive tool to process synchrotron data for a set of simulations.

# Description
`MOOSE` is an interactive Julia-based tool designed to guide the user through processing synchrotron data for a set of simulations. The function allows the user to configure various parameters, including:
- The base directory of simulations.
- Specific simulations and lines of sight to process.
- Synchrotron data processing options, such as Faraday rotation and filtering.
- Conversion factors for physical units.

It processes the chosen simulations and lines of sight, performing calculations and saving results for later analysis.

# Pre-requisites
- The function assumes the presence of certain dependencies such as `CSV` and `DataFrames` packages.
- Ensure the interpolation file path points to a valid file, typically located in the Synchrotron data folder.

# Interactive Prompts
The function interacts with the user via the terminal to gather the following information:
### 1. Base Directory
- Prompt: `Enter the base directory for simulations`
- Default: Current working directory (`pwd()`)

### 2. Simulation Selection
- Option to process all simulations or select specific ones by their indices.
- Prompt: `Do you want to process all simulations or choose specific ones?`
- Input format for specific simulations: Comma-separated indices, e.g., `1,3,5`.

### 3. Unit Conversions
- For Magnetic Field (`B`): Default unit is μG. Provide a conversion factor if in a different unit.
- For Number Density (`n`): Default unit is cm^-3. Provide a conversion factor if in a different unit.
- For Temperature (`T`): Default unit is K. Provide a conversion factor if in a different unit.

### 4. Synchrotron Data Processing
- Option to include Faraday rotation.
  - Prompt: `Do you want to include Faraday rotation in the computation of Q and U?`
- Option to perform filtering.
  - Prompt: `Do you want to perform filtering for Synchrotron data?`
  - Specify kernel size for filtering.

### 5. Lines of Sight
- Option to process all lines of sight (`x`, `y`, `z`) or choose specific ones.
- Input format for specific lines: Comma-separated values, e.g., `x,y`.

### 6. Interpolation File Path
- Prompt: `Enter the path to the interpolation file`
- Default: `Synchrotron/emissivity.dat` in the user's home directory.

# Returns
The function outputs the following computed data:
- **RM Map**: Rotation Measure map.
- **I, Q, U**: Stokes parameters for synchrotron data.
- **FDF**: Faraday Dispersion Function.

# Example
```julia
# To run the interactive workflow, simply call:
run_moose()

# Sample interaction
Enter the base directory for simulations: /path/to/simulations
Do you want to process all simulations or choose specific ones? (Enter 'all' or 'choose'): choose
Enter the indices of the simulations you want to process, separated by commas (e.g., 1,3,5): 1,2
Do you want to include Faraday rotation in the computation of Q and U? (Y/N): Y
Enter the lines of sight you want to process, separated by commas (e.g., x,y): x,z
```
"""


function print_logo()
   cols = displaysize(stdout)[2]
   rainbow = [:red, :yellow, :green, :cyan, :blue, :magenta]

   rainbow_text_line(s::AbstractString) = begin
       colors = Iterators.cycle(rainbow)
       join([Crayon(foreground=color, bold=true)(string(c)) for (c, color) in zip(collect(s), colors)], "")
   end

   if cols < 60
       colors = Iterators.cycle(rainbow)
       for c in collect("MOOSE")
           color = iterate(colors)[1]
           print(Crayon(foreground=color, bold=true)(string(c)))
           sleep(0.1)
       end
       println("│")
       println(Crayon(foreground = :light_green, bold = true)("Synchrotron Data Tool -- dev. by Jack Berat"))
   else
       println("\n" ^ 2)
       logo_text = raw"""
        ____    ____   ___      ___     ______   ________ 
       |_   \  /   _|.'   `.  .'   `. .' ____ \ |_   __  |
         |   \/   | /  .-.  \/  .-.  \| (___ \_|  | |_ \_|
         | |\  /| | | |   | || |   | | _.____`.   |  _| _ 
        _| |_\/_| |_\  `-'  /\  `-'  /| \____) | _| |__/ |
       |_____||_____|`.___.'  `.___.'  \______.'|________|
       """
       logo_lines = filter(line -> !isempty(strip(line)), split(logo_text, "\n"))
       max_len = maximum(length.(logo_lines))
       pad_left = max(0, (cols - max_len) ÷ 2)
       border = repeat("─", max_len)
       border = repeat("─", max_len)
       println(" "^pad_left * "╭" * border * "╮")
       colors = Iterators.cycle(rainbow)
       color_state = iterate(colors)
       for line in logo_lines
           print(" "^pad_left * "│")
           for c in collect(line)
               color = color_state[1]
               print(Crayon(foreground=color, bold=true)(string(c)))
               sleep(0.0001)
               color_state = iterate(colors, color_state[2])
           end
           println("│")
       end

       println(" "^pad_left * "╰" * border * "╯")

       elk_lines = [
           (:light_yellow, raw"     \\_//"),
           (:white,    raw"   __/ \" ."),
           (:white,    raw"  /__ |"),
           (:white,    raw"  || ||")
       ]
       for (color, line) in elk_lines
           for c in collect(line)
               print(Crayon(foreground=color)(string(c)))
               sleep(0.01)
           end
           println()
       end
       println(Crayon(foreground = :light_green, bold = true)("Mock Observation Of Synchrotron Emission -- dev. by Jack Berat"))
       println(Crayon(foreground = :light_red, bold = true)("Version 1.0"))
   end
end

function load_previous_config(config_path="moose_config.json")
    if isfile(config_path)
        return JSON.parsefile(config_path)
    else
        println("[Info] No existing config found at $(config_path). Starting with defaults.")
        return Dict{String, Any}()
    end
end

function save_config(config::Dict, config_path="moose_config.json")
   open(config_path, "w") do io
       write(io, JSON.json(config))
   end
end

function format_duration(elapsed)
    total_ms = Dates.value(elapsed)
    hours = total_ms ÷ 3_600_000
    minutes = (total_ms % 3_600_000) ÷ 60_000
    seconds = (total_ms % 60_000) ÷ 1_000
    milliseconds = total_ms % 1_000
    padded_hours = lpad(hours, 2, "0")
    padded_minutes = lpad(minutes, 2, "0")
    padded_seconds = lpad(seconds, 2, "0")
    padded_ms = lpad(milliseconds, 3, "0")
    return string(padded_hours, ":", padded_minutes, ":", padded_seconds, ".", padded_ms)
end

function write_summary_log(base_dir, chosen_simu, chosen_LOS, elapsed; config_path="moose_config.json", faraday="N", responseSynchrotron="N", add_noise="N", interpolation_file_path=nothing, conversionB=nothing, conversionn=nothing, conversionT=nothing, ne_option=nothing)
    log_path = joinpath(base_dir, "MOOSE_summary.log")
    open(log_path, "a") do io
        println(io)
        println(io, "MOOSE Summary Log")
        println(io, "=================")
        println(io, "Run completed at: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "Simulations processed:")
        for simu in chosen_simu
            println(io, simu)
        end
        println(io, "Lines of sight: $(join(chosen_LOS, ", "))")
        println(io, "Output directory: $base_dir")
        println(io, "Total execution time: $(format_duration(elapsed))")
        println(io, "Config file: $(config_path)")
        conversionB !== nothing && println(io, "Conversion B (to μG): $(conversionB)")
        conversionn !== nothing && println(io, "Conversion n (to cm^-3): $(conversionn)")
        conversionT !== nothing && println(io, "Conversion T (to K): $(conversionT)")
        println(io, "Faraday rotation: $(faraday)")
        println(io, "Synchrotron filtering: $(responseSynchrotron)")
        println(io, "Noise added: $(add_noise)")
        interpolation_file_path !== nothing && println(io, "Interpolation file: $(interpolation_file_path)")
        ne_option !== nothing && println(io, "Electron density option: $(ne_option)")
    end
end


struct RunConfig
    base_dir::String
    simulations::Vector{String}
    chosen_LOS::Vector{String}
    conversionB::Float64
    conversionn::Float64
    conversionT::Float64
    faraday_rotation::String
    phimin::Float64
    phimax::Float64
    dphi::Float64
    responseSynchrotron::String
    kernel_size_synchrotron::Union{Nothing, Int}
    add_noise::String
    SNR_nu::Union{Nothing, Float64}
    interpolation_file_path::String
    ne_option::String
    IonizationFraction::Union{Nothing, Float64}
    nustart::Float64
    nuend::Float64
    dnu::Float64
    BoxLength_pc::Float64
    BoxLength_pix::Int
    config_path::String

    function RunConfig(
        base_dir,
        simulations,
        chosen_LOS,
        conversionB,
        conversionn,
        conversionT,
        faraday_rotation,
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
    )
        faraday_flag = uppercase(faraday_rotation)
        response_flag = uppercase(responseSynchrotron)
        noise_flag = uppercase(add_noise)
        sim_paths = map(String, simulations)
        los_list = map(String, chosen_LOS)

        kernel_size_synchrotron = kernel_size_synchrotron === nothing ? nothing : Int(kernel_size_synchrotron)
        SNR_nu = SNR_nu === nothing ? nothing : Float64(SNR_nu)
        IonizationFraction = IonizationFraction === nothing ? nothing : Float64(IonizationFraction)

        return new(
            String(base_dir),
            sim_paths,
            los_list,
            Float64(conversionB),
            Float64(conversionn),
            Float64(conversionT),
            faraday_flag,
            Float64(phimin),
            Float64(phimax),
            Float64(dphi),
            response_flag,
            kernel_size_synchrotron,
            noise_flag,
            SNR_nu,
            String(interpolation_file_path),
            String(ne_option),
            IonizationFraction,
            Float64(nustart),
            Float64(nuend),
            Float64(dnu),
            Float64(BoxLength_pc),
            Int(BoxLength_pix),
            String(config_path),
        )
    end
end

function config_dict_from_struct(cfg::RunConfig)
    return Dict(
        "base_dir" => cfg.base_dir,
        "chosen_simu" => cfg.simulations,
        "chosen_LOS" => cfg.chosen_LOS,
        "conversionB" => cfg.conversionB,
        "conversionn" => cfg.conversionn,
        "conversionT" => cfg.conversionT,
        "FaradayRotation" => cfg.faraday_rotation,
        "phimin" => cfg.phimin,
        "phimax" => cfg.phimax,
        "dphi" => cfg.dphi,
        "responseSynchrotron" => cfg.responseSynchrotron,
        "kernel_size_synchrotron" => cfg.kernel_size_synchrotron,
        "add_noise" => cfg.add_noise,
        "SNR_nu" => cfg.SNR_nu,
        "interpolation_file_path" => cfg.interpolation_file_path,
        "ne_option" => cfg.ne_option,
        "IonizationFraction" => cfg.IonizationFraction,
        "nustart" => cfg.nustart,
        "nuend" => cfg.nuend,
        "dnu" => cfg.dnu,
        "BoxLength_pc" => cfg.BoxLength_pc,
        "BoxLength_pix" => cfg.BoxLength_pix,
        "config_path" => cfg.config_path,
    )
end

function run_moose_processing(cfg::RunConfig; quiet::Bool = false, persisted_config::Union{Nothing, Dict} = nothing)
    nuArray = range(start = cfg.nustart, stop = cfg.nuend, step = cfg.dnu)
    PhiArray = cfg.faraday_rotation == "Y" ? range(start = cfg.phimin, stop = cfg.phimax, step = cfg.dphi) : nothing
    PixelLength_pc, PixelLength_cm, DistanceArray = los_pixel_scale(cfg.BoxLength_pc, cfg.BoxLength_pix)

    isfile(cfg.interpolation_file_path) || error("The interpolation file $(cfg.interpolation_file_path) was not found.")
    df = CSV.File(cfg.interpolation_file_path) |> DataFrame

    if !quiet
        print_logo()
    end

    start_time = now()

    if cfg.ne_option == "3"
        missing_cubes = [simu for simu in cfg.simulations if !isfile(joinpath(simu, "densityHp.fits"))]
        !isempty(missing_cubes) && error("Electron density cube 'densityHp.fits' is missing for: $(join(missing_cubes, ", ")).")
    end

    for (i, simu) in enumerate(cfg.simulations)
        println("Processing Simulation: $(simu)")

        for LOS in cfg.chosen_LOS
            println(Crayon(foreground=:yellow, bold=true)("→ Processing LOS: $(LOS)"))

            if cfg.ne_option == "1"
                zeta, Geff, omegaPAH, XC = WolfireConstants()
                ProcessSynchrotron(simu, LOS, cfg.faraday_rotation, cfg.responseSynchrotron, df, cfg.add_noise, cfg.SNR_nu,
                    cfg.kernel_size_synchrotron, zeta, Geff, omegaPAH, XC, nuArray, PhiArray, PixelLength_pc, PixelLength_cm,
                    cfg.BoxLength_pc, DistanceArray, cfg.conversionn, cfg.conversionT, cfg.conversionB)
            elseif cfg.ne_option == "2"
                ion_fraction = cfg.IonizationFraction === nothing ? 0.01 : cfg.IonizationFraction
                ProcessSynchrotron(simu, LOS, cfg.faraday_rotation, cfg.responseSynchrotron, df, cfg.add_noise, cfg.SNR_nu,
                    cfg.kernel_size_synchrotron, ion_fraction, nuArray, PhiArray, PixelLength_pc, PixelLength_cm,
                    cfg.BoxLength_pc, DistanceArray, cfg.conversionn, cfg.conversionT, cfg.conversionB)
            else
                ProcessSynchrotron(simu, LOS, cfg.faraday_rotation, cfg.responseSynchrotron, df, cfg.add_noise, cfg.SNR_nu,
                    cfg.kernel_size_synchrotron, nuArray, PhiArray, PixelLength_pc, PixelLength_cm,
                    cfg.BoxLength_pc, DistanceArray, cfg.conversionn, cfg.conversionT, cfg.conversionB)
            end
        end

        if length(cfg.simulations) > 1
            println("Finished processing all chosen LOS for simulation: $simu")
            print_progress(i, length(cfg.simulations))
        end
    end

    println("Finished processing all simulations.")

    elapsed = now() - start_time
    println(Crayon(foreground=:green, bold=true)("Summary:"))
    println(Crayon(foreground=:green)("Simulations processed: $(join(cfg.simulations, ", "))"))
    println(Crayon(foreground=:green)("Lines of sight: $(join(cfg.chosen_LOS, ", "))"))
    println(Crayon(foreground=:green)("Output directory: $(cfg.base_dir)"))
    println(Crayon(foreground=:green)("Total execution time: $(format_duration(elapsed))"))

    config_to_save = persisted_config === nothing ? config_dict_from_struct(cfg) : persisted_config
    save_config(config_to_save, cfg.config_path)
    write_summary_log(cfg.base_dir, map(basename, cfg.simulations), cfg.chosen_LOS, elapsed; config_path=cfg.config_path,
        faraday=cfg.faraday_rotation, responseSynchrotron=cfg.responseSynchrotron, add_noise=cfg.add_noise,
        interpolation_file_path=cfg.interpolation_file_path, conversionB=cfg.conversionB, conversionn=cfg.conversionn,
        conversionT=cfg.conversionT, ne_option=cfg.ne_option)
end


function run_moose(; quiet::Bool = false, reset_config::Bool = true, help::Bool = false)
   if help
       println("""
MOOSE v1.0 — Mock Observation Of Synchrotron Emission

Usage:
 run_moose(; quiet=false, reset_config=true, help=false)

Options:
 --quiet             Disable the rainbow logo at startup.
 --reset-config      Ignore previous config and prompt user again.
 --help              Show this help message and exit.

Description:
 Author: Jack Berat
 MOOSE is an interactive tool to process synchrotron mock observation data.

 Non-interactive CLI:
  For scripted or assistive use, run `julia src/MOOSE_cli.jl --help` to see the
  full set of command-line flags (e.g., --config, --base-dir, --simu, --los,
  --quiet) for configuring MOOSE without prompts.

 It computes the synchrotron Stokes parameters Q and U for a set of simulation outputs,
 optionally applying Faraday rotation, noise and primary beam filtering.

 The tool interacts with the user to define simulation directories, select simulations and
 lines of sight, and configure physical unit conversions and data processing options.

 Previous configuration is loaded from `moose_config.json`, unless reset.
Outputs are saved to the base simulation directory, with logs in `MOOSE_summary.log`.

Flow:
 ┌─────────────────────────────┐
 │ Interactive config (or JSON)│
 └────────────┬────────────────┘
              │
              ▼
   Read simulation + unit setup
              │
              ▼
   Select LOS + process Q/U/I
              │
              ▼
   Apply noise + Gaussian Primary Beam
              │
              ▼
   Process RM synthesis on Q and U
              │
              ▼
   Save FITS in the simulation file + logs

""")
       return
   end
    cfg, persisted_config = run_moose_interactive(; quiet = quiet, reset_config = reset_config)
    run_moose_processing(cfg; quiet = true, persisted_config = persisted_config)
end

function run_moose_interactive(; quiet::Bool = false, reset_config::Bool = true)
    if !quiet
        print_logo()
    end

    default_config_path = joinpath(pwd(), "moose_config.json")
    if reset_config
        println("[Info] Previous configuration ignored (reset_config=true)")
        config = Dict{String, Any}()
        config_path = default_config_path
    else
        config_path = ask_user("Enter the path to the configuration file to load", default_config_path)
        config = load_previous_config(config_path)
    end

    base_dir = ""
    while true
        candidate_dir = ask_user("Enter the base directory for simulations", get(config, "base_dir", pwd()))
        if isdir(candidate_dir)
            base_dir = candidate_dir
            break
        else
            println("[Error] Base directory $(candidate_dir) does not exist. Please provide a valid folder.")
        end
    end
    config["base_dir"] = base_dir

    simu_list = get_simulation_list(base_dir)
    isempty(simu_list) && error("No simulations containing FITS files were found in $(base_dir).")
    display_simulations(simu_list)

    simu_prompt = "Enter 'all' to process all simulations or provide comma-separated indices (e.g., 1,3,5):"
    simu_choice = ask_user(simu_prompt, get(config, "simu_choice", "all"))
    chosen_simu = begin
        parsed_indices = Int[]
        while true
            if uppercase(strip(simu_choice)) == "ALL"
                empty!(parsed_indices)
                append!(parsed_indices, 1:length(simu_list))
                break
            end

            raw_indices = split(simu_choice, ",")
            empty!(parsed_indices)
            for raw_idx in raw_indices
                candidate = strip(raw_idx)
                isempty(candidate) && continue
                parsed = tryparse(Int, candidate)
                if parsed === nothing
                    println("[Warning] Ignoring invalid simulation index: $(candidate). Enter comma-separated integers like 1,3,5.")
                elseif parsed < 1 || parsed > length(simu_list)
                    println("[Warning] Simulation index $(parsed) is out of range (1-$(length(simu_list))).")
                else
                    push!(parsed_indices, parsed)
                end
            end

            if !isempty(parsed_indices)
                break
            end

            println("[Error] No valid simulation indices provided. Use 'all' or comma-separated integers like 1,3,5 and try again.")
            simu_choice = ask_user(simu_prompt, simu_choice)
        end

        config["simu_choice"] = simu_choice
        config["simu_indices"] = uppercase(strip(simu_choice)) == "ALL" ? "" : join(parsed_indices, ",")
        map(i -> simu_list[i], parsed_indices)
    end
    config["chosen_simu"] = chosen_simu

    conversionB = ask_user("Enter the conversion factor for magnetic field B to μG (microGauss):", Float64(get(config, "conversionB", 1.0)))
    conversionn = ask_user("Enter the conversion factor for number density n to cm^-3:", Float64(get(config, "conversionn", 1.0)))
    conversionT = ask_user("Enter the conversion factor for temperature T to K:", Float64(get(config, "conversionT", 1.0)))
    config["conversionB"] = conversionB
    config["conversionn"] = conversionn
    config["conversionT"] = conversionT

    BoxLength_pc = ask_user("Side of the Box size (pc), please give a Float", Float64(get(config, "BoxLength_pc", 50.0)))
    BoxLength_pix = ask_user("Number of pixels along the line of sight", Int(get(config, "BoxLength_pix", 256)))
    config["BoxLength_pc"] = BoxLength_pc
    config["BoxLength_pix"] = BoxLength_pix

    nustart = ask_user("Frequency range start (MHz)", Float64(get(config, "nustart", 120)))
    nuend = ask_user("Frequency range end (MHz)", Float64(get(config, "nuend", 167)))
    dnu = ask_user("Frequency resolution (MHz)", Float64(get(config, "dnu", 0.098)))
    config["nustart"] = nustart
    config["nuend"] = nuend
    config["dnu"] = dnu

    FaradayRotation = ask_user("Do you want to include Faraday rotation in the computation of Q and U? (Y/N)", get(config,"FaradayRotation", "N"))
    faraday_flag = uppercase(FaradayRotation)
    config["FaradayRotation"] = FaradayRotation
    phimin = get(config, "phimin", -10.0)
    phimax = get(config, "phimax", 10.0)
    dphi = get(config, "dphi", 0.25)
    if faraday_flag == "Y"
        phimin = ask_user("Faraday depth range start (rad/m^2)", Float64(phimin))
        phimax = ask_user("Faraday depth range end (rad/m^2)", Float64(phimax))
        dphi = ask_user("Faraday depth resolution (rad/m^2)", Float64(dphi))
    end
    config["phimin"] = phimin
    config["phimax"] = phimax
    config["dphi"] = dphi

    responseSynchrotron = ask_user("Do you want to perform filtering (primary beam) for Synchrotron data? (Y/N)", get(config, "responseSynchrotron", "N"))
    kernel_size_synchrotron = uppercase(responseSynchrotron) == "Y" ? ask_user("What kernel size (in pix) do you want for Synchrotron filtering?", get(config, "kernel_size_synchrotron", 11)) : nothing
    config["responseSynchrotron"] = responseSynchrotron
    config["kernel_size_synchrotron"] = kernel_size_synchrotron

    add_noise = ask_user("Do you want to add noise to Q and U? (Y/N)", get(config, "add_noise", "N"))
    config["add_noise"] = add_noise

    SNR_nu = uppercase(add_noise) == "Y" ? ask_user("Enter the desired SNR in the frequency space:", get(config, "SNR_nu", 0.9)) : nothing
    config["SNR_nu"] = SNR_nu

    list_LOS = ["x", "y", "z"]
    los_prompt = "Enter 'all' to process all lines of sight (x, y, z) or specify comma-separated ones (e.g., x,y):"

    function parse_los_choice(los_choice)
        valid_los = String[]
        if uppercase(strip(los_choice)) == "ALL"
            append!(valid_los, list_LOS)
            return valid_los
        end

        for los in split(los_choice, ",")
            candidate = lowercase(strip(los))
            isempty(candidate) && continue
            if candidate in list_LOS
                push!(valid_los, candidate)
            else
                println("[Warning] Ignoring invalid line of sight: $(candidate). Valid options are x, y, z.")
            end
        end

        return valid_los
    end

    los_choice = get(config, "chosen_LOS_input", "all")
    valid_los = String[]
    while isempty(valid_los)
        los_choice = ask_user(
            los_prompt,
            los_choice;
            validate = choice -> !isempty(parse_los_choice(choice)),
            error_message = "[Error] No valid lines of sight provided. Use 'all' or comma-separated values like x,y and try again.",
        )
        valid_los = parse_los_choice(los_choice)
    end

    config["chosen_LOS_input"] = uppercase(strip(los_choice)) == "ALL" ? "all" : join(valid_los, ",")
    chosen_LOS = valid_los
    config["chosen_LOS"] = chosen_LOS

    interpolation_default = get(config, "interpolation_file_path", joinpath(homedir(), "emissivity.dat"))
    interpolation_file_path = interpolation_default
    while true
        interpolation_file_path = ask_user("Enter the path to the interpolation file", interpolation_file_path)
        if isfile(interpolation_file_path)
            break
        else
            println("[Error] The interpolation file $(interpolation_file_path) was not found. Please provide a valid path.")
        end
    end
    config["interpolation_file_path"] = interpolation_file_path

    ne_option = ""
    while true
        ne_option = ask_user("Choose electron density prescription: (1) Wolfire et al. 2003, (2) Proportional to nH, (3) Provide ne cube", get(config, "ne_option", "1"))
        ne_option in ("1", "2", "3") && break
        println("[Warning] Please choose 1, 2, or 3 for the electron density prescription.")
    end
    config["ne_option"] = ne_option

    IonizationFraction = get(config, "IonizationFraction", nothing)
    if ne_option == "2"
        IonizationFraction = ask_user("Enter the ionization fraction for the alternative prescription:", get(config, "IonizationFraction", 0.01))
        config["IonizationFraction"] = IonizationFraction
    elseif ne_option == "3"
        missing_cubes = [simu for simu in chosen_simu if !isfile(joinpath(simu, "densityHp.fits"))]
        if !isempty(missing_cubes)
            error("Electron density cube 'densityHp.fits' is missing for: $(join(missing_cubes, ", ")).")
        end
    end

    save_path_default = get(config, "config_path", joinpath(base_dir, "moose_config.json"))
    config_path = ask_user("Enter the path where the configuration should be saved", save_path_default)
    config["config_path"] = config_path

    cfg = RunConfig(
        base_dir,
        chosen_simu,
        chosen_LOS,
        conversionB,
        conversionn,
        conversionT,
        faraday_flag,
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
    )

    return cfg, config
end
