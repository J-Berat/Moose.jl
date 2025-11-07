"""
    MOOSE()

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
# To run the MOOSE function, simply call it:
MOOSE()

# Sample interaction
Enter the base directory for simulations: /path/to/simulations
Do you want to process all simulations or choose specific ones? (Enter 'all' or 'choose'): choose
Enter the indices of the simulations you want to process, separated by commas (e.g., 1,3,5): 1,2
Do you want to include Faraday rotation in the computation of Q and U? (Y/N): Y
Enter the lines of sight you want to process, separated by commas (e.g., x,y): x,z
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
       logo_lines = split(logo_text, "\n")
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
    isfile(config_path) ? JSON.parsefile(config_path) : Dict{String, Any}()
end

function save_config(config::Dict, config_path="moose_config.json")
   open(config_path, "w") do io
       JSON.print(io, config)
   end
end

function write_summary_log(base_dir, chosen_simu, chosen_LOS, elapsed)
    log_path = joinpath(base_dir, "MOOSE_summary.log")
    open(log_path, "w") do io
        println(io, "MOOSE Summary Log")
        println(io, "=================")
        println(io, "Simulations processed:")
        for simu in chosen_simu
            println(io, simu)
        end
        println(io, "Lines of sight: $(join(chosen_LOS, ", "))")
        println(io, "Output directory: $base_dir")
        println(io, "Execution completed at: $(now())")
        println(io, "Total execution time: $(Dates.value(elapsed) ÷ 1_000) seconds")
        println(io, "Using config file: moose_config.json")
    end
end


function MOOSE(; quiet::Bool = false, reset_config::Bool = true, help::Bool = false)
   if help
       println("""
MOOSE v1.0 — Mock Observation Of Synchrotron Emission

Usage:
 MOOSE(; quiet=false, reset_config=true, help=false)

Options:
 --quiet             Disable the rainbow logo at startup.
 --reset-config      Ignore previous config and prompt user again.
 --help              Show this help message and exit.

Description:
 Author: Jack Berat
 MOOSE is an interactive tool to process synchrotron mock observation data.

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

   if !quiet
       print_logo()
   end

   if reset_config
       println("[Info] Previous configuration ignored (reset_config=true)")
   end
   config = reset_config ? Dict{String, Any}() : load_previous_config()
   start_time = now()

   base_dir = get(config, "base_dir", ask_user("Enter the base directory for simulations", pwd()))
   config["base_dir"] = base_dir

   simu_list = get_simulation_list(base_dir)
   display_simulations(simu_list)

   simu_choice = ask_user("Do you want to process all simulations or choose specific ones? (Enter 'all' or 'choose')", "all")
   chosen_simu = if uppercase(simu_choice) == "CHOOSE"
       indices = split(ask_user("Enter the indices of the simulations you want to process, separated by commas (e.g., 1,3,5): ", ""), ",")
       map(i -> simu_list[parse(Int, i)], indices)
   else
       simu_list
   end
   config["chosen_simu"] = chosen_simu

   conversionB = ask_user("Enter the conversion factor for magnetic field B to μG (microGauss):", get(config, "conversionB", 1.0))
   conversionn = ask_user("Enter the conversion factor for number density n to cm^-3:", get(config, "conversionn", 1.0))
   conversionT = ask_user("Enter the conversion factor for temperature T to K:", get(config, "conversionT", 1.0))
   config["conversionB"] = conversionB
   config["conversionn"] = conversionn
   config["conversionT"] = conversionT

   PixelLength_pc, PixelLength_cm, BoxLength_pc, DistanceArray = DistanceParameters()
   nuArray = FrequencyParameters()
   FaradayRotation = ask_user("Do you want to include Faraday rotation in the computation of Q and U? (Y/N)", get(config,"FaradayRotation", "N"))
   config["FaradayRotation"] = FaradayRotation
   PhiArray = uppercase(FaradayRotation) == "Y" ? FaradayParameters() : nothing

   responseSynchrotron = ask_user("Do you want to perform filtering (primary beam) for Synchrotron data? (Y/N)", get(config, "responseSynchrotron", "N"))
   kernel_size_synchrotron = uppercase(responseSynchrotron) == "Y" ? ask_user("What kernel size (in pix) do you want for Synchrotron filtering?", get(config, "kernel_size_synchrotron", 11)) : nothing
   config["responseSynchrotron"] = responseSynchrotron
   config["kernel_size_synchrotron"] = kernel_size_synchrotron

   add_noise = ask_user("Do you want to add noise to Q and U? (Y/N)", get(config, "add_noise", "N"))
   config["add_noise"] = add_noise

   SNR_nu = uppercase(add_noise) == "Y" ? ask_user("Enter the desired SNR in the frequency space:", get(config, "SNR_nu", 0.9)) : nothing
   config["SNR_nu"] = SNR_nu

   list_LOS = ["x", "y", "z"]
   LOS_choice = ask_user("Do you want to process all lines of sight (x, y, z), or choose specific ones? (Enter 'all' or 'choose')", get(config, "LOS_choice", "All"))
   chosen_LOS = uppercase(LOS_choice) == "CHOOSE" ? split(ask_user("Enter the lines of sight you want to process, separated by commas (e.g., x,y): ", ""), ",") : list_LOS
   config["chosen_LOS"] = chosen_LOS

   interpolation_file_path = ask_user("Enter the path to the interpolation file", get(config, "interpolation_file_path", joinpath(homedir(), "emissivity.dat")))
   config["interpolation_file_path"] = interpolation_file_path
   df = CSV.File(interpolation_file_path) |> DataFrame

   ne_option = ask_user("Choose electron density prescription: (1) Wolfire et al. 2003, (2) Proportional to nH, (3) Provide ne cube", get(config, "ne_option", "1"))
   config["ne_option"] = ne_option

   if ne_option == "1"
       zeta, Geff, omegaPAH, XC = WolfireConstants()
   elseif ne_option == "2"
       IonizationFraction = ask_user("Enter the ionization fraction for the alternative prescription:", get(config, "IonizationFraction", 0.01))
       config["IonizationFraction"] = IonizationFraction
   else
       println("The electron density cube must be present in the simulation directory and named 'densityHp.fits'")
   end

   total_simu = length(chosen_simu)
   for (i, simu) in enumerate(chosen_simu)

       println("Processing Simulation: $(simu)")
       
       for LOS in chosen_LOS
               println(Crayon(foreground=:yellow, bold=true)("→ Processing LOS: $(LOS)"))

           if ne_option == "1"
               ProcessSynchrotron(simu, LOS, FaradayRotation, responseSynchrotron, df, add_noise, SNR_nu,
                   kernel_size_synchrotron, zeta, Geff, omegaPAH, XC, nuArray, PhiArray, PixelLength_pc, PixelLength_cm, 
                   BoxLength_pc, DistanceArray, conversionn, conversionT, conversionB)
           elseif ne_option == "2"
               ProcessSynchrotron(simu, LOS, FaradayRotation, responseSynchrotron, df, add_noise, SNR_nu,
                   kernel_size_synchrotron, IonizationFraction, nuArray, PhiArray, PixelLength_pc, PixelLength_cm,
                   BoxLength_pc, DistanceArray, conversionn, conversionT, conversionB)
           else
               ProcessSynchrotron(simu, LOS, FaradayRotation, responseSynchrotron, df, add_noise, SNR_nu, kernel_size_synchrotron, 
                   nuArray, PhiArray, PixelLength_pc, PixelLength_cm, 
                   BoxLength_pc, DistanceArray, conversionn, conversionT, conversionB)
           end
       end
       if length(chosen_simu) > 1
           println("Finished processing all chosen LOS for simulation: $simu")
           print_progress(i, total_simu) 
       end
   end

   println("Finished processing all simulations.")

   elapsed = now() - start_time
   println(Crayon(foreground=:green, bold=true)("Summary:"))
   println(Crayon(foreground=:green)("Simulations processed: $(join(chosen_simu, ", "))"))
   println(Crayon(foreground=:green)("Lines of sight: $(join(chosen_LOS, ", "))"))
   println(Crayon(foreground=:green)("Output directory: $base_dir"))
   println(Crayon(foreground=:green)("Total execution time: $(Dates.value(elapsed) ÷ 1_000) seconds"))
   reset_config == true ? println(Crayon(foreground=:green)("Using config file: moose_config.json")) : println(Crayon(foreground=:green)("No config file used."))
       
    save_config(config, base_dir * "/moose_config.json")
    write_summary_log(base_dir, chosen_simu, chosen_LOS, elapsed)
end