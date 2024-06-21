"""
    MOOSE()

Interactive tool to process synchrotron data for a set of simulations.

# Description
The `MOOSE` function is an interactive tool that guides the user through processing synchrotron data for a set of simulations. It prompts the user for various inputs, including the base directory of simulations, specific simulations to process, lines of sight, and parameters for synchrotron data processing. The function then processes the chosen simulations and lines of sight, performing calculations and saving results as required.

# Interactive Prompts
1. **Base Directory**: The user is prompted to enter the base directory for simulations.
2. **Simulation Selection**: The user can choose to process all simulations or select specific ones.
3. **Synchrotron Data Processing**: The user is prompted to decide whether to include Faraday rotation and perform filtering for synchrotron data.
4. **Lines of Sight**: The user can choose to process all lines of sight or select specific ones.
5. **Interpolation File Path**: The user is prompted to provide the path to the interpolation file.

# Returns
- `Nothing`: The function does not return any value. It performs data processing and prints progress and results to the console.

# Example
```julia
# To run the MOOSE function, simply call it:
MOOSE()
"""
function MOOSE()
    # Ask user for the base directory of simulations
    base_dir = ask_user("Enter the base directory for simulations", pwd())
    simu_list = get_simulation_list(base_dir)
    
    # Display available simulations
    display_simulations(simu_list)

    # Ask user to select specific simulations or all
    simu_choice = ask_user("Do you want to process all simulations or choose specific ones? (Enter 'all' or 'choose')", "all")
    chosen_simu = if uppercase(simu_choice) == "CHOOSE"
        indices = split(ask_user("Enter the indices of the simulations you want to process, separated by commas (e.g., 1,3,5): ", ""), ",")
        map(i -> simu_list[parse(Int, i)], indices)
    else
        simu_list
    end

    # Check if the units for each parameter are correct and ask for conversion if necessary
    unit_response_B = ask_user("Is the unit of magnetic field B in μG (microGauss)? (Y/N)", "N")
    if uppercase(unit_response_B) == "N"
        conversionB = ask_user("Enter the conversion factor for magnetic field B to μG (microGauss):", 1e3)
    else
        conversionB = 1
    end
    
    unit_response_n = ask_user("Is the unit of number density n in cm^-3? (Y/N)", "N")
    if uppercase(unit_response_n) == "N"
        conversionn = ask_user("Enter the conversion factor for number density n to cm^-3:", 1.0)
    else
        conversionn = 1
    end
    
    unit_response_T = ask_user("Is the unit of temperature T in K? (Y/N)", "N")
    if uppercase(unit_response_T) == "N"
        conversionT = ask_user("Enter the conversion factor for temperature T to K:", 1.0)
    else
        conversionT = 1
    end
    
    unit_response_V = ask_user("Is the unit of velocity V in km/s? (Y/N)", "N")
    if uppercase(unit_response_V) == "N"
        conversionV = ask_user("Enter the conversion factor for velocity V to km/s:", 1.0)
    else
        conversionV = 1
    end
    
    # Ask user for Synchrotron data processing
    FaradayRotation = ""
    responseSynchrotron = ""
    kernel_size_synchrotron = 5

    zeta, Geff, omegaPAH, XC = WolfireConstants()
    nuArray, PhiArray, PixelLength_pc, PixelLength_cm, BoxLength_pc = SynchrotronInstrumentalParameters()
    DistanceArray = range(start = 0, stop = BoxLength_pc, step = PixelLength_pc)

    FaradayRotation = ask_user("Do you want to include Faraday rotation in the computation of Q and U?", "N")
    responseSynchrotron = ask_user("Do you want to perform filtering for Synchrotron data?", "N")
    if uppercase(responseSynchrotron) == "Y"
        kernel_size_synchrotron = ask_user("What kernel size (in pix) do you want for Synchrotron filtering?", 5)
    end

    # Ask user for lines of sight to process
    list_LOS = ["x", "y", "z"]
    LOS_choice = ask_user("Do you want to process all lines of sight (x, y, z), or choose specific ones? (Enter 'all' or 'choose')", "All")
    chosen_LOS = uppercase(LOS_choice) == "CHOOSE" ? split(ask_user("Enter the lines of sight you want to process, separated by commas (e.g., x,y): ", ""), ",") : list_LOS

    # Ask user for the interpolation file path
    interpolation_file_path = ask_user("Enter the path to the interpolation file", joinpath(homedir(), "Synchrotron/emissivity.dat"))
    df = CSV.File(interpolation_file_path) |> DataFrame

    # Loop through each chosen simulation and process the chosen lines of sight
    total_simu = length(chosen_simu)
    for (i, simu) in enumerate(chosen_simu)
        println("------------------------------------------------------------------------------------------------")
        println("Processing simulation: $simu")
        for LOS in chosen_LOS
            println("Processing LOS: $LOS")
            ProcessSynchrotron(simu, LOS, FaradayRotation, responseSynchrotron, df, kernel_size_synchrotron, zeta, Geff, omegaPAH, XC, nuArray, PhiArray, PixelLength_pc, PixelLength_cm, BoxLength_pc, DistanceArray, conversionn, conversionT, conversionV, conversionB)
        end
        if length(chosen_simu) > 1
            println("Finished processing all chosen LOS for simulation: $simu")
            print_progress(i, total_simu)  # Update progress bar
        end
    end
    println("\nFinished processing all simulations.")
end