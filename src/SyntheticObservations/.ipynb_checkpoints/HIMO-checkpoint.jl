"""
    HIMO()

Interactive tool to process HI data for a set of simulations.

# Description
The `HIMO` function is an interactive tool that guides the user through processing HI data for a set of simulations. It prompts the user for various inputs, including the base directory of simulations, specific simulations to process, lines of sight, and parameters for HI data processing. The function then processes the chosen simulations and lines of sight, performing calculations and saving results as required.

# Interactive Prompts
1. **Base Directory**: The user is prompted to enter the base directory for simulations.
2. **Simulation Selection**: The user can choose to process all simulations or select specific ones.
3. **HI Data Processing**: The user is prompted to decide whether to perform filtering for HI data and to provide temperatures for CNM and WNM.
4. **Lines of Sight**: The user can choose to process all lines of sight or select specific ones.

# Returns
- `Nothing`: The function does not return any value. It performs data processing and prints progress and results to the console.

# Example
```julia
# To run the HIMO function, simply call it:
HIMO()
"""
function HIMO()
    # Ask user for the base directory of simulations
    base_dir = ask_user("Enter the base directory for simulations", pwd())
    simu_list = get_simulation_list(base_dir)
    
    # Display available simulations
    display_simulations(simu_list)

    # Ask user to select specific simulations or all
    simu_choice = ask_user("Do you want to process all simulations or choose specific ones? (Enter 'all' or 'choose')", "all")
    chosen_simu = if simu_choice == "choose"
        indices = split(ask_user("Enter the indices of the simulations you want to process, separated by commas (e.g., 1,3,5): ", ""), ",")
        map(i -> simu_list[parse(Int, i)], indices)
    else
        simu_list
    end

    # Check if the units for each parameter are correct and ask for conversion if necessary

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

    # Ask user for HI data processing
    kernel_size_hi = 5
    responseHI = ask_user("Do you want to perform filtering for HI data?", "N")
    if uppercase(responseHI) == "Y"
        kernel_size_hi = ask_user("What kernel size do you want for HI filtering?", 5)
    end
    TCNM = ask_user("What is your CNM temperature?", 200)
    TWNM = ask_user("What is your WNM temperature?", 2000)
    velArray = VelocityParameters()
    PixelLength_pc, PixelLength_cm, BoxLength_pc, DistanceArray = DistanceParameters()
    list_LOS = ["x", "y", "z"]

    # Ask user for lines of sight to process
    LOS_choice = ask_user("Do you want to process all lines of sight (x, y, z), or choose specific ones? (Enter 'all' or 'choose')", "All")
    chosen_LOS = uppercase(LOS_choice) == "CHOOSE" ? split(ask_user("Enter the lines of sight you want to process, separated by commas (e.g., x,y): ", ""), ",") : list_LOS

    # Loop through each chosen simulation and process the chosen lines of sight
    total_simu = length(chosen_simu)
    for (i, simu) in enumerate(chosen_simu)
        println("------------------------------------------------------------------------------------------------")
        println("Processing simulation: $simu")
        for LOS in chosen_LOS
            println("Processing LOS: $LOS")
            ProcessHI(simu, LOS, responseHI, kernel_size_hi, TCNM, TWNM, PixelLength_cm, BoxLength_pc, velArray, conversionn, conversionT, conversionV, 1)
        end
        if length(chosen_simu) > 1
            println("Finished processing all chosen LOS for simulation: $simu")
            print_progress(i, total_simu)  # Update progress bar
        end
    end
    println("\nFinished processing all simulations.")
end