"""
Helpers to locate simulations and report progress while iterating over them.
"""

function contains_fits_files(dir)
    return any(f -> endswith(f, ".fits"), readdir(dir))
end

function get_simulation_list(base_dir)
    simulation_dirs = String[]
    dirs_to_check = [base_dir]
    while !isempty(dirs_to_check)
        current_dir = popfirst!(dirs_to_check)
        if contains_fits_files(current_dir)
            push!(simulation_dirs, current_dir)
        else
            subdirs = [joinpath(current_dir, d) for d in readdir(current_dir) if isdir(joinpath(current_dir, d))]
            append!(dirs_to_check, subdirs)
        end
    end
    return sort!(simulation_dirs)
end

function display_simulations(simu_list)
    println("Available simulations:")
    for (i, simu) in enumerate(simu_list)
        println("[$i] $simu")
    end
end
