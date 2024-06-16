"""
    max(cube::AbstractArray) -> AbstractArray

Calculate the maximum values along the third dimension of a 3D array.

# Arguments
- `cube::AbstractArray`: 3D array where the maximum values will be calculated along the third dimension.

# Returns
- `AbstractArray`: 2D array representing the maximum values along the third dimension for each (x, y) position in `cube`.

# Description
This function calculates the maximum values of the input 3D array `cube` along the third dimension and returns a 2D array. Each element of the returned array represents the maximum value over the third dimension for each corresponding (x, y) position in `cube`.

# Example
```julia
# Example input array
cube = rand(100, 100, 50)  # Example data cube with dimensions 100x100x50

# Function call
max_values = max(cube)
"""
max(cube::AbstractArray) = dropdims(maximum(cube, dims=3), dims=3)

"""
    intLOS(cube::AbstractArray, PixelLength_cm::Float64) -> AbstractArray

Calculate the integrated line-of-sight (LOS) values for a 3D array.

# Arguments
- `cube::AbstractArray`: 3D array representing the data to be integrated along the third dimension.
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `AbstractArray`: 2D array representing the integrated values along the third dimension of `cube`.

# Description
This function calculates the integrated line-of-sight (LOS) values for the input 3D array `cube` by summing along the third dimension after multiplying each element by `PixelLength_cm`. The result is a 2D array where each element represents the integrated value along the third dimension for the corresponding (x, y) position.

# Example
```julia
# Example input array
cube = rand(100, 100, 50)  # Example data cube with dimensions 100x100x50
PixelLength_cm = 1.0

# Function call
integrated_los = intLOS(cube, PixelLength_cm)
"""
intLOS(cube::AbstractArray, PixelLength_cm::Float64) = dropdims(sum(cube .* PixelLength_cm, dims=3), dims=3)

"""
    MeanSpectrum(cube::AbstractArray) -> AbstractArray

Calculate the mean spectrum along the first and second dimensions of a 3D array.

# Arguments
- `cube::AbstractArray`: 3D array where the mean will be calculated along the first and second dimensions.

# Returns
- `AbstractArray`: 1D array representing the mean values along the first and second dimensions for each element in the third dimension of `cube`.

# Description
This function calculates the mean of the input 3D array `cube` along the first and second dimensions and returns a 1D array. Each element of the returned array represents the mean value over the first and second dimensions for each corresponding element in the third dimension of `cube`.

# Example
```julia
# Example input array
cube = rand(100, 100, 50)  # Example data cube with dimensions 100x100x50

# Function call
mean_spectrum = MeanSpectrum(cube)
"""
MeanSpectrum(cube::AbstractArray) = dropdims(mean(cube, dims=(1, 2)), dims=(1, 2))

"""
    MaxIndicesMap(cube::AbstractArray, ValueArray::AbstractArray) -> AbstractArray

Generate a map of the values from `ValueArray` corresponding to the maximum indices in the third dimension of `cube`.

# Arguments
- `cube::AbstractArray`: 3D array where the maximum value indices will be calculated along the third dimension.
- `ValueArray::AbstractArray`: 1D array of values corresponding to the third dimension of `cube`.

# Returns
- `AbstractArray`: 2D array where each element represents the value from `ValueArray` corresponding to the maximum intensity in `cube` along the third dimension.

# Description
This function calculates the indices of the maximum values in the third dimension of `cube` and maps these indices to values from `ValueArray`. The result is a 2D array where each element represents the value from `ValueArray` corresponding to the maximum intensity in `cube` for each (x, y) position.

# Example
```julia
# Example input arrays
cube = rand(100, 100, 50)  # Example data cube with dimensions 100x100x50
ValueArray = range(-500, stop=500, length=50)

# Function call
MapMaxIndices = MaxIndicesMap(cube, ValueArray)
"""
function MaxIndicesMap(cube::AbstractArray, ValueArray::AbstractArray)
    MapMaxIndices = zeros((size(cube,1),size(cube,2)))
    MaxIndices =  dropdims(argmax(cube, dims=3), dims=3)
    
    for i in 1:size(cube,1)
        for j in 1:size(cube,2)
            MapMaxIndices[i,j] = ValueArray[MaxIndices[i,j][3]]
        end
    end
    return MapMaxIndices
end

"""
    ask_user(prompt::String, default::Float64) -> Float64

Prompt the user for input with a default float value.

# Arguments
- `prompt::String`: The message to display to the user.
- `default::Float64`: The default float64 value to return if the user provides no input.
OR
- `default::Int`: The default integer value to return if the user provides no input.
OR
- `default::String`: The default string to return if the user provides no input.

# Returns
- `Float64`: The value entered by the user, or the default value if no input is provided.

# Description
This function prompts the user for input and returns the entered value as a `Float64`. If the user provides no input, the function returns the specified default value.

# Example
```julia
value = ask_user("Enter a float value", 3.14)
"""
function ask_user(prompt::String, default::Float64)
    println(prompt, "(default: ",default, "): ")
    val = readline()
    isempty(val) ? default : parse(Float64, val)
end
function ask_user(prompt::String, default::Int)
    println(prompt, "(default: ",default, "): ")
    val = readline()
    isempty(val) ? default : parse(Int, val)
end
function ask_user(prompt::String, default::String)
    println(prompt, "(default: ",default, "): ")
    response = readline()
    isempty(response) ? default : parse(String, response) 
end

"""
    contains_fits_files(dir::String) -> Bool

Check if a directory contains any FITS files.

# Arguments
- `dir::String`: The directory path to check for FITS files.

# Returns
- `Bool`: Returns `true` if the directory contains any files with the ".fits" extension, otherwise `false`.

# Description
This function checks if the specified directory contains any files with the ".fits" extension. It reads the directory contents and returns `true` if any of the files end with ".fits", otherwise it returns `false`.

# Example
```julia
# Example usage
directory_path = "path/to/directory"
has_fits_files = contains_fits_files(directory_path)
println("Contains FITS files: ", has_fits_files)
"""
function contains_fits_files(dir)
    return any(f -> endswith(f, ".fits"), readdir(dir))
end

"""
    get_simulation_list(base_dir::String) -> Vector{String}

Retrieve a list of directories containing FITS files within a base directory.

# Arguments
- `base_dir::String`: The base directory to search for simulation directories containing FITS files.

# Returns
- `Vector{String}`: A sorted vector of directory paths that contain FITS files.

# Description
This function searches recursively within the specified base directory for subdirectories that contain FITS files. It returns a sorted list of these directories.

# Example
```julia
# Example usage
base_directory = "path/to/base_directory"
simulation_list = get_simulation_list(base_directory)
println("Simulation directories containing FITS files: ", simulation_list)
"""
function get_simulation_list(base_dir)
    simulation_dirs = []
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

"""
    display_simulations(simu_list::Vector{String})

Display a list of available simulations.

# Arguments
- `simu_list::Vector{String}`: A vector of directory paths for available simulations.

# Returns
- `Nothing`: This function does not return a value. It prints the list of available simulations to the console.

# Description
This function takes a list of directory paths for available simulations and prints them to the console in a numbered format.

# Example
```julia
# Example usage
simulation_list = ["simu1", "simu2", "simu3"]
display_simulations(simulation_list)
# Output:
# Available simulations:
# [1] simu1
# [2] simu2
# [3] simu3
"""
function display_simulations(simu_list)
    println("Available simulations:")
    for (i, simu) in enumerate(simu_list)
        println("[$i] $simu")
    end
end

"""
    extract_category(dir::String) -> String

Extract the category from a directory name based on a specific pattern.

# Arguments
- `dir::String`: The directory name from which to extract the category.

# Returns
- `String`: The extracted category if the pattern matches, otherwise "unknown".

# Description
This function attempts to extract a category from the directory name using a regex pattern. If a match is found, it returns the captured category. If no match is found, it prints a message indicating that the simulation name was not found and returns "unknown".

# Example
```julia
# Example usage
dir_name = "d1cf123bx456rms789"
category = extract_category(dir_name)
println(category)  # Output: "123bx456rms789"

dir_name_invalid = "invalid_directory_name"
category = extract_category(dir_name_invalid)
# Output: "Nom de simulation non trouvé pour : invalid_directory_name"
println(category)  # Output: "unknown"
"""
function extract_category(dir)
    m = match(r"d1cf(\d+bx\d+rms\d+)", dir)
    if m !== nothing
        return m.captures[1]
    else
        println("Nom de simulation non trouvé pour : $dir")
        return "unknown"
    end
end

"""
    clean_path(path::String) -> String

Clean a file path by removing specific substrings based on predefined patterns.

# Arguments
- `path::String`: The file path to be cleaned.

# Returns
- `String`: The cleaned file path with specific substrings removed.

# Description
This function cleans a file path by removing specific substrings based on predefined patterns. It sequentially applies a series of replacements to remove unwanted parts of the path.

# Example
```julia
# Example usage
original_path = "/Users/jb270005/Desktop/simu_RAMSES/nograv1024/some_path/256/bx10d1otherparts"
cleaned_path = clean_path(original_path)
println(cleaned_path)  # Output: "some_path/otherparts"
"""
function clean_path(path)
    cleaned_path = replace(path, r"/Users/jb270005/Desktop/simu_RAMSES/" => "")
    cleaned_path = replace(cleaned_path, r"nograv1024/.+kyr/256" => "")
    cleaned_path = replace(cleaned_path, r"bx10" => "")
    cleaned_path = replace(cleaned_path, r"d1" => "") 
    return cleaned_path
end

"""
    print_progress(progress::Int, total::Int)

Print a progress bar to the console.

# Arguments
- `progress::Int`: The current progress value.
- `total::Int`: The total value indicating 100% completion.

# Description
This function prints a progress bar to the console to visually indicate the progress of a task. The progress bar is 50 characters wide and shows the ratio of `progress` to `total`.

# Example
```julia
# Example usage
for i in 1:100
    sleep(0.1)  # Simulate work
    print_progress(i, 100)
end
println()  # Move to the next line after completion
"""
function print_progress(progress::Int, total::Int)
    bar_width = 50
    progress_ratio = progress / total
    filled_length = Int(round(bar_width * progress_ratio))
    bar = "█" ^ filled_length * " " ^ (bar_width - filled_length)
    print("\rProgress: |$bar| $progress/$total")
end

function create_color_palette()
    return [
        "blue", "green", "red", "purple", "orange", 
        "cyan", "magenta", "brown", "pink", "gray", 
        "lime", "navy", "teal", "violet"
    ]
end