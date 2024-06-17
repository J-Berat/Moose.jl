"""
    CreateFreqFile(start_freq_Hz::Float64, end_freq_Hz::Float64, num_freq::Int, repertory::String="")

Create a text file containing a frequency vector within the specified range.

# Arguments
- `start_freq_Hz`: Starting frequency in Hertz.
- `end_freq_Hz`: Ending frequency in Hertz.
- `num_freq`: Number of frequencies in the vector.
- `repertory`: Directory path where the file will be created. If empty, the file will be created in the current working directory.
"""
function CreateFreqFile(start_freq_Hz::Float64, end_freq_Hz::Float64, num_freq::Int, repertory::String="")
    # Input validation
    if start_freq_Hz < 0 || end_freq_Hz < 0
        throw(ArgumentError("Frequencies must be non-negative."))
    end

    if num_freq <= 0
        throw(ArgumentError("Number of frequencies must be a positive integer."))
    end

    # Set the directory path
    directory_path = isempty(repertory) ? "." : repertory

    # Generate the frequency vector
    freq_vector = range(start_freq_Hz, stop=end_freq_Hz, length=num_freq)
    
    # Define the full file path
    file_path = joinpath(directory_path, "FreqHz.txt")

    # Open the file and write the frequency vector
    open(file_path, "w") do file
        for freq in freq_vector
            println(file, freq)
        end
    end

    println("Your file $file_path has been created.")
end