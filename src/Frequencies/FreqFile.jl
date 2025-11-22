"""
    CreateFreqFile(start_freq_Hz::Number, end_freq_Hz::Number, num_freq::Int, repertory::AbstractString)

Creates a text file named `FreqHz.txt` containing a range of frequencies.

# Arguments
- `start_freq_Hz::Number`: The starting frequency in Hertz.
- `end_freq_Hz::Number`: The ending frequency in Hertz.
- `num_freq::Int`: The number of frequency points to generate between `start_freq_Hz` and `end_freq_Hz`.
- `repertory::AbstractString`: The directory where the file will be saved. If empty, the file is saved in the current directory.

# Description
This function generates a range of frequencies between the specified starting and ending frequencies. The number of points in the range is determined by `num_freq`. The function writes these frequencies line by line into a text file named `FreqHz.txt` in the specified directory. If no directory is provided, the file is saved in the current working directory.

# Example
```julia
# Create a file with 10 frequencies ranging from 100 Hz to 200 Hz (inclusive)
CreateFreqFile(100.0, 200.0, 10, "output_directory")

# File `output_directory/FreqHz.txt` will contain:
# 100.0
# 111.11111111111111
# 122.22222222222223
# ...
# 200.0
"""

function CreateFreqFile(start_freq_Hz, end_freq_Hz, num_freq, repertory)
    
    repertory = isempty(repertory) ? "" : repertory * "/"

    freq_vector = range(start_freq_Hz, stop=end_freq_Hz, length=num_freq)
    
    path = joinpath(repertory, "FreqHz.txt")

    open(path, "w") do file
        for freq in freq_vector
            println(file, freq)
        end
    end

    println("Your file $path has been created.")
    
end