"""
    CreateFreqFile(start_freq_Hz, end_freq_Hz, num_freq, repertory)

Creates a text file named `FreqHz.txt` containing a range of frequencies.

# Arguments
- `start_freq_Hz::Number`: The starting frequency in Hertz.
- `end_freq_Hz::Number`: The ending frequency in Hertz.
- `num_freq::Int`: The number of frequency points to generate between `start_freq_Hz` and `end_freq_Hz`.
- `repertory::AbstractString`: The directory where the file will be saved. If empty, the file is saved in the current directory.

# Example
```julia
CreateFreqFile(100.0, 200.0, 10, "output_directory")
"""

function CreateFreqFile(start_freq_Hz, end_freq_Hz, num_freq, repertory)
    
    repertory = isempty(repertory) ? "" : repertory * "/"

    freq_vector = range(start_freq_Hz, end_freq_Hz, num_freq)
    
    path = joinpath(repertory, "FreqHz.txt")

    open(path, "w") do file
        for freq in freq_vector
            println(file, freq)
        end
    end

    println("Your file $path has been created.")
    
end