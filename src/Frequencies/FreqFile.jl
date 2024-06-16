"""
    create_freq_file(start_freq_Hz, end_freq_Hz, num_freq, repertory)

Create a text file containing a frequency vector within the specified range.

# Arguments
- `start_freq_Hz`: Starting frequency in Hertz.
- `end_freq_Hz`: Ending frequency in Hertz.
- `num_freq`: Number of frequencies in the vector.
- `repertory`: Directory path where the file will be created. If empty, the file will be created in the current working directory.
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