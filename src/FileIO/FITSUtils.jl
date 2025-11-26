"""
Helpers for reading FITS files and aligning data cubes with the chosen line of sight.
"""

read_FITS_file(file) = read(FITS(file)[1])
read_file(file, conversion) = read_FITS_file(file) .* conversion

function permute_dims(array, LOS)
    if LOS == "x"
        permutedims(array, [2, 3, 1])
    elseif LOS == "y"
        permutedims(array, [3, 1, 2])
    else
        array
    end
end

read_optional_file(file, conversion, LOS) = isfile(file) ? permute_dims(read_file(file, conversion), LOS) : nothing
