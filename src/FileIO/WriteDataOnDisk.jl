"""
    WriteData3D(resultspath::String, data::AbstractArray, DataName::String, specarray::AbstractArray)

Write a 3D data array to a FITS file with appropriate headers.

# Arguments
- `resultspath::String`: The directory path where the FITS file will be saved.
- `data::AbstractArray`: The 3D data array to be written to the FITS file.
- `DataName::String`: The name of the data, used for naming the FITS file and retrieving header information.
- `specarray::AbstractArray`: An array containing spectral information to be included in the FITS file header.

# Description
This function writes a 3D data array to a FITS file in the specified directory with headers built from predefined metadata. It ensures the target directory exists, constructs the appropriate header using `buildHeader3D`, and writes the data to the FITS file.

# Example
```julia
# Example usage
resultspath = "path/to/results"
data = rand(100, 100, 50)  # Example 3D data array
DataName = "ExampleData"
specarray = [1.0, 2.0, 3.0]  # Example spectral array

WriteData3D(resultspath, data, DataName, specarray)
```
"""
function WriteData3D(resultspath::String, data::AbstractArray, DataName::String, specarray::AbstractArray)

    # Path
    mkpath(resultspath)

    params = header_params(
        naxis=DictHeader[DataName]["naxis"],
        ctype1=DictHeader[DataName]["ctype1"],
        ctype2=DictHeader[DataName]["ctype2"],
        ctype3=DictHeader[DataName]["ctype3"],
        cunit1=DictHeader[DataName]["cunit1"],
        cunit2=DictHeader[DataName]["cunit2"],
        cunit3=DictHeader[DataName]["cunit3"],
        bunit=DictHeader[DataName]["bunit"],
    )

    header = buildHeader3D(
        params["naxis"],
        size(data),
        params["ctype1"],
        params["ctype2"],
        params["ctype3"],
        params["cunit1"],
        params["cunit2"],
        params["cunit3"],
        params["bunit"],
        specarray
    )

    fits_path = joinpath(resultspath, "$DataName.fits")
    FITS(fits_path, "w") do f
        write(f, data; header=header)
    end

    println("The FITS file of $DataName has been written in this directory: $fits_path")
end

"""
    WriteData2D(resultspath::String, data::AbstractArray, DataName::String)

Write a 2D data array to a FITS file with appropriate headers.

# Arguments
- `resultspath::String`: The directory path where the FITS file will be saved.
- `data::AbstractArray`: The 2D data array to be written to the FITS file.
- `DataName::String`: The name of the data, used for naming the FITS file and retrieving header information.

# Description
This function writes a 2D data array to a FITS file in the specified directory with headers built from predefined metadata. It ensures the target directory exists, constructs the appropriate header using `buildHeader2D`, and writes the data to the FITS file.

# Example
```julia
# Example usage
resultspath = "path/to/results"
data = rand(100, 100)  # Example 2D data array
DataName = "ExampleData"

WriteData2D(resultspath, data, DataName)
```
"""
function WriteData2D(resultspath::String, data::AbstractArray, DataName::String)

    # Path
    mkpath(resultspath)

    params = header_params(
        naxis=DictHeader[DataName]["naxis"],
        ctype1=DictHeader[DataName]["ctype1"],
        ctype2=DictHeader[DataName]["ctype2"],
        cunit1=DictHeader[DataName]["cunit1"],
        cunit2=DictHeader[DataName]["cunit2"],
        bunit=DictHeader[DataName]["bunit"],
    )

    header = buildHeader2D(
        params["naxis"],
        size(data),
        params["ctype1"],
        params["ctype2"],
        params["cunit1"],
        params["cunit2"],
        params["bunit"],
    )

    fits_path = joinpath(resultspath, "$DataName.fits")
    FITS(fits_path, "w") do f
        write(f, data; header=header)
    end

    println("The FITS file of $DataName has been written in this directory: $fits_path")

end
