"""
    buildHeader3D(naxis::Int, size::Tuple, ctype1::String, ctype2::String, ctype3::String, 
                  cunit1::String, cunit2::String, cunit3::String, bunit::String, specarray::AbstractArray) -> FITSHeader

Build a FITS header for a 3D data array.

# Arguments
- `naxis::Int`: The number of data axes.
- `size::Tuple`: The size of the 3D data array.
- `ctype1::String`: The type of the first axis.
- `ctype2::String`: The type of the second axis.
- `ctype3::String`: The type of the third axis.
- `cunit1::String`: The unit of the first axis.
- `cunit2::String`: The unit of the second axis.
- `cunit3::String`: The unit of the third axis.
- `bunit::String`: The unit of the data values.
- `specarray::AbstractArray`: An array containing spectral information.

# Returns
- `FITSHeader`: A FITS header object populated with the provided metadata.

# Description
This function constructs a FITS header for a 3D data array using the specified metadata. The header includes information about the data axes, their units, and the data value unit. The spectral information is derived from `specarray`.

# Example
```julia
# Example usage
naxis = 3
size = (100, 100, 50)
ctype1 = "RA---TAN"
ctype2 = "DEC--TAN"
ctype3 = "FREQ"
cunit1 = "deg"
cunit2 = "deg"
cunit3 = "Hz"
bunit = "Jy/beam"
specarray = [1.0, 2.0, 3.0]

header = buildHeader3D(naxis, size, ctype1, ctype2, ctype3, cunit1, cunit2, cunit3, bunit, specarray)
"""
function buildHeader3D(naxis, size, ctype1, ctype2, ctype3, cunit1, cunit2, cunit3, bunit, specarray)
    
    header = FITSHeader(["NAXIS"], [naxis], [""])

    header["NAXIS1"] = size[1]
    header["NAXIS2"] = size[2]

    header["CTYPE1"] = ctype1
    header["CRVAL1"] = 0
    header["CRPIX1"] = 1
    header["CDELT1"] = 1
    header["CUNIT1"] = cunit1

    header["CTYPE2"] = ctype2
    header["CRVAL2"] = 0
    header["CRPIX2"] = 1
    header["CDELT2"] = 1
    header["CUNIT2"] = cunit2

    header["NAXIS3"] = size[3]
    header["CTYPE3"] = ctype3
    header["CRVAL3"] = specarray[1]
    header["CRPIX3"] = 1
    header["CDELT3"] = specarray[2] - specarray[1]
    header["CUNIT3"] = cunit3
    header["BLENGTH"] = length(specarray)

    header["BUNIT"] = bunit

    return header
end

"""
    buildHeader2D(naxis::Int, size::Tuple, ctype1::String, ctype2::String, 
                  cunit1::String, cunit2::String, bunit::String) -> FITSHeader

Build a FITS header for a 2D data array.

# Arguments
- `naxis::Int`: The number of data axes.
- `size::Tuple`: The size of the 2D data array.
- `ctype1::String`: The type of the first axis.
- `ctype2::String`: The type of the second axis.
- `cunit1::String`: The unit of the first axis.
- `cunit2::String`: The unit of the second axis.
- `bunit::String`: The unit of the data values.

# Returns
- `FITSHeader`: A FITS header object populated with the provided metadata.

# Description
This function constructs a FITS header for a 2D data array using the specified metadata. The header includes information about the data axes, their units, and the data value unit.

# Example
```julia
# Example usage
naxis = 2
size = (100, 100)
ctype1 = "RA---TAN"
ctype2 = "DEC--TAN"
cunit1 = "deg"
cunit2 = "deg"
bunit = "Jy/beam"

header = buildHeader2D(naxis, size, ctype1, ctype2, cunit1, cunit2, bunit)
"""
function buildHeader2D(naxis, size, ctype1, ctype2, cunit1, cunit2, bunit)
    
    header = FITSHeader(["NAXIS"], [naxis], [""])

    header["NAXIS1"] = size[1]
    header["NAXIS2"] = size[2]

    header["CTYPE1"] = ctype1
    header["CRVAL1"] = 0
    header["CRPIX1"] = 1
    header["CDELT1"] = 1
    header["CUNIT1"] = cunit1

    header["CTYPE2"] = ctype2
    header["CRVAL2"] = 0
    header["CRPIX2"] = 1
    header["CDELT2"] = 1
    header["CUNIT2"] = cunit2

    header["BUNIT"] = bunit

    return header
end 