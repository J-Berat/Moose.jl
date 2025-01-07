"""
    Get3rdAxisValue(fitsfile::String, frame::Int) -> Float64

Calculate the value along the third axis of a FITS file for a given frame.

# Arguments
- `fitsfile::String`: The path to the FITS file.
- `frame::Int`: The frame index for which to calculate the value along the third axis.

# Returns
- `Float64`: The calculated value along the third axis for the specified frame.

# Description
This function calculates the value along the third axis of a FITS file for a specified frame. FITS files use a reference position (`CRPIX`), a reference value (`CRVAL`), and a step size (`CDELT`) to determine the pixel coordinate values. The function reads these values from the FITS header and computes the value for the given frame index.

# Example
```julia
# Example usage
fitsfile = "path/to/your/file.fits"
frame = 10
value = Get3rdAxisValue(fitsfile, frame)
"""

function Get3rdAxisValue(fitsfile,frame)
    # Open file
    i = frame
    f = FITS(fitsfile)
    header = read_header(f[1])
    #Calculate value. FITS uses a reference position (CRPIX) 
    #and value (CRVAL), and a step size (CDELT) to determine the pixel coordinate values.
    FDEP_slice = header["CRVAL3"]+header["CDELT3"]*(i-header["CRPIX3"]+1)
end