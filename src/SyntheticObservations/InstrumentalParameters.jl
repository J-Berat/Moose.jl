"""
    FrequencyParameters() -> AbstractArray

Prompt the user to enter the frequency range and resolution, and compute the frequency array.

# Returns
- `AbstractArray`: An array of frequency values in MHz.

# Description
This function prompts the user to provide:
- The starting frequency (`nustart`) in MHz.
- The ending frequency (`nuend`) in MHz.
- The frequency resolution (`dnu`) in MHz.

It computes a range of frequencies based on these values.

# Example
```julia
# Example usage
nuArray = FrequencyParameters()

# Sample interaction
Frequency range start (MHz): 115
Frequency range end (MHz): 175
Frequency resolution (MHz): 0.2

# Result
nuArray = [115.0, 115.2, ..., 174.8, 175.0]
"""
function FrequencyParameters()
    println("Please enter the values for the parameters:")
    nustart = ask_user("Frequency range start (MHz)", 120)
    nuend = ask_user("Frequency range end (MHz)", 167)
    dnu = ask_user("Frequency resolution (MHz)", 0.098)
    nuArray = range(start=nustart, stop=nuend, step=dnu)
    
    return nuArray
end

"""
    FaradayParameters() -> AbstractArray

Prompt the user to enter the Faraday depth range and resolution, and compute the Faraday depth array.

# Returns
- `AbstractArray`: An array of Faraday depth values in rad/m².

# Description
This function prompts the user to provide:
- The starting Faraday depth (`Phistart`) in rad/m².
- The ending Faraday depth (`Phiend`) in rad/m².
- The Faraday depth resolution (`dPhi`) in rad/m².

It computes a range of Faraday depth values based on these inputs.

# Example
```julia
# Example usage
PhiArray = FaradayParameters()

# Sample interaction
Faraday depth range start (rad/m^2): -20
Faraday depth range end (rad/m^2): 20
Faraday depth resolution (rad/m^2): 0.1

# Result
PhiArray = [-20.0, -19.9, ..., 19.9, 20.0]
"""
function FaradayParameters()
    println("Please enter the values for the parameters:")
    Phistart = ask_user("Faraday depth range start (rad/m^2)", -10)
    Phiend = ask_user("Faraday depth range end (rad/m^2)", 10)
    dPhi = ask_user("Faraday depth resolution (rad/m^2)", 0.25)
    PhiArray = range(start=Phistart, stop=Phiend, step=dPhi)

    return PhiArray
end

"""
    DistanceParameters() -> Tuple{Float64, Float64, Float64, AbstractArray}

Prompt the user to enter box size parameters and compute the related distances and pixel lengths.

# Returns
- `Tuple{Float64, Float64, Float64, AbstractArray}`: A tuple containing:
  - `PixelLength_pc::Float64`: Pixel length in parsecs.
  - `PixelLength_cm::Float64`: Pixel length in centimeters.
  - `BoxLength_pc::Float64`: Side length of the simulation box in parsecs.
  - `DistanceArray::AbstractArray`: Array of distances from the box center in parsecs.

# Description
This function prompts the user to provide:
- The side length of the box in parsecs (`BoxLength_pc`).
- The resolution of the box in pixels (`BoxLength_pix`).

It computes:
1. The pixel length in parsecs and centimeters.
2. An array of distances across the box.

# Example
```julia
# Example usage
PixelLength_pc, PixelLength_cm, BoxLength_pc, DistanceArray = DistanceParameters()

# Sample interaction
Side of the Box size (pc), please give a Float: 50.0
Side of the Box size (pixel): 256

# Result
PixelLength_pc = 0.1953125
PixelLength_cm = 6.03281e17
BoxLength_pc = 50.0
DistanceArray = [0.0, 0.1953125, ..., 49.8046875, 50.0]
"""
function DistanceParameters()   
    println("Please enter the values for the parameters:")
    BoxLength_pc = ask_user("Side of the Box size (pc), please give a Float", 50.)
    BoxLength_pix = ask_user("Side of the Box size (pixel)", 256)
    PixelLength_pc = BoxLength_pc / BoxLength_pix
    PixelLength_cm = PixelLength_pc * PARSEC_TO_CM
    
    Dstart = 0
    Dend = BoxLength_pc
    dD = PixelLength_pc
    DistanceArray = range(start=Dstart, stop=Dend, step=dD)
    
    return PixelLength_pc, PixelLength_cm, BoxLength_pc, DistanceArray
end

"""
    VelocityParameters() -> Tuple{AbstractArray, Float64, Float64, Float64}

Prompt the user to enter HI instrumental parameters and compute related arrays.

# Returns
- `Tuple{AbstractArray, Float64, Float64, Float64}`: A tuple containing:
  - `velArray::AbstractArray`: Array of velocity values in km/s.
  - `PixelLength_pc::Float64`: Pixel length in parsecs.
  - `PixelLength_cm::Float64`: Pixel length in centimeters.
  - `BoxLength_pc::Float64`: Side length of the simulation box in parsecs.

# Description
This function prompts the user to enter values for various HI instrumental parameters, including velocity range, velocity resolution, box size in parsecs, and box size in pixels. It then computes and returns the velocity array, pixel length in parsecs, pixel length in centimeters, and the box length in parsecs.

# Example
```julia
# Example usage
velArray, PixelLength_pc, PixelLength_cm, BoxLength_pc = HIInstrumentalParameters()
"""
function VelocityParameters()
    println("Please enter the values for the parameters:")
    velstart = ask_user("Velocity range start (km/s)", -30)
    velend = ask_user("Velocity range end (km/s)", 30)
    dvel = ask_user("Velocity resolution (km/s)", 1.29)
    velArray = range(start=velstart, stop=velend, step=dvel)

    return velArray
end
