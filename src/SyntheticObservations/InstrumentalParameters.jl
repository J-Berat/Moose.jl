"""
    SynchrotronInstrumentalParameters() -> Tuple{AbstractArray, AbstractArray, Float64, Float64, Float64}

Prompt the user to enter synchrotron instrumental parameters and compute related arrays.

# Returns
- `Tuple{AbstractArray, AbstractArray, Float64, Float64, Float64}`: A tuple containing:
  - `nuArray::AbstractArray`: Array of frequency values in MHz.
  - `PhiArray::AbstractArray`: Array of Faraday depth values in rad/m².
  - `PixelLength_pc::Float64`: Pixel length in parsecs.
  - `PixelLength_cm::Float64`: Pixel length in centimeters.
  - `BoxLength_pc::Float64`: Side length of the simulation box in parsecs.

# Description
This function prompts the user to enter values for various synchrotron instrumental parameters, including frequency range, frequency resolution, Faraday depth range, Faraday depth resolution, box size in parsecs, and box size in pixels. It then computes and returns the frequency array, Faraday depth array, pixel length in parsecs, pixel length in centimeters, and the box length in parsecs.

# Example
```julia
# Example usage
nuArray, PhiArray, PixelLength_pc, PixelLength_cm, BoxLength_pc = SynchrotronInstrumentalParameters()
"""
function SynchrotronInstrumentalParameters()

    println("Please enter the values for the parameters:")
    nustart = ask_user("Frequency range start (MHz)", 115)
    nuend = ask_user("Frequency range end (MHz)", 175)
    dnu = ask_user("Frequency resolution (MHz)", 0.2)
    Phistart = ask_user("Faraday depth range start (rad/m^2)", -20)
    Phiend = ask_user("Faraday depth range end (rad/m^2)", 20)
    dPhi = ask_user("Faraday depth resolution (rad/m^2)", 0.1)
    BoxLength_pc = ask_user("Side of the Box size (pc), please give a Float", 50.)
    BoxLength_pix = ask_user("Side of the Box size (pixel)", 256)

    nuArray = range(start=nustart, stop=nuend, step=dnu)
    PhiArray = range(start=Phistart, stop=Phiend, step=dPhi)

    PixelLength_pc = BoxLength_pc / BoxLength_pix
    PixelLength_cm = PixelLength_pc * PARSEC_TO_CM
    
    return nuArray, PhiArray, PixelLength_pc, PixelLength_cm, BoxLength_pc
end

"""
    HIInstrumentalParameters() -> Tuple{AbstractArray, Float64, Float64, Float64}

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
function HIInstrumentalParameters()
    
    println("Please enter the values for the parameters:")
    velstart = ask_user("Velocity range start (km/s)", -40)
    velend = ask_user("Velocity range end (km/s)", 40)
    dvel = ask_user("Velocity resolution (km/s)", 0.4)
    BoxLength_pc = ask_user("Side of the Box size (pc), please give a Float", 50.)
    BoxLength_pix = ask_user("Side of the Box size (pixel)", 256)
    velArray = range(start=velstart, stop=velend, step=dvel)

    PixelLength_pc = BoxLength_pc / BoxLength_pix
    PixelLength_cm = PixelLength_pc * PARSEC_TO_CM

    return velArray, PixelLength_pc, PixelLength_cm, BoxLength_pc
end
