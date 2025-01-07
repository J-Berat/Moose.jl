"""
    Bperp(B1::AbstractArray, B2::AbstractArray) -> AbstractArray

Calculate the perpendicular component of the magnetic field for given magnetic field components B1 and B2.

# Arguments
- `B1::AbstractArray`: An array representing the first component of the magnetic field.
- `B2::AbstractArray`: An array representing the second component of the magnetic field.

# Returns
- `AbstractArray`: An array representing the perpendicular component of the magnetic field.

# Description
This function calculates the perpendicular component of the magnetic field using the formula:
B_perp = sqrt(B1^2 + B2^2)
where B1 and B2 are the magnetic field components.

# Example
```julia
# Example usage
B1 = rand(100, 100)  # Example magnetic field component B1
B2 = rand(100, 100)  # Example magnetic field component B2
B_perp = Bperp(B1, B2)
println(B_perp)
"""
Bperp(B1::AbstractArray, B2::AbstractArray) = @.  sqrt(B1^2 + B2^2)

"""
    Btot(Bx::AbstractArray, By::AbstractArray, Bz::AbstractArray) -> AbstractArray

Calculate the total magnetic field strength for given magnetic field components Bx, By, and Bz.

# Arguments
- `Bx::AbstractArray`: An array representing the x-component of the magnetic field.
- `By::AbstractArray`: An array representing the y-component of the magnetic field.
- `Bz::AbstractArray`: An array representing the z-component of the magnetic field.

# Returns
- `AbstractArray`: An array representing the total magnetic field strength.

# Description
This function calculates the total magnetic field strength using the formula:
Btot = sqrt(Bx^2 + By^2 + Bz^2)
where Bx, By, and Bz are the magnetic field components.

# Example usage
Bx = rand(100, 100)  # Example x-component of the magnetic field
By = rand(100, 100)  # Example y-component of the magnetic field
Bz = rand(100, 100)  # Example z-component of the magnetic field
B_total = Btot(Bx, By, Bz)
println(B_total)
"""
Btot(Bx::AbstractArray, By::AbstractArray, Bz::AbstractArray) = @. sqrt(Bx^2 + By^2 + Bz^2)

"""
    Bpulsar(RM::AbstractArray, DM::AbstractArray) -> AbstractArray

This function calculates the Bfield from pulsar rotation measure to dispersion measure ratio.

# Arguments
- `RM::AbstractArray`: An array representing the rotation measure.
- `DM::AbstractArray`: An array representing the dispersion measure.

# Returns
- `AbstractArray`: An array representing the pulsar magnetic field strength.

# Description
This function calculates the Bfield from pulsar rotation measure to dispersion measure ratio using the formula:
Bpulsar = 1.232 * RM ./ DM
where RM is the rotation measure and DM is the dispersion measure.

# Example usage
RM = rand(100)  # Example rotation measure data
DM = rand(100)  # Example dispersion measure data
B_pulsar = Bpulsar(RM, DM)
println(B_pulsar)

"""
Bpulsar(RM::AbstractArray, DM::AbstractArray) = 1.232 * RM ./ DM