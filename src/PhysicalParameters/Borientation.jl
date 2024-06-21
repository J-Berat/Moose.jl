"""
This function computes the angle between the line of sight magnetic field (BLOS) and the total magnetic field (Btot).

# Parameters:
- `BLOS`: Line of sight magnetic field component.
- `Btot`: Total magnetic field magnitude.

# Returns:
- Angle between BLOS and Btot in radians.
"""
Borientation(BLOS,Btot) = @. acos(BLOS / Btot)