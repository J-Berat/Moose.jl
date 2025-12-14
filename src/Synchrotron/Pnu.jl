"""
    Pnu(Q::AbstractArray, U::AbstractArray) -> AbstractArray

Calculate the magnitude of the Stokes parameters `Q` and `U` elementwise,
using `hypot` for improved numerical stability.

# Arguments
- `Q::AbstractArray`: An array of any size representing the Stokes parameter `Q`.
- `U::AbstractArray`: An array of any size representing the Stokes parameter `U`.

# Returns
- `AbstractArray`: An array of the same size as the input arrays representing the magnitude of the Stokes parameters.

# Example
```julia
Q = [1.0, 2.0, 3.0]
U = [4.0, 5.0, 6.0]
P = Pnu(Q, U)
# P should be [4.123105625617661, 5.385164807134504, 6.708203932499369]
```
"""
Pnu(Q, U) = hypot.(Q, U)

