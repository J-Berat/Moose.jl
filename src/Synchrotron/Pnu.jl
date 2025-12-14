"""
    Pnu(Q, U) -> AbstractArray

Calculate the magnitude of the Stokes parameters `Q` and `U` elementwise,
using `hypot` for improved numerical stability.

# Arguments
- `Q`: A scalar or array representing the Stokes parameter `Q`.
- `U`: A scalar or array representing the Stokes parameter `U`.

`Q` and `U` only need to be broadcastable with one another; they need not be
the same shape.

# Returns
- `AbstractArray`: An array representing the magnitude of the Stokes parameters.

# Example
```julia
Q = [1.0, 2.0, 3.0]
U = [4.0, 5.0, 6.0]
P = Pnu(Q, U)
# P should be [4.123105625617661, 5.385164807134504, 6.708203932499369]
```
"""
Pnu(Q, U) = hypot.(Q, U)

