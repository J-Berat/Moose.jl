"""
    moments(y::AbstractArray; x::AbstractArray=1:length(y)) -> Tuple{Float64, Float64, Float64}

Calculate the zeroth, first, and second moments of a given array.

# Arguments
- `y::AbstractArray`: The input array for which the moments are to be calculated.
- `x::AbstractArray=1:length(y)`: The x-values corresponding to the y-values. Defaults to an array from 1 to the length of `y`.

# Returns
- `Tuple{Float64, Float64, Float64}`: A tuple containing:
  - `m0::Float64`: The zeroth moment (sum of `y`).
  - `m1::Float64`: The first moment (mean of `x` weighted by `y`).
  - `m2::Float64`: The second moment (standard deviation of `x` weighted by `y`).

# Description
This function calculates the zeroth, first, and second moments of the input array `y`. The zeroth moment (`m0`) is the sum of `y`. The first moment (`m1`) is the mean of `x` weighted by `y`. The second moment (`m2`) is the standard deviation of `x` weighted by `y`, ensuring non-negative values.

# Example
```julia
# Example usage
y = [1, 2, 3, 4, 5]
x = [1, 2, 3, 4, 5]
m0, m1, m2 = moments(y, x=x)
println("Zeroth moment: ", m0)
println("First moment: ", m1)
println("Second moment: ", m2)
"""

function moments(y; x=1:length(y), threshold=0.0)
    mask = y .> threshold
    y_masked = y[mask]
    x_masked = x[mask]

    if isempty(y_masked)
        return NaN, NaN, NaN
    end

    m0 = sum(y_masked)
    m1 = sum(y_masked .* x_masked) / m0
    m2_2 = sum(@. y_masked * (x_masked - m1)^2) / m0
    m2 = m2_2 >= 0 ? sqrt(m2_2) : NaN

    return m0, m1, m2
end
