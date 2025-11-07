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

function faraday_moments(phi, P, noise_level; upsample=10)
    dphi = phi[2] - phi[1]

    bias = length(P) * noise_level
    M0 = (sum(P) - bias) * dphi

    Weff = M0 / (maximum(P) - noise_level)
    M2 = Weff / 2.354

    nP = length(P)
    itp = LinearInterpolation(1:nP, P)
    dx = 1 / upsample
    xup = 1:dx:nP
    Pup = itp(xup)

    k = Kernel.gaussian((M2 / dphi / dx,))
    Pconv = imfilter(Pup, k)
    imax = argmax(Pconv) * dx

    itp_phi = LinearInterpolation(1:nP, phi, extrapolation_bc=Line())
    M1 = itp_phi(imax)

    return M0, M1, Weff
end

function moments_map(data, array; threshold=0.0)
    M0 = zeros(size(data, 1), size(data, 2))
    M1 = zeros(size(data, 1), size(data, 2))
    M2 = zeros(size(data, 1), size(data, 2))
    for i in 1:size(data, 1)
        for j in 1:size(data, 2)
            m0, m1, m2 = moments(data[i,j,:], x=array,threshold=0.0)
            M0[i,j] = m0
            M1[i,j] = m1
            M2[i,j] = m2
        end
    end
    return M0, M1, M2
end