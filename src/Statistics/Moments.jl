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

function moments(y::AbstractArray; x::AbstractArray=1:length(y)) 
    
    m0 = sum(y)
    m1 = sum(y .* x) / m0
    m2_2 = sum(@. y * (x-m1)^2) / m0
    m2 = m2_2 >= 0 ? sqrt(m2_2) : NaN
    
    return m0, m1, m2    
end

"""
    MomentsofMomentMap(cube::AbstractArray; x::AbstractArray=1:length(cube[1,1,:])) -> Tuple{Float64, Float64, Float64}

Calculate statistical moments of the second moment map of a data cube.

# Arguments
- `cube::AbstractArray`: A 3D array representing the data cube.
- `x::AbstractArray=1:length(cube[1,1,:])`: The x-values corresponding to the third dimension of the cube. Defaults to an array from 1 to the length of the third dimension of the cube.

# Returns
- `Tuple{Float64, Float64, Float64}`: A tuple containing:
  - `m1_M2::Float64`: The mean of the second moment map.
  - `m2_M2::Float64`: The standard deviation of the second moment map.
  - `max_M2::Float64`: The maximum value of the second moment map.

# Description
This function calculates the second moment (variance) for each spatial position (i, j) in a 3D data cube along the third dimension. It creates a 2D map of these second moments (M2map) and then computes the mean, standard deviation, and maximum value of this map.

The second moment is calculated using the `moments` function, which returns the zeroth, first, and second moments of a given array.

# Example
```julia
# Example usage
cube = rand(100, 100, 50)  # Example 3D data cube
m1_M2, m2_M2, max_M2 = MomentsofMomentMap(cube)
println("Mean of M2 map: ", m1_M2)
println("Standard deviation of M2 map: ", m2_M2)
println("Maximum value of M2 map: ", max_M2)
"""
function MomentsofMomentMap(cube; x=1:length(cube[1,1,:]))
    M2map = zeros((size(cube,1),size(cube,2)))
    for i in 1:size(M2map,1)
        for j in 1:size(M2map,2)
            M2map[i,j] = moments(cube[i,j,:]; x=x)[3]
        end
    end
    m1_M2, m2_M2, max_M2 = mean(M2map), std(M2map), maximum(M2map)
    return m1_M2, m2_M2, max_M2
end
