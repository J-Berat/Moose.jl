"""
    RMSynthesis(Q::AbstractArray, U::AbstractArray, nuArray::AbstractArray, PhiArray::AbstractArray) -> Tuple{AbstractArray, AbstractArray, AbstractArray}

Perform Rotation Measure (RM) synthesis on Stokes Q and U parameters.

# Arguments
- `Q::AbstractArray`: Array representing the Stokes Q parameter. Can be 1D, 2D, or 3D.
- `U::AbstractArray`: Array representing the Stokes U parameter. Can be 1D, 2D, or 3D.
- `nuArray::AbstractArray`: Array of frequency values in Hz.
- `PhiArray::AbstractArray`: Array of Faraday depths in rad/m².

# Returns
- `Tuple{AbstractArray, AbstractArray, AbstractArray}`: Three arrays representing the absolute value, real part, and imaginary part of the Faraday dispersion function.

# Description
This function performs RM synthesis, a technique used in radio astronomy to study the Faraday rotation effect. The input Stokes Q and U parameters are combined to form the complex polarization P. The function then calculates the Faraday dispersion function F for each value in `PhiArray`.

# Example
```julia
# Example input arrays
Q = [0.1, 0.2, 0.3]
U = [0.4, 0.5, 0.6]
nuArray = [1e9, 1.1e9, 1.2e9]
PhiArray = [-100, 0, 100]

# Function call
absF, realF, imagF = RMSynthesis(Q, U, nuArray, PhiArray)
"""

function RMSynthesis(Q::AbstractArray, U::AbstractArray, nuArray::AbstractArray, PhiArray::AbstractArray)
    
    LambdaSqArray = @. (C_m/nuArray)^2
    
    nPhi = length(PhiArray)
    nLambda = length(LambdaSqArray)
    nDims = length(size(Q))
    
    WeightArray = ones(nLambda)
    K = 1.0 / sum(WeightArray)

    if nDims == 1
        Q = reshape(Q, (1,1,size(Q,1)))
        U = reshape(U, (1,1,size(U,1)))
    elseif nDims == 2
        Q = reshape(Q, (1,size(Q,1), size(Q,2)))
        U = reshape(U, (1,size(U,1), size(U,2)))
    end

    P = @. (Q + 1im * U) * WeightArray[[CartesianIndex()],[CartesianIndex()],:]

    nx, ny = size(Q,1), size(Q,2)
    
    Lambda0Sq = K .* sum(WeightArray .* LambdaSqArray)
    a = (LambdaSqArray .- Lambda0Sq)

    F = complex(zeros(nx,ny,nPhi))
    arg = complex(zeros(nPhi))
    
    for i in 1:nPhi
        arg = exp.((-2.0im .* PhiArray[i]) .* a)[[CartesianIndex()],[CartesianIndex()],:]
        F[:,:,i] = K .* sum(P .* arg, dims=3)
    end
    
    if nDims == 1
        F = dropdims(dropdims(F,dims=1),dims=1)
    elseif nDims == 2
        F = dropdims(F,dims=1)
    end
     
    return(abs.(F),real.(F),imag.(F))
end

"""
    getRMSF(nuArray::AbstractArray, PhiArray::AbstractArray) -> Tuple{AbstractArray, Float64}

Calculate the Rotation Measure Spread Function (RMSF) and its full width at half maximum (FWHM).

# Arguments
- `nuArray::AbstractArray`: Array of frequency values in Hz.
- `PhiArray::AbstractArray`: Array of Faraday depths in rad/m².

# Returns
- `Tuple{AbstractArray, Float64}`: 
  - `AbstractArray`: An array representing the absolute value of the RMSF.
  - `Float64`: The full width at half maximum (FWHM) of the RMSF.

# Description
This function calculates the Rotation Measure Spread Function (RMSF) for a given set of frequency values and Faraday depths. The RMSF describes the response of an RM synthesis to a single Faraday depth component. The function also computes the full width at half maximum (FWHM) of the RMSF.

# Example
```julia
# Example input arrays
nuArray = [1e9, 1.1e9, 1.2e9]
PhiArray = [-100, 0, 100]

# Function call
absRMSF, fwhmRMSF = getRMSF(nuArray, PhiArray)
"""
function getRMSF(nuArray::AbstractArray, PhiArray::AbstractArray)
    
    LambdaSqArray = @. (C_m/nuArray)^2
    
    nPhi = length(PhiArray)
    nLambda = length(LambdaSqArray)
    
    WeightArray = ones(nLambda)
    K = 1.0 / sum(WeightArray)
    
    Lambda0Sq = K .* sum(WeightArray .* LambdaSqArray)
    a = (LambdaSqArray .- Lambda0Sq)
    
    fwhmRMSF = 3.8 / (maximum(LambdaSqArray) - minimum(LambdaSqArray))

    RMSF = complex(zeros(nPhi))
    arg = complex(zeros(nPhi))
    for i in 1:nPhi
        arg = exp.((-2.0im .* PhiArray[i]) .* a)
        RMSF[i] = K .* sum(WeightArray .* arg)
    end
     
    return(abs.(RMSF),fwhmRMSF)
end