"""
    fbm2d(exponent::Float64, nx::Int, ny::Int; sigma::Float64=1.0, avg::Float64=0.0, positive::Bool=false) -> AbstractArray

Generate a 2D fractional Brownian motion (fBm) image.

# Arguments
- `exponent::Float64`: The exponent controlling the roughness of the fBm.
- `nx::Int`: The number of pixels along the x-axis.
- `ny::Int`: The number of pixels along the y-axis.
- `sigma::Float64=1.0`: The standard deviation of the output image. Default is 1.0.
- `avg::Float64=0.0`: The average value of the output image. Default is 0.0.
- `positive::Bool=false`: If true, ensures the output image has only positive values.

# Returns
- `AbstractArray`: A 2D array representing the fractional Brownian motion image.

# Description
This function generates a 2D fractional Brownian motion (fBm) image based on the given parameters. It starts with a white noise image, applies a Fourier transform, filters it using a power-law filter determined by the exponent, and then performs an inverse Fourier transform to obtain the fBm image. The image is then normalized and adjusted to have the specified standard deviation (`sigma`) and average value (`avg`). If `positive` is set to true and `avg` is greater than 0, the image is transformed to ensure all values are positive.

# Example
```julia
# Example usage
exponent = 2.0
nx, ny = 256, 256
sigma = 1.0
avg = 0.0
positive = false

image = fbm2d(exponent, nx, ny, sigma=sigma, avg=avg, positive=positive)
using Plots
heatmap(image, color=:viridis)
"""

function fbm2d(exponent::Float64, nx::Int, ny::Int; sigma::Float64=1.0, avg::Float64=0.0, positive::Bool=false)

    # initial white noise image and its FFT
    im0 = rand(nx,ny)
    im0f = fft(im0)

    # KMAT
    kx = fftfreq(nx,1)
    ky = fftfreq(ny,1)
    kmat = sqrt.((kx .^ 2) .+ (ky' .^ 2))
   
    # AMPLITUDE
    amplitude = kmat .^ (exponent / 2.0)
    pos_center = findall(@. (kmat == 0))
    amplitude[pos_center] .= 0
    
    # Filter the white noise image
    im0f = amplitude .* im0f

    # BACK TO REAL SPACE - do the Inverse FFT (ifft)
    image = real(ifft(im0f))
    
    # Normalization
    image /= std(image)
    
    if positive && avg > 0
        bval = log10.(range(0.1, stop=3.0, length=10))
        val = [std(exp(b * image)) / avg(exp(b * image)) for b in bval]
        B = LinearInterpolation(bval, val)(sigma / avg)
        A = sigma / std(exp(B .* image))
        image = A .* exp.(B .* image)
    else
        image *= sigma
        image .+= avg
    end
    
    return image
end