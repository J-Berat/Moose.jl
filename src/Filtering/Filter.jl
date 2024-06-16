"""
    HighPass(image, kernel)

Applies a high-pass filter to the input image.

# Arguments
- `image::AbstractArray`: The input image to be filtered.
- `kernel::AbstractArray`: The kernel used for filtering the image.

# Returns
- `imagef::AbstractArray`: The high-pass filtered image.

# Description
This function performs high-pass filtering by first applying a low-pass filter to the image using the provided `kernel` and then subtracting the filtered image from the original image. The result is an image where high-frequency components (edges, fine details) are enhanced.

# Example
```julia
using Images

image = rand(256, 256)  # Example image
kernel = ones(3, 3) / 9  # Example kernel for low-pass filtering

filtered_image = HighPass(image, kernel)
"""
function HighPass(image, kernel)
    filtered_image = imfilter(image, kernel)
    imagef = image .- filtered_image
    return imagef
end