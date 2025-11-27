"""
Collection of array-level helpers used throughout the pipeline.
"""

maxCube(cube::AbstractArray) = dropdims(maximum(cube, dims = 3), dims = 3)
intLOS(cube::AbstractArray, PixelLength_cm::Float64) = dropdims(sum(cube .* PixelLength_cm, dims = 3), dims = 3)
sigmaLOS(cube::AbstractArray) = dropdims(std(cube, dims = 3), dims = 3)
