"""
Collection of array-level helpers used throughout the pipeline.
"""

maxCube(cube::AbstractArray) = dropdims(maximum(cube, dims = 3), dims = 3)
intLOS(cube::AbstractArray, PixelLength_cm::Float64) = dropdims(sum(cube .* PixelLength_cm, dims = 3), dims = 3)
sigmaLOS(cube::AbstractArray) = dropdims(std(cube, dims = 3), dims = 3)

MeanSpectrum(cube::AbstractArray) = dropdims(mean(cube, dims = (1, 2)), dims = (1, 2))

function nearest_index(a, L)
    distances = abs.(L .- a)
    return argmin(distances)
end

function HOGbox(MeanSpectrumHI, MeanSpectrumFDF, VelArray, PhiArray)
    HImean = nearest_index(moments(MeanSpectrumHI, x = VelArray)[2], VelArray)
    FDFmean = nearest_index(moments(MeanSpectrumFDF, x = PhiArray)[2], PhiArray)

    HIboxline = nearest_index(moments(MeanSpectrumHI, x = VelArray)[3], VelArray)
    FDFboxline = nearest_index(moments(MeanSpectrumFDF, x = PhiArray)[3], PhiArray)

    return HImean, FDFmean, HIboxline, FDFboxline
end

function MaxIndicesMap(cube::AbstractArray, ValueArray::AbstractArray)
    MapMaxIndices = zeros((size(cube, 1), size(cube, 2)))
    MaxIndices = dropdims(argmax(cube, dims = 3), dims = 3)

    for i in 1:size(cube, 1)
        for j in 1:size(cube, 2)
            MapMaxIndices[i, j] = ValueArray[MaxIndices[i, j][3]]
        end
    end
    return MapMaxIndices
end

logindgen(nb, minv, maxv) = 10 .^ range(log10(minv), log10(maxv), length = nb)

function AreaUnderCurve(x, y)
    length(x) == length(y) || error("Vectors x and y must have the same length.")
    return sum(diff(x) .* ((y[1:end-1] + y[2:end]) ./ 2))
end

function AreaBetweenCurves(x, y1, y2)
    if length(x) != length(y1) || length(x) != length(y2)
        error("Vectors x, y1, and y2 must have the same length.")
    end

    sum(diff(x) .* abs.((y1[2:end] .+ y1[1:end-1] .- y2[2:end] .- y2[1:end-1]) ./ 2))
end

vectornorm(Ax, Ay, Az) = @. sqrt(Ax^2 + Ay^2 + Az^2)

function centraldiff(x, dims)
    ∇x = diff(x, dims = dims)

    if dims == 1
        a = cat(∇x[1:1, :, :], ∇x, dims = dims)
        a .+= cat(∇x, ∇x[end:end, :, :], dims = dims)
    elseif dims == 2
        a = cat(∇x[:, 1:1, :], ∇x, dims = dims)
        a .+= cat(∇x, ∇x[:, end:end, :], dims = dims)
    else
        a = cat(∇x[:, :, 1:1], ∇x, dims = dims)
        a .+= cat(∇x, ∇x[:, :, end:end], dims = dims)
    end

    return a
end

function centraldiff3D(cube)
    ∇x = centraldiff(cube, 1)
    ∇y = centraldiff(cube, 2)
    ∇z = centraldiff(cube, 3)

    return ∇x, ∇y, ∇z
end

function dot_product3D(Ax, Ay, Az, Bx, By, Bz)
    return Ax .* Bx .+ Ay .* By .+ Az .* Bz
end

function cross_product3D(Ax, Ay, Az, Bx, By, Bz)
    Cx = Ay .* Bz .- Az .* By
    Cy = Az .* Bx .- Ax .* Bz
    Cz = Ax .* By .- Ay .* Bx

    return Cx, Cy, Cz
end

function calculate_angletan(Bx, By, Bz, ∇n_x, ∇n_y, ∇n_z)
    norm∇n = vectornorm(∇n_x, ∇n_y, ∇n_z)
    normB = vectornorm(Bx, By, Bz)

    dot_product = dot_product3D(Bx, By, Bz, ∇n_x, ∇n_y, ∇n_z)

    Cx, Cy, Cz = cross_product3D(Bx, By, Bz, ∇n_x, ∇n_y, ∇n_z)
    normC = vectornorm(Cx, Cy, Cz)

    theta = atan.(normC, dot_product)
    cos_theta = cos.(theta)

    return cos_theta, theta
end

function calculate_anglecos(Ax, Ay, Az, Bx, By, Bz)
    norm = vectornorm(Ax, Ay, Az) .* vectornorm(Bx, By, Bz)
    dot_product = dot_product3D(Ax, Ay, Az, Bx, By, Bz)

    cos_theta = dot_product ./ norm
    theta = acos.(cos_theta)

    return cos_theta, theta
end
