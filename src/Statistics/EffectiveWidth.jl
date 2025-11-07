function EffectiveWidth(spectrum, xarray)
    dx = mean(diff(xarray))
    maxspectrum = maximum(spectrum)
    if maxspectrum == 0
        return 0.0
    end
    return sum(spectrum) * dx / maxspectrum
end