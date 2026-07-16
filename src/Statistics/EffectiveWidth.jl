function EffectiveWidth(spectrum, xarray)
    length(spectrum) == length(xarray) ||
        throw(ArgumentError("spectrum and xarray must have the same length."))
    length(xarray) >= 2 || throw(ArgumentError("EffectiveWidth needs at least two samples."))
    all(isfinite, xarray) || throw(ArgumentError("xarray must contain only finite values."))

    finite_spectrum = filter(isfinite, spectrum)
    isempty(finite_spectrum) && return NaN
    dx = mean(diff(xarray))
    isfinite(dx) || return NaN
    maxspectrum = maximum(finite_spectrum)
    if maxspectrum == 0
        return 0.0
    end
    return sum(finite_spectrum) * dx / maxspectrum
end
