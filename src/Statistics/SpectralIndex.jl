"""
    spectral_index_map(cube::AbstractArray{<:Real,3}, nuArray::AbstractVector; min_channels::Integer=3)
        -> (index::Matrix{Float64}, index_err::Matrix{Float64})

Compute a per-pixel spectral index map from a multi-frequency intensity cube
by fitting `log10(I) = index * log10(nu) + const` with ordinary least squares
along the spectral (third) axis.

# Arguments
- `cube`: intensity cube of shape `(nx, ny, nchan)`. Any intensity-like
  quantity works (brightness temperature, flux density, ...); the fitted
  slope is the power-law index of that quantity versus frequency.
- `nuArray`: frequency of each channel. Any single consistent unit is fine
  (MHz, Hz, ...): the log-log slope is invariant under a rescaling of the
  frequency axis. All values must be positive and finite.

# Keywords
- `min_channels`: minimum number of valid channels required to fit a pixel
  (default `3`). Pixels with fewer valid channels are set to `NaN`.

# Returns
A tuple `(index, index_err)` of `(nx, ny)` matrices:
- `index`: fitted power-law index per pixel.
- `index_err`: 1σ standard error of the slope from the fit residuals.
  `NaN` when only two channels are used (zero degrees of freedom).

Channels with non-positive or non-finite intensities are excluded pixel by
pixel, so noisy channels that scatter below zero do not poison the fit.
Pixels where fewer than `min_channels` channels remain are set to `NaN` in
both outputs.

Note on conventions: fitting a brightness-temperature cube `T_nu` yields the
temperature index `beta` with `T ∝ nu^beta`; the flux-density index is
`alpha = beta + 2` (`S_nu ∝ nu^alpha`).
"""
function spectral_index_map(cube::AbstractArray{<:Real,3}, nuArray::AbstractVector; min_channels::Integer=3)
    nchan = size(cube, 3)
    length(nuArray) == nchan ||
        error("Frequency axis mismatch: cube has $(nchan) channels but nuArray has $(length(nuArray)) entries.")
    min_channels >= 2 ||
        error("min_channels must be >= 2 to fit a slope, got $(min_channels).")
    nchan >= min_channels ||
        error("Cube has $(nchan) channels but at least $(min_channels) are required (min_channels).")
    all(nu -> isfinite(nu) && nu > 0, nuArray) ||
        error("All frequencies must be positive and finite to take log10(nu).")

    lognu = log10.(Float64.(collect(nuArray)))

    nx, ny = size(cube, 1), size(cube, 2)
    index = fill(NaN, nx, ny)
    index_err = fill(NaN, nx, ny)

    Threads.@threads for idx in CartesianIndices((1:nx, 1:ny))
        i, j = Tuple(idx)
        slope, err = _loglog_slope(view(cube, i, j, :), lognu, min_channels)
        index[i, j] = slope
        index_err[i, j] = err
    end

    return index, index_err
end

"""
    _loglog_slope(values, lognu, min_channels) -> (slope, stderr)

Least-squares slope of `log10(values)` versus `lognu`, using only entries
where `values` is positive and finite. Returns `(NaN, NaN)` when fewer than
`min_channels` valid entries remain, and `stderr = NaN` when the fit has no
residual degrees of freedom (exactly two points).
"""
function _loglog_slope(values::AbstractVector, lognu::AbstractVector, min_channels::Integer)
    # Two-pass centered computation: log10 of small emissivities carries a
    # large offset (e.g. y ≈ -40), and one-pass sums of squares would lose the
    # centered variance to catastrophic cancellation.
    n = 0
    sx = 0.0
    sy = 0.0

    @inbounds for k in eachindex(values, lognu)
        v = values[k]
        (isfinite(v) && v > 0) || continue
        n += 1
        sx += lognu[k]
        sy += log10(Float64(v))
    end

    n >= min_channels || return (NaN, NaN)

    mx = sx / n
    my = sy / n
    sxx = 0.0
    sxy = 0.0
    syy = 0.0

    @inbounds for k in eachindex(values, lognu)
        v = values[k]
        (isfinite(v) && v > 0) || continue
        dx = lognu[k] - mx
        dy = log10(Float64(v)) - my
        sxx += dx * dx
        sxy += dx * dy
        syy += dy * dy
    end

    sxx > 0 || return (NaN, NaN)

    slope = sxy / sxx

    n > 2 || return (slope, NaN)

    # Residual sum of squares of the fit, clamped at zero against round-off.
    ssr = max(syy - slope * sxy, 0.0)
    stderr = sqrt(ssr / (n - 2) / sxx)

    return (slope, stderr)
end
