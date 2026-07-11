"""
RM-CLEAN deconvolution of Faraday dispersion functions (Heald 2009).

RM synthesis produces a *dirty* Faraday dispersion function (FDF) that is the
true Faraday spectrum convolved with the Rotation Measure Spread Function
(RMSF). Because the `λ²` sampling is incomplete, the RMSF has strong sidelobes
that contaminate the recovered FDF. RM-CLEAN removes them with a Högbom-style
loop: it iteratively locates the brightest pixel of the residual FDF, records a
fraction of it as a clean component, and subtracts the correspondingly shifted
and scaled RMSF. The clean components are then restored with a Gaussian beam
whose width matches the RMSF main lobe and added back to the residual.

The routines here build on [`RMSynthesis`](@ref) and
[`rmsf_diagnostics`](@ref) and mirror their array conventions (Faraday depth is
the last array dimension; leading dimensions are spatial).
"""

"""
    RMCleanResult

Result of an RM-CLEAN deconvolution.

# Fields
- `cleanFDF`: `|restored FDF|`, same spatial shape as the input with Faraday
  depth as the last axis.
- `realCleanFDF`, `imagCleanFDF`: real and imaginary parts of the restored FDF.
- `model`: complex clean-component model (delta functions on the `φ` grid).
- `residual`: complex residual FDF left after cleaning.
- `phi::Vector{Float64}`: Faraday-depth axis of the FDF (rad/m²).
- `rmsf::RMSFDiagnostics`: the RMSF and resolution metrics used for cleaning.
- `niter::Int`, `gain::Float64`, `threshold::Float64`: the loop settings used.
"""
struct RMCleanResult{A,M,R}
    cleanFDF::A
    realCleanFDF::A
    imagCleanFDF::A
    model::M
    residual::R
    phi::Vector{Float64}
    rmsf::RMSFDiagnostics
    niter::Int
    gain::Float64
    threshold::Float64
end

"""
    rmclean_1d(dirty, rmsf, rmsf_center, dphi, fwhm; gain, threshold, niter)
        -> (restored, model, residual, used_iter)

Deconvolve a single complex dirty FDF spectrum.

# Arguments
- `dirty::AbstractVector{<:Complex}`: dirty FDF sampled on the Faraday-depth grid.
- `rmsf::AbstractVector{<:Complex}`: complex RMSF sampled on a symmetric lag grid
  that shares `dirty`'s spacing.
- `rmsf_center::Int`: index of zero lag in `rmsf`.
- `dphi::Real`: Faraday-depth spacing (rad/m²).
- `fwhm::Real`: restoring-beam FWHM (rad/m²).

# Keywords
- `gain` (default `0.1`): loop gain (fraction of the peak subtracted per iteration).
- `threshold` (default `0.0`): stop once the residual peak `|F|` drops to or below
  this absolute level.
- `niter` (default `1000`): maximum number of clean iterations.
"""
function rmclean_1d(dirty::AbstractVector{<:Complex}, rmsf::AbstractVector{<:Complex},
                    rmsf_center::Integer, dphi::Real, fwhm::Real;
                    gain::Real = 0.1, threshold::Real = 0.0, niter::Integer = 1000)
    n = length(dirty)
    nrmsf = length(rmsf)
    residual = ComplexF64.(dirty)          # fresh copy; never mutates the input
    model = zeros(ComplexF64, n)

    used = 0
    for _ in 1:niter
        # Locate the brightest residual pixel.
        pmax = 1
        amax = abs(residual[1])
        @inbounds for k in 2:n
            ak = abs(residual[k])
            if ak > amax
                amax = ak
                pmax = k
            end
        end
        (isfinite(amax) && amax > threshold) || break

        comp = gain * residual[pmax]
        model[pmax] += comp

        # Subtract the shifted, scaled RMSF.
        @inbounds for i in 1:n
            ridx = rmsf_center + (i - pmax)
            (1 <= ridx <= nrmsf) || continue
            residual[i] -= comp * rmsf[ridx]
        end
        used += 1
    end

    # Restore with a Gaussian clean beam matched to the RMSF main lobe.
    sigma = fwhm / (2 * sqrt(2 * log(2)))
    restored = copy(residual)
    if isfinite(sigma) && sigma > 0
        half = max(1, ceil(Int, 3 * sigma / dphi))
        beam = [exp(-((k * dphi)^2) / (2 * sigma^2)) for k in -half:half]
        @inbounds for p in 1:n
            m = model[p]
            m == 0 && continue
            for (j, k) in enumerate(-half:half)
                i = p + k
                (1 <= i <= n) || continue
                restored[i] += m * beam[j]
            end
        end
    else
        restored .+= model
    end

    return restored, model, residual, used
end

# Reshape an FDF array whose last dimension is Faraday depth into an
# (npix, nPhi) matrix and remember the original shape for reconstruction.
function _fdf_to_matrix(F::AbstractArray)
    sz = size(F)
    nPhi = sz[end]
    npix = length(F) ÷ nPhi            # 1 for a plain vector
    return reshape(F, npix, nPhi), sz
end

"""
    rmclean(realFDF, imagFDF, PhiArray, diag; gain, threshold, niter, log_progress)
        -> RMCleanResult

Run RM-CLEAN on a *precomputed* dirty FDF (its real and imaginary parts, as
returned by [`RMSynthesis`](@ref)) using the RMSF in `diag`
([`RMSFDiagnostics`](@ref)). `realFDF`/`imagFDF` may be 1D, 2D, or 3D with
Faraday depth as the last axis.
"""
function rmclean(realFDF::AbstractArray, imagFDF::AbstractArray, PhiArray::AbstractArray,
                 diag::RMSFDiagnostics; gain::Real = 0.1, threshold::Real = 0.0,
                 niter::Integer = 1000, log_progress::Bool = false)
    size(realFDF) == size(imagFDF) || error("realFDF and imagFDF must have the same shape.")
    size(realFDF)[end] == length(PhiArray) ||
        error("The last FDF dimension ($(size(realFDF)[end])) must match length(PhiArray) ($(length(PhiArray))).")

    dphi = _uniform_step(PhiArray; name = "PhiArray")
    rmsf_center = argmin(abs.(diag.phi))

    realMat, sz = _fdf_to_matrix(realFDF)
    imagMat, _ = _fdf_to_matrix(imagFDF)
    npix, nPhi = size(realMat)

    log_progress && @info "Starting RM-CLEAN" n_pixels = npix n_phi = nPhi gain = gain niter = niter

    restoredMat = Matrix{ComplexF64}(undef, npix, nPhi)
    modelMat = Matrix{ComplexF64}(undef, npix, nPhi)
    residualMat = Matrix{ComplexF64}(undef, npix, nPhi)

    spectrum = Vector{ComplexF64}(undef, nPhi)
    nan_spectrum = ComplexF64(NaN, NaN)
    for p in 1:npix
        # Masked pixels (NaN FDF from HEALPix UNSEEN or partial-sky maps):
        # skip the clean loop and propagate NaN.
        masked = false
        @inbounds for k in 1:nPhi
            re, im_ = realMat[p, k], imagMat[p, k]
            if !(isfinite(re) && isfinite(im_))
                masked = true
                break
            end
            spectrum[k] = ComplexF64(re, im_)
        end
        if masked
            @inbounds for k in 1:nPhi
                restoredMat[p, k] = nan_spectrum
                modelMat[p, k] = nan_spectrum
                residualMat[p, k] = nan_spectrum
            end
            log_progress && print_progress(p, npix; label="RM-CLEAN")
            continue
        end
        restored, model, residual, _ = rmclean_1d(spectrum, diag.rmsf, rmsf_center, dphi, diag.fwhm;
                                                  gain = gain, threshold = threshold, niter = niter)
        @inbounds for k in 1:nPhi
            restoredMat[p, k] = restored[k]
            modelMat[p, k] = model[k]
            residualMat[p, k] = residual[k]
        end
        log_progress && print_progress(p, npix; label="RM-CLEAN")
    end

    restored = reshape(restoredMat, sz)
    model = reshape(modelMat, sz)
    residual = reshape(residualMat, sz)

    log_progress && @info "RM-CLEAN complete" output_size = size(restored)

    return RMCleanResult(
        abs.(restored),
        real.(restored),
        imag.(restored),
        model,
        residual,
        Float64.(collect(PhiArray)),
        diag,
        Int(niter),
        Float64(gain),
        Float64(threshold),
    )
end

"""
    RMClean(Q, U, nuArray, PhiArray; gain, threshold, niter, diagnostics, log_progress)
        -> RMCleanResult

Perform RM synthesis on Stokes `Q`/`U` and deconvolve the resulting Faraday
dispersion function with RM-CLEAN.

# Arguments
- `Q`, `U`: Stokes cubes with frequency as the last axis (1D/2D/3D), matching
  [`RMSynthesis`](@ref).
- `nuArray::AbstractArray`: frequencies in **Hz**.
- `PhiArray::AbstractArray`: uniformly spaced Faraday depths (rad/m²).

# Keywords
- `gain` (default `0.1`), `threshold` (default `0.0`), `niter` (default `1000`):
  RM-CLEAN loop settings.
- `diagnostics`: a precomputed [`RMSFDiagnostics`](@ref) to reuse; computed from
  `nuArray`/`PhiArray` when `nothing`.
- `log_progress` (default `false`): emit progress logging.

# Example
```julia
result = RMClean(Q, U, nuArray_Hz, PhiArray; gain = 0.1, niter = 2000)
result.cleanFDF        # |restored FDF|
result.rmsf.fwhm       # restoring-beam FWHM (rad/m²)
```
"""
function RMClean(Q::AbstractArray, U::AbstractArray, nuArray::AbstractArray, PhiArray::AbstractArray;
                 gain::Real = 0.1, threshold::Real = 0.0, niter::Integer = 1000,
                 diagnostics::Union{Nothing,RMSFDiagnostics} = nothing, log_progress::Bool = false)
    _, realFDF, imagFDF = RMSynthesis(Q, U, nuArray, PhiArray; log_progress = log_progress)
    diag = diagnostics === nothing ? rmsf_diagnostics(nuArray, PhiArray) : diagnostics
    return rmclean(realFDF, imagFDF, PhiArray, diag;
                   gain = gain, threshold = threshold, niter = niter, log_progress = log_progress)
end

"""
    write_rmsf(resultspath, diag; filename="RMSF.fits", ensure_path=true, metadata=nothing)
        -> String

Write the RMSF and its resolution metrics to a FITS file. The image has shape
`(nφ, 3)` holding `|R|`, `Re R`, and `Im R`; the Faraday-depth axis is described
by `CRVAL1`/`CDELT1` and the scalar metrics are stored in the header keywords
`RMSFFWHM`, `RMSFTHEO`, `PHIMAX`, `MAXSCALE`, and `LAMBDA02`.
"""
function write_rmsf(resultspath::AbstractString, diag::RMSFDiagnostics;
                    filename::AbstractString = "RMSF.fits", ensure_path::Bool = true, metadata = nothing)
    ensure_path && mkpath(resultspath)

    nphi = length(diag.phi)
    data = Matrix{Float64}(undef, nphi, 3)
    @inbounds for k in 1:nphi
        data[k, 1] = abs(diag.rmsf[k])
        data[k, 2] = real(diag.rmsf[k])
        data[k, 3] = imag(diag.rmsf[k])
    end

    header = FITSHeader(["NAXIS"], [2], [""])
    header["NAXIS1"] = nphi
    header["NAXIS2"] = 3
    header["CTYPE1"] = "FARADAY"
    header["CRVAL1"] = diag.phi[1]
    header["CRPIX1"] = 1
    header["CDELT1"] = nphi >= 2 ? diag.phi[2] - diag.phi[1] : 0.0
    header["CUNIT1"] = "rad/m^2"
    header["CTYPE2"] = "RMSFCOMP"
    header["BUNIT"] = ""
    header["RMSFFWHM"] = diag.fwhm
    header["RMSFTHEO"] = diag.fwhm_theoretical
    header["PHIMAX"] = diag.phi_max
    header["MAXSCALE"] = diag.max_scale
    header["LAMBDA02"] = diag.lambda0_2
    annotate_header!(header, metadata)

    fits_path = joinpath(resultspath, filename)
    atomic_write_path(fits_path) do tmp_path
        FITS(tmp_path, "w") do f
            write(f, data; header = header)
        end
    end

    @info "Wrote FITS file" data = "RMSF" path = fits_path
    return fits_path
end

"""
    RMCleanHealpix(Q, U, nuArray, PhiArray; nside, order, gain, threshold, niter,
                   diagnostics, log_progress) -> HealpixRMResult

HEALPix-aware RM-CLEAN. Accepts `Q`/`U` as `HealpixStack`s, `Npix × Nfreq`
matrices, or vectors of `Healpix.jl` maps (mirroring
[`RMSynthesisHealpix`](@ref)), runs [`RMClean`](@ref), and returns a
[`HealpixRMResult`](@ref) holding the **restored** Faraday dispersion function,
ready for [`write_healpix_rm_result`](@ref).
"""
function RMCleanHealpix(Q, U, nuArray::AbstractArray, PhiArray::AbstractArray;
                        nside::Union{Nothing,Integer} = nothing,
                        order::HealpixOrderName = :ring,
                        gain::Real = 0.1, threshold::Real = 0.0, niter::Integer = 1000,
                        diagnostics::Union{Nothing,RMSFDiagnostics} = nothing,
                        log_progress::Bool = false)
    q_stack = _healpix_stack_from_input(Q; nside = nside, order = order)
    u_stack = _healpix_stack_from_input(U; nside = q_stack.nside, order = q_stack.order)

    _check_healpix_stack_consistency(q_stack, u_stack)
    length(nuArray) == size(q_stack.pixels, 2) ||
        error("nuArray length ($(length(nuArray))) must match the number of HEALPix frequency maps ($(size(q_stack.pixels, 2))).")

    result = RMClean(q_stack.pixels, u_stack.pixels, nuArray, PhiArray;
                     gain = gain, threshold = threshold, niter = niter,
                     diagnostics = diagnostics, log_progress = log_progress)

    return HealpixRMResult(
        _as_matrix(result.cleanFDF),
        _as_matrix(result.realCleanFDF),
        _as_matrix(result.imagCleanFDF),
        Float64.(collect(PhiArray)),
        q_stack.nside,
        q_stack.order,
        q_stack.coordsys,
    )
end

"""
    RMCleanAuto(Q, U, nuArray, PhiArray; kwargs...)

Run RM-CLEAN while automatically dispatching regular FITS images/cubes to
[`RMClean`](@ref) and HEALPix FITS maps/stacks to [`RMCleanHealpix`](@ref).
`Q` and `U` accept the same inputs as [`RMSynthesisAuto`](@ref).
"""
function RMCleanAuto(Q, U, nuArray::AbstractArray, PhiArray::AbstractArray;
                     conversion::Real = 1.0,
                     column = 1,
                     T::Type = Float64,
                     nside::Union{Nothing,Integer} = nothing,
                     order::HealpixOrderName = :ring,
                     gain::Real = 0.1, threshold::Real = 0.0, niter::Integer = 1000,
                     diagnostics::Union{Nothing,RMSFDiagnostics} = nothing,
                     log_progress::Bool = false,
                     allow_nonfinite::Bool = false)
    q_input = _read_auto_grid_input(Q; conversion = conversion, column = column, T = T, allow_nonfinite = allow_nonfinite)
    u_input = _read_auto_grid_input(U; conversion = conversion, column = column, T = T, allow_nonfinite = allow_nonfinite)
    q_is_healpix = _is_healpix_input(q_input; nside = nside)
    u_is_healpix = _is_healpix_input(u_input; nside = nside)

    q_is_healpix == u_is_healpix ||
        error("Q and U must both be HEALPix grids or both be regular image/cube grids.")

    if q_is_healpix
        return RMCleanHealpix(q_input, u_input, nuArray, PhiArray;
            nside = nside,
            order = order,
            gain = gain,
            threshold = threshold,
            niter = niter,
            diagnostics = diagnostics,
            log_progress = log_progress,
        )
    end

    return RMClean(q_input, u_input, nuArray, PhiArray;
        gain = gain,
        threshold = threshold,
        niter = niter,
        diagnostics = diagnostics,
        log_progress = log_progress,
    )
end
