using FFTW

# =============================================================================
# Application of a Fourier-domain transfer function
# =============================================================================

"""
    apply_instrument_2d(img, H)

Apply a 2D Fourier-domain instrumental transfer function. `H` is given in
FFT order (zero frequency at index [1,1]). It can be the binary mask of
`instrument_bandpass_L` or the smooth transfer function of
`instrument_bandpass_aperture`.
"""
function apply_instrument_2d(img::AbstractMatrix, H::AbstractMatrix)
    size(img) == size(H) || error("Filter shape mismatch: image size=$(size(img)) filter size=$(size(H))")
    all(isfinite, img) || throw(ArgumentError(
        "apply_instrument_2d requires finite image values; fill, crop, or inpaint masked pixels before the FFT."))
    all(isfinite, H) || throw(ArgumentError("The Fourier filter H contains non-finite values."))
    return real.(ifft(fft(img) .* H))
end

"""
    apply_to_array_xy(data, H; n=size(H, 1), m=size(H, 2))

Apply a Fourier-domain transfer function to a 2D image or to every sky-plane
slice of a Stokes cube. Supported cube layouts are `(n, m, nν)` and
`(nν, n, m)`.

`H` is either
  * a matrix, applied identically to every channel (achromatic filter), or
  * a function `ic -> H_ic` returning the matrix for channel index `ic`
    (chromatic filter, e.g. an aperture edge width kD ∝ ν). In that case
    `n` and `m` must be given explicitly.
"""
function apply_to_array_xy(data, H::AbstractMatrix; n::Int=size(H, 1), m::Int=size(H, 2))
    size(H) == (n, m) || error("Filter mask H must have size ($n,$m), got $(size(H))")
    return _apply_to_array_xy(data, _ -> H, n, m)
end

function apply_to_array_xy(data, Hfun::Function; n::Int, m::Int)
    return _apply_to_array_xy(data, Hfun, n, m)
end

function _apply_to_array_xy(data, Hfun, n::Int, m::Int)
    nd = ndims(data)

    if nd == 2
        size(data) == (n, m) || error("2D input must have size ($n,$m), got $(size(data))")
        return apply_instrument_2d(data, Hfun(1))
    elseif nd == 3
        sz = size(data)
        Tout = float(eltype(data))
        out = similar(data, Tout, sz)

        if sz[1] == n && sz[2] == m
            @views for k in axes(data, 3)
                out[:, :, k] = apply_instrument_2d(data[:, :, k], Hfun(k))
            end
            return out
        elseif sz[2] == n && sz[3] == m
            @views for k in axes(data, 1)
                out[k, :, :] = apply_instrument_2d(data[k, :, :], Hfun(k))
            end
            return out
        else
            error("Unsupported 3D shape $(sz). Expected (n,m,nν) or (nν,n,m) with n=$n m=$m.")
        end
    else
        error("Unsupported ndims(data)=$nd")
    end
end

# =============================================================================
# Shared input checks
# =============================================================================

function _check_bandpass_inputs(n, m, Δx, Δy, Lcut_small, Llarge, fNy)
    n > 0 || error("n must be positive, got $n")
    m > 0 || error("m must be positive, got $m")
    isfinite(Δx) && Δx > 0 || error("Δx must be positive and finite, got $Δx")
    isfinite(Δy) && Δy > 0 || error("Δy must be positive and finite, got $Δy")
    isfinite(Lcut_small) && Lcut_small > 0 || error("Lcut_small must be positive and finite, got $Lcut_small")
    isfinite(Llarge) && Llarge > 0 || error("Llarge must be positive and finite, got $Llarge")
    isfinite(fNy) && fNy > 0 || error("fNy must be positive and finite, got $fNy")
    Llarge > Lcut_small || error("Llarge ($Llarge) must be larger than Lcut_small ($Lcut_small)")
    return nothing
end

# =============================================================================
# 1) Hard 0/1 band-pass (original model, Eq. 8)
# =============================================================================

"""
    instrument_bandpass_L(n, m; Δx, Δy=Δx, Lcut_small, Llarge, fNy)

Build a hard 0/1 spatial-frequency band-pass mask. The mask removes scales
larger than `Llarge` and smaller than `Lcut_small`, capped at the Nyquist
frequency `fNy`. Returns `(H, fftshift(H))`, with `H` in FFT order.
"""
function instrument_bandpass_L(n::Int, m::Int;
                               Δx::Real, Δy::Real=Δx,
                               Lcut_small::Real,
                               Llarge::Real,
                               fNy::Real)
    _check_bandpass_inputs(n, m, Δx, Δy, Lcut_small, Llarge, fNy)

    # fftfreq(n, fs) expects the *sampling frequency* fs = 1/Δx, not the step Δx.
    # Spatial frequencies are then in cycles per unit of Δx (same unit as
    # Lcut_small/Llarge, which must be expressed in that same length unit).
    fx = FFTW.fftfreq(n, 1 / Δx)
    fy = FFTW.fftfreq(m, 1 / Δy)

    flo = 1 / Llarge
    fhi_raw = 1 / Lcut_small
    fhi = min(fhi_raw, fNy)

    @debug "Band-pass filter frequencies" Lcut_small=Lcut_small Llarge=Llarge flo=flo fhi=fhi fhi_raw=fhi_raw

    flo2 = flo^2
    fhi2 = fhi^2
    H = Matrix{Float32}(undef, n, m)
    @inbounds for j in 1:m
        fy2 = float(fy[j])^2
        for i in 1:n
            f2 = float(fx[i])^2 + fy2
            H[i, j] = (f2 >= flo2 && f2 <= fhi2) ? 1f0 : 0f0
        end
    end

    return H, fftshift(H)
end

# =============================================================================
# 2) Band-pass with edges smoothed by the antenna aperture autocorrelation
# =============================================================================
#
# A baseline b does not sample a single point of the uv plane but a patch
# given by the autocorrelation of the station aperture, of radius D/λ.
# The transfer function is therefore the annulus [flo, fhi] convolved with
# the normalised autocorrelation of a uniformly illuminated circular
# aperture ("Chinese hat"):
#
#     A(ρ) = (2/π) [acos(x) - x sqrt(1 - x²)],   x = ρ / kD ≤ 1,
#
# with kD = D/λ expressed in the same frequency units as flo, fhi.
#
# Because both the annulus and A are radially symmetric, the 2D convolution
# reduces to a 1D radial profile, computed here semi-analytically:
#
#     disk_R ∗ Â (k) = ∫ Â(ρ) ρ φ(k, ρ, R) dρ / ∫ Â(ρ) ρ dρ,
#
# where φ(k, ρ, R) is the fraction of the circle of radius ρ centred at
# distance k from the origin that lies inside the disk of radius R.
# This is exact (up to quadrature) and does NOT require kD to be resolved
# by the k grid: each pixel receives the true value of H at its |k|.
# =============================================================================

"""
    aperture_acf(ρ, kD)

Autocorrelation of a uniformly illuminated circular aperture, normalised to
1 at ρ = 0 and vanishing for ρ ≥ kD.
"""
function aperture_acf(ρ::Real, kD::Real)
    isfinite(kD) && kD > 0 || error("kD must be positive and finite")
    isfinite(ρ) && ρ >= 0 || error("ρ must be finite and nonnegative")
    ρ >= kD && return 0.0
    x = ρ / kD
    return (2 / π) * (acos(x) - x * sqrt(1 - x^2))
end

# Fraction of the circle of radius ρ, centred at distance k from the origin,
# lying inside the disk of radius R centred at the origin.
function _circle_in_disk_fraction(k::Float64, ρ::Float64, R::Float64)
    k + ρ <= R && return 1.0            # circle entirely inside
    (ρ >= k + R || k >= ρ + R) && return 0.0  # entirely outside (or disk inside circle)
    c = (k^2 + ρ^2 - R^2) / (2k * ρ)
    return acos(clamp(c, -1.0, 1.0)) / π
end

# Radial profile of (disk of radius R) ∗ (normalised aperture ACF of radius kD).
function _smoothed_disk(k::Float64, R::Float64, ρq::Vector{Float64}, wq::Vector{Float64})
    k + last(ρq) <= R && return 1.0
    k >= R + last(ρq) && return 0.0
    s = 0.0
    @inbounds for q in eachindex(ρq)
        s += wq[q] * _circle_in_disk_fraction(k, ρq[q], R)
    end
    return s
end

"""
    aperture_edge_profile(k, flo, fhi, kD; nquad=256, smooth_high=true)

Radial transfer function H(|k|) of an annulus [flo, fhi] convolved with the
autocorrelation of a circular aperture of radius `kD` (all in the same
frequency units). H ≈ 0.5 at flo, H = 1 on the plateau, and the transition
spans [flo - kD, flo + kD]. If `smooth_high=false`, only the inner edge is
smoothed and the outer edge stays sharp at `fhi`.
"""
function aperture_edge_profile(k::AbstractVector, flo::Real, fhi::Real, kD::Real;
                               nquad::Int=256, smooth_high::Bool=true)
    isfinite(kD) && kD > 0 || error("kD must be positive and finite")
    nquad > 0 || error("nquad must be positive")
    isfinite(flo) && isfinite(fhi) && 0 <= flo <= fhi || error("Require 0 <= flo <= fhi, both finite")
    all(x -> isfinite(x) && x >= 0, k) || error("Radial frequencies must be finite and nonnegative")
    kD = Float64(kD)
    # midpoint quadrature in ρ, weights ∝ A(ρ) ρ, normalised to unit sum
    ρq = [(q - 0.5) * kD / nquad for q in 1:nquad]
    wq = [aperture_acf(ρ, kD) * ρ for ρ in ρq]
    wq ./= sum(wq)

    flo_, fhi_, kD_ = Float64(flo), Float64(fhi), Float64(kD)
    prof = similar(k, Float64)
    @inbounds for i in eachindex(k)
        ki = Float64(k[i])
        outer = smooth_high ? _smoothed_disk(ki, fhi_, ρq, wq) : (ki <= fhi_ ? 1.0 : 0.0)
        inner = if ki <= flo_ - kD_
            1.0
        elseif ki >= flo_ + kD_
            0.0
        else
            _smoothed_disk(ki, flo_, ρq, wq)
        end
        prof[i] = clamp(outer - inner, 0.0, 1.0)
    end
    return prof
end

"""
    instrument_bandpass_aperture(n, m; Δx, Δy=Δx, Lcut_small, Llarge, fNy,
                                 kD, smooth_high=true, nquad=256, nlut=8192)

Smooth band-pass transfer function: the hard annulus of
`instrument_bandpass_L` convolved with the autocorrelation of a circular
station aperture of radius `kD` in the uv plane.

`kD` is in cycles per unit length (same unit as 1/Δx and 1/Llarge). For a
station of diameter D observed at wavelength λ and a source at distance d,
`kD = D / (λ d)`; see `aperture_kD`. For a fixed physical baseline cut the
ratio kD/flo = D/b_min; for a uv-cut in wavelengths only kD varies with ν.

The high cut is capped at the Nyquist frequency `fNy`. With
`smooth_high=true` the outer edge is also smoothed; the part of that edge
beyond `fNy` is simply not represented on the grid.

H is evaluated through a fine 1D look-up table in |k| (`nlut` points) with
linear interpolation. Returns `(H, fftshift(H))`, with `H` in FFT order.
"""
function instrument_bandpass_aperture(n::Int, m::Int;
                                      Δx::Real, Δy::Real=Δx,
                                      Lcut_small::Real,
                                      Llarge::Real,
                                      fNy::Real,
                                      kD::Real,
                                      smooth_high::Bool=true,
                                      nquad::Int=256,
                                      nlut::Int=8192)
    _check_bandpass_inputs(n, m, Δx, Δy, Lcut_small, Llarge, fNy)
    isfinite(kD) && kD > 0 || error("kD must be positive and finite, got $kD")

    nquad > 0 || error("nquad must be positive")
    nlut >= 2 || error("nlut must be at least 2")
    fx = FFTW.fftfreq(n, 1 / Δx)
    fy = FFTW.fftfreq(m, 1 / Δy)

    flo = 1 / Llarge
    fhi = min(1 / Lcut_small, fNy)
    flo > fhi && return zeros(Float32, n, m), zeros(Float32, n, m)
    kD > flo && @warn "kD ≥ flo (D ≥ b_min): the inner edge extends down to k = 0, so H(0) > 0 and part of the mean is kept" kD flo
    Δk = min(1 / (n * Δx), 1 / (m * Δy))
    kD < Δk && @debug "kD is smaller than the k-grid spacing; the edge is sampled pointwise but spans less than one Fourier cell" kD Δk

    @debug "Aperture band-pass" flo=flo fhi=fhi kD=kD Δk=Δk

    # look-up table in |k|
    kmax_grid = sqrt(maximum(abs2, fx) + maximum(abs2, fy))
    klut = collect(range(0.0, kmax_grid; length=nlut))
    # Keep a sharp outer cutoff out of the interpolated LUT so interpolation
    # cannot attenuate the last retained Fourier cell.
    profile_hi = smooth_high ? fhi : max(fhi, kmax_grid + kD)
    Hlut = aperture_edge_profile(klut, flo, profile_hi, kD; nquad=nquad, smooth_high=smooth_high)
    dk = klut[2] - klut[1]

    H = Matrix{Float32}(undef, n, m)
    @inbounds for j in 1:m
        fy2 = float(fy[j])^2
        for i in 1:n
            f = sqrt(float(fx[i])^2 + fy2)
            t = dk == 0 ? 0.0 : f / dk
            i0 = min(floor(Int, t), nlut - 2)
            w = t - i0
            h = (1 - w) * Hlut[i0 + 1] + w * Hlut[i0 + 2]
            # keep the Nyquist cap hard: frequencies above fNy are not sampled
            H[i, j] = f > fNy || (!smooth_high && f > fhi) ? 0f0 : Float32(h)
        end
    end

    return H, fftshift(H)
end

"""
    aperture_kD(ν; d, D=30.75, c=2.99792458e8)

Radius of the aperture autocorrelation in the uv plane, converted to
cycles per unit length at the source distance `d`: kD = D ν / (c d).
`ν` in Hz, `D` in metres, `d` in the length unit of the maps (e.g. pc);
the result is in cycles per that unit. Default `D = 30.75` m is the
effective Dutch HBA station diameter in LoTSS (HBA_DUAL_INNER).
"""
function aperture_kD(ν::Real; d::Real, D::Real=30.75, c::Real=2.99792458e8)
    all(x -> isfinite(x) && x > 0, (ν, d, D, c)) || error("ν, d, D and c must be positive and finite")
    return D * ν / (c * d)
end

# =============================================================================
# Single entry point
# =============================================================================

"""
    instrument_bandpass(n, m; edge=:hard, kwargs...)

Dispatch between the two instrument models:
  * `edge = :hard`     → `instrument_bandpass_L` (0/1 mask, Eq. 8)
  * `edge = :aperture` → `instrument_bandpass_aperture` (requires `kD`)
"""
function instrument_bandpass(n::Int, m::Int; edge::Symbol=:hard, kwargs...)
    edge === :hard     && return instrument_bandpass_L(n, m; kwargs...)
    edge === :aperture && return instrument_bandpass_aperture(n, m; kwargs...)
    error("Unknown edge model :$edge (use :hard or :aperture)")
end

# =============================================================================
# Usage examples
# =============================================================================
#
# # (a) original hard mask, achromatic
# H, _ = instrument_bandpass(n, m; edge=:hard, Δx=Δx, Lcut_small=1.0,
#                            Llarge=Llarge, fNy=fNy)
# Qf = apply_to_array_xy(Q, H)
#
# # (b) aperture-smoothed edges, single kD at the reference frequency 144 MHz
# kD = aperture_kD(144e6; d=d_pc)
# H, _ = instrument_bandpass(n, m; edge=:aperture, Δx=Δx, Lcut_small=1.0,
#                            Llarge=Llarge, fNy=fNy, kD=kD)
# Qf = apply_to_array_xy(Q, H)
#
# # (c) aperture-smoothed edges, chromatic kD(ν) (LoTSS: uv-cut in λ fixes
# #     flo, only kD ∝ ν changes). Cube layout (n, m, nν), frequencies νs in Hz.
# Hs = Dict{Int,Matrix{Float32}}()
# Hof(ic) = get!(Hs, ic) do
#     first(instrument_bandpass(n, m; edge=:aperture, Δx=Δx, Lcut_small=1.0,
#                               Llarge=Llarge, fNy=fNy,
#                               kD=aperture_kD(νs[ic]; d=d_pc)))
# end
# Qf = apply_to_array_xy(Q, Hof; n=n, m=m)
# Uf = apply_to_array_xy(U, Hof; n=n, m=m)   # reuses the cached filters


"""Normalize the pipeline's `filter` configuration. Spatial frequencies use cycles/pixel."""
function normalize_filter_options(options=nothing)
    options === nothing && (options = Dict{String, Any}())
    options isa AbstractDict || throw_config_error("`filter` must be an object."; code=:invalid_filter)
    defaults = Dict{String, Any}("edge" => "hard", "Lcut_small" => 2.0,
        "kD" => nothing, "chromatic" => false, "reference_frequency_mhz" => 144.0,
        "smooth_high" => true, "nquad" => 256, "nlut" => 8192)
    for (key, value) in options
        haskey(defaults, key) || throw_config_error("Unknown filter option: $key"; code=:invalid_filter)
        defaults[key] = value
    end
    defaults["edge"] = lowercase(string(defaults["edge"]))
    defaults["edge"] in ("hard", "aperture") || throw_config_error("filter.edge must be hard or aperture"; code=:invalid_filter)
    for key in ("Lcut_small", "reference_frequency_mhz", "kD")
        value = defaults[key]
        key == "kD" && value === nothing && continue
        value isa Real && !(value isa Bool) && isfinite(value) && value > 0 ||
            throw_config_error("filter.$key must be positive and finite"; code=:invalid_filter)
        defaults[key] = Float64(value)
    end
    for key in ("chromatic", "smooth_high")
        defaults[key] isa Bool || throw_config_error("filter.$key must be boolean"; code=:invalid_filter)
    end
    for (key, minimum) in (("nquad", 1), ("nlut", 2))
        value = defaults[key]
        value isa Integer && !(value isa Bool) && value >= minimum ||
            throw_config_error("filter.$key must be an integer >= $minimum"; code=:invalid_filter)
        defaults[key] = Int(value)
    end
    defaults["edge"] == "aperture" && defaults["kD"] === nothing &&
        throw_config_error("filter.kD (cycles/pixel) is required for aperture filtering"; code=:invalid_filter)
    defaults["edge"] == "hard" && defaults["chromatic"] &&
        throw_config_error("filter.chromatic requires edge=aperture"; code=:invalid_filter)
    return defaults
end
