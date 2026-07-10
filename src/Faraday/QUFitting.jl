"""
QU fitting: direct model fitting of Stokes q(λ²) and u(λ²) spectra.

While RM synthesis ([`RMSynthesis`](@ref)) recovers the Faraday dispersion
function non-parametrically, QU fitting adjusts explicit physical models of
the Faraday rotation and depolarization directly to the observed `q = Q/I`
and `u = U/I` spectra as a function of `λ²` (O'Sullivan et al. 2012;
Sokoloff et al. 1998). Both approaches are complementary: QU fitting does not
suffer from RMSF sidelobes and yields uncertainties on physical parameters,
but requires choosing a model. Model selection between the built-in models is
provided through the Akaike and Bayesian information criteria.

Available models (`QU_FIT_MODELS`):
- `:screen` — external Faraday screen (no depolarization):
  `p(λ²) = p₀ · exp(2i(χ₀ + RM·λ²))`; parameters `[p0, chi0, RM]`.
- `:burn_slab` — uniform slab with internal Faraday rotation (Burn 1966):
  `p(λ²) = p₀ · sin(φλ²)/(φλ²) · exp(2i(χ₀ + φλ²/2))`;
  parameters `[p0, chi0, phi]`.
- `:external_dispersion` — external screen with turbulent RM dispersion
  (Burn 1966): `p(λ²) = p₀ · exp(−2σ²_RM λ⁴) · exp(2i(χ₀ + RM·λ²))`;
  parameters `[p0, chi0, RM, sigma_rm]`.
- `:internal_dispersion` — internal Faraday dispersion (Sokoloff et al. 1998):
  `p(λ²) = p₀ · e^{2iχ₀} · (1 − e^{−S})/S` with `S = 2σ²_RM λ⁴ − 2iφλ²`;
  parameters `[p0, chi0, phi, sigma_rm]`.

Angles are in radians, `RM`/`phi`/`sigma_rm` in rad/m². The polarization can
be fractional (`Q/I`, `U/I`) or in any consistent flux unit; `p0` is returned
in the same unit as the input.
"""

const QU_FIT_MODELS = (:screen, :burn_slab, :external_dispersion, :internal_dispersion)

const QU_FIT_PARAM_NAMES = Dict(
    :screen => ["p0", "chi0", "RM"],
    :burn_slab => ["p0", "chi0", "phi"],
    :external_dispersion => ["p0", "chi0", "RM", "sigma_rm"],
    :internal_dispersion => ["p0", "chi0", "phi", "sigma_rm"],
)

"""
    qu_model(model::Symbol, params::AbstractVector, lambda2::AbstractVector) -> Vector{ComplexF64}

Evaluate the complex polarization `p(λ²) = q + iu` of a QU-fitting `model`
(one of `QU_FIT_MODELS`) for parameter vector `params` on the `λ²` grid
`lambda2` (m²). See the module docstring of `QUFitting.jl` for the parameter
conventions of each model.
"""
function qu_model(model::Symbol, params::AbstractVector, lambda2::AbstractVector)
    n = length(lambda2)
    p = Vector{ComplexF64}(undef, n)

    if model === :screen
        p0, chi0, rm = params
        @inbounds for k in 1:n
            p[k] = p0 * cis(2 * (chi0 + rm * lambda2[k]))
        end
    elseif model === :burn_slab
        p0, chi0, phi = params
        @inbounds for k in 1:n
            x = phi * lambda2[k]
            depol = abs(x) < 1e-9 ? 1.0 : sin(x) / x
            p[k] = p0 * depol * cis(2 * chi0 + x)
        end
    elseif model === :external_dispersion
        p0, chi0, rm, sigma = params
        @inbounds for k in 1:n
            l2 = lambda2[k]
            p[k] = p0 * exp(-2 * sigma^2 * l2^2) * cis(2 * (chi0 + rm * l2))
        end
    elseif model === :internal_dispersion
        p0, chi0, phi, sigma = params
        @inbounds for k in 1:n
            l2 = lambda2[k]
            S = complex(2 * sigma^2 * l2^2, -2 * phi * l2)
            frac = abs(S) < 1e-9 ? complex(1.0) : (1 - exp(-S)) / S
            p[k] = p0 * cis(2 * chi0) * frac
        end
    else
        throw(ArgumentError("Unknown QU-fitting model :$model. " *
                            "Available models: $(join(QU_FIT_MODELS, ", "))."))
    end

    return p
end

"""
    QUFitResult

Result of a QU model fit (see [`QUFit`](@ref)).

# Fields
- `model::Symbol`: the fitted model (one of `QU_FIT_MODELS`).
- `param_names::Vector{String}`: names of the fitted parameters.
- `params::Vector{Float64}`: best-fit parameter values.
- `stderr::Vector{Float64}`: 1σ uncertainties from the fit covariance
  (`NaN` when the normal-equation matrix is singular).
- `chi2::Float64`, `dof::Int`, `chi2_red::Float64`: goodness of fit.
- `aic::Float64`, `bic::Float64`: information criteria for model selection.
- `converged::Bool`, `niter::Int`: optimizer status.
- `lambda2::Vector{Float64}`: `λ²` grid of the fit (m²).
- `q_model::Vector{Float64}`, `u_model::Vector{Float64}`: best-fit spectra.
"""
struct QUFitResult
    model::Symbol
    param_names::Vector{String}
    params::Vector{Float64}
    stderr::Vector{Float64}
    chi2::Float64
    dof::Int
    chi2_red::Float64
    aic::Float64
    bic::Float64
    converged::Bool
    niter::Int
    lambda2::Vector{Float64}
    q_model::Vector{Float64}
    u_model::Vector{Float64}
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

function _qu_residuals!(r::Vector{Float64}, model::Symbol, params::AbstractVector,
                        lambda2::Vector{Float64}, q::Vector{Float64}, u::Vector{Float64},
                        wq::Vector{Float64}, wu::Vector{Float64})
    p = qu_model(model, params, lambda2)
    n = length(lambda2)
    @inbounds for k in 1:n
        r[k] = (real(p[k]) - q[k]) * wq[k]
        r[n+k] = (imag(p[k]) - u[k]) * wu[k]
    end
    return r
end

function _qu_jacobian!(J::Matrix{Float64}, r0::Vector{Float64}, model::Symbol,
                       params::Vector{Float64}, lambda2, q, u, wq, wu)
    np = length(params)
    rp = similar(r0)
    trial = copy(params)
    for j in 1:np
        h = max(abs(params[j]), 1.0) * 1e-7
        trial[j] = params[j] + h
        _qu_residuals!(rp, model, trial, lambda2, q, u, wq, wu)
        @inbounds @. J[:, j] = (rp - r0) / h
        trial[j] = params[j]
    end
    return J
end

# Coarse grid search for the rotation measure that maximizes the derotated
# polarization — equivalent to locating the peak of the FDF, and a robust
# initial guess even for wrapped polarization angles.
function _initial_rm_chi0_p0(lambda2::Vector{Float64}, q::Vector{Float64}, u::Vector{Float64})
    p = complex.(q, u)
    dl2 = diff(sort(lambda2))
    dl2min = maximum([minimum(dl2[dl2 .> 0]; init=Inf), eps()])
    rm_max = isfinite(dl2min) ? pi / (2 * dl2min) : 1000.0
    best_rm = 0.0
    best_amp = -Inf
    best_phase = 0.0
    for rm in range(-rm_max, rm_max; length=2001)
        f = zero(ComplexF64)
        @inbounds for k in eachindex(lambda2)
            f += p[k] * cis(-2 * rm * lambda2[k])
        end
        a = abs(f)
        if a > best_amp
            best_amp = a
            best_rm = rm
            best_phase = angle(f)
        end
    end
    p0 = best_amp / length(lambda2)
    chi0 = best_phase / 2
    return best_rm, chi0, p0
end

function _initial_params(model::Symbol, lambda2, q, u)
    rm, chi0, p0 = _initial_rm_chi0_p0(lambda2, q, u)
    p0 = max(p0, 1e-12)
    sigma0 = max(abs(rm) / 4, 0.1 / max(maximum(lambda2), eps()))
    model === :screen && return [p0, chi0, rm]
    model === :burn_slab && return [p0, chi0, 2 * rm]
    model === :external_dispersion && return [p0, chi0, rm, sigma0]
    model === :internal_dispersion && return [p0, chi0, 2 * rm, sigma0]
    throw(ArgumentError("Unknown QU-fitting model :$model. " *
                        "Available models: $(join(QU_FIT_MODELS, ", "))."))
end

# Map parameters back to their canonical ranges: p0 ≥ 0, sigma_rm ≥ 0 (the
# models only depend on σ²) and χ0 ∈ (−π/2, π/2] (χ0 enters as 2χ0).
function _canonicalize!(model::Symbol, params::Vector{Float64})
    if params[1] < 0
        params[1] = -params[1]
        params[2] += pi / 2
    end
    if model === :external_dispersion || model === :internal_dispersion
        params[4] = abs(params[4])
    end
    params[2] = mod(params[2] + pi / 2, pi) - pi / 2
    return params
end

"""
    QUFit(q, u, nuArray; model=:screen, sigma_q=nothing, sigma_u=nothing,
          initial_params=nothing, max_iterations=200, tolerance=1e-10) -> QUFitResult

Fit a Faraday rotation/depolarization `model` (one of `QU_FIT_MODELS`) to the
Stokes spectra `q(ν)` and `u(ν)` by Levenberg–Marquardt least squares on the
joint `[q; u]` residuals.

# Arguments
- `q::AbstractVector`, `u::AbstractVector`: polarization spectra, typically
  fractional (`Q/I`, `U/I`) but any consistent unit works.
- `nuArray::AbstractVector`: channel frequencies in Hz.

# Keywords
- `model::Symbol`: one of `:screen`, `:burn_slab`, `:external_dispersion`,
  `:internal_dispersion` (default `:screen`).
- `sigma_q`, `sigma_u`: per-channel 1σ uncertainties (scalar or vector). When
  omitted the fit is unweighted and `chi2` is in spectrum units; `stderr`
  is then rescaled by `√chi2_red` as customary.
- `initial_params`: optional starting parameter vector; by default a robust
  guess is derived from a coarse RM grid search (FDF peak).
- `max_iterations`, `tolerance`: optimizer settings.

Channels where `q`, `u` or the uncertainties are not finite are ignored.

See also [`QUFitCompare`](@ref) to fit and rank several models, and
[`QUFitCube`](@ref) for per-pixel maps.
"""
function QUFit(q::AbstractVector, u::AbstractVector, nuArray::AbstractVector;
               model::Symbol = :screen,
               sigma_q = nothing, sigma_u = nothing,
               initial_params = nothing,
               max_iterations::Int = 200,
               tolerance::Float64 = 1e-10)

    model in QU_FIT_MODELS ||
        throw(ArgumentError("Unknown QU-fitting model :$model. " *
                            "Available models: $(join(QU_FIT_MODELS, ", "))."))
    length(q) == length(u) == length(nuArray) ||
        throw(ArgumentError("q, u and nuArray must have the same length " *
                            "(got $(length(q)), $(length(u)), $(length(nuArray)))."))

    lambda2_all = @. (C_m / float(nuArray))^2
    sq = sigma_q === nothing ? fill(NaN, length(q)) :
         (sigma_q isa Number ? fill(float(sigma_q), length(q)) : float.(sigma_q))
    su = sigma_u === nothing ? fill(NaN, length(u)) :
         (sigma_u isa Number ? fill(float(sigma_u), length(u)) : float.(sigma_u))

    weighted = sigma_q !== nothing || sigma_u !== nothing
    keep = [isfinite(q[k]) && isfinite(u[k]) && isfinite(lambda2_all[k]) &&
            (!weighted || (isfinite(sq[k]) && sq[k] > 0 && isfinite(su[k]) && su[k] > 0))
            for k in eachindex(nuArray)]

    lambda2 = Float64.(lambda2_all[keep])
    qv = Float64.(q[keep])
    uv = Float64.(u[keep])
    np = length(QU_FIT_PARAM_NAMES[model])
    nchan = length(lambda2)
    2 * nchan > np ||
        throw(ArgumentError("Not enough valid channels ($nchan) to fit model :$model ($np parameters)."))

    wq = weighted ? 1.0 ./ sq[keep] : ones(nchan)
    wu = weighted ? 1.0 ./ su[keep] : ones(nchan)

    params = initial_params === nothing ? _initial_params(model, lambda2, qv, uv) :
             Float64.(collect(initial_params))
    length(params) == np ||
        throw(ArgumentError("initial_params must have $np elements for model :$model."))

    # --- Levenberg–Marquardt ---
    m = 2 * nchan
    r = Vector{Float64}(undef, m)
    J = Matrix{Float64}(undef, m, np)
    _qu_residuals!(r, model, params, lambda2, qv, uv, wq, wu)
    chi2 = sum(abs2, r)

    damping = 1e-3
    converged = false
    niter = 0
    for iter in 1:max_iterations
        niter = iter
        _qu_jacobian!(J, r, model, params, lambda2, qv, uv, wq, wu)
        JtJ = transpose(J) * J
        g = transpose(J) * r

        accepted = false
        for _ in 1:12
            A = JtJ + damping * Diagonal(max.(diag(JtJ), 1e-12))
            delta = try
                -(A \ g)
            catch
                break
            end
            all(isfinite, delta) || break
            trial = params .+ delta
            rt = similar(r)
            _qu_residuals!(rt, model, trial, lambda2, qv, uv, wq, wu)
            chi2_trial = sum(abs2, rt)
            if isfinite(chi2_trial) && chi2_trial < chi2
                if chi2 - chi2_trial < tolerance * max(chi2, 1e-30)
                    converged = true
                end
                params = trial
                r .= rt
                chi2 = chi2_trial
                damping = max(damping / 3, 1e-12)
                accepted = true
                break
            else
                damping *= 10
            end
        end

        if !accepted
            converged = true  # no descent direction improves the fit further
        end
        converged && break
    end

    _canonicalize!(model, params)
    _qu_residuals!(r, model, params, lambda2, qv, uv, wq, wu)
    chi2 = sum(abs2, r)
    dof = m - np
    chi2_red = dof > 0 ? chi2 / dof : NaN

    # Parameter covariance from the final Jacobian.
    _qu_jacobian!(J, r, model, params, lambda2, qv, uv, wq, wu)
    stderr = fill(NaN, np)
    try
        cov = inv(transpose(J) * J)
        scale = weighted ? 1.0 : (dof > 0 ? chi2 / dof : NaN)
        for j in 1:np
            v = cov[j, j] * scale
            stderr[j] = v >= 0 ? sqrt(v) : NaN
        end
    catch
    end

    aic = chi2 + 2 * np
    bic = chi2 + np * log(m)

    pbest = qu_model(model, params, lambda2)
    return QUFitResult(model, copy(QU_FIT_PARAM_NAMES[model]), params, stderr,
                       chi2, dof, chi2_red, aic, bic, converged, niter,
                       lambda2, real.(pbest), imag.(pbest))
end

"""
    QUFitCompare(q, u, nuArray; models=QU_FIT_MODELS, kwargs...)
        -> (best::QUFitResult, results::Dict{Symbol,QUFitResult})

Fit each model in `models` with [`QUFit`](@ref) and return the model with the
lowest Bayesian information criterion together with all individual results.
Models that fail to fit are skipped (a warning is logged); at least one model
must succeed. Keyword arguments are forwarded to `QUFit`.
"""
function QUFitCompare(q::AbstractVector, u::AbstractVector, nuArray::AbstractVector;
                      models = QU_FIT_MODELS, kwargs...)
    results = Dict{Symbol,QUFitResult}()
    for model in models
        try
            results[model] = QUFit(q, u, nuArray; model=model, kwargs...)
        catch err
            err isa InterruptException && rethrow()
            @warn "QU fit failed for model, skipping" model err
        end
    end
    isempty(results) && error("QUFitCompare: all QU model fits failed.")
    best = argmin(m -> results[m].bic, collect(keys(results)))
    return results[best], results
end

"""
    QUFitCube(Q, U, nuArray; model=:screen, sigma_q=nothing, sigma_u=nothing,
              log_progress=false, kwargs...)
        -> (params::Array{Float64,3}, stderr::Array{Float64,3}, chi2_red::Matrix{Float64})

Run [`QUFit`](@ref) independently on every spatial pixel of the Stokes cubes
`Q` and `U` (shape `(nx, ny, nchan)`, frequency last, matching
[`RMSynthesis`](@ref) conventions).

Returns parameter maps `params[nx, ny, np]` and `stderr[nx, ny, np]` (ordered
as `QU_FIT_PARAM_NAMES[model]`) and the reduced chi-square map. Pixels whose
spectra contain no valid channel, or whose fit fails, are set to `NaN` — this
preserves HEALPix UNSEEN/partial-sky masks.
"""
function QUFitCube(Q::AbstractArray{<:Real,3}, U::AbstractArray{<:Real,3},
                   nuArray::AbstractVector;
                   model::Symbol = :screen,
                   sigma_q = nothing, sigma_u = nothing,
                   log_progress::Bool = false, kwargs...)

    size(Q) == size(U) ||
        throw(ArgumentError("Q and U cubes must have the same size (got $(size(Q)) and $(size(U)))."))
    size(Q, 3) == length(nuArray) ||
        throw(ArgumentError("The last axis of Q/U ($(size(Q, 3))) must match nuArray ($(length(nuArray)))."))

    nx, ny = size(Q, 1), size(Q, 2)
    np = length(QU_FIT_PARAM_NAMES[model])
    params = fill(NaN, nx, ny, np)
    stderr = fill(NaN, nx, ny, np)
    chi2_red = fill(NaN, nx, ny)

    log_progress && @info "Starting QU fitting" model npix = nx * ny nchan = length(nuArray)

    for j in 1:ny, i in 1:nx
        qv = vec(Q[i, j, :])
        uv = vec(U[i, j, :])
        any(isfinite, qv) && any(isfinite, uv) || continue
        result = try
            QUFit(qv, uv, nuArray; model=model, sigma_q=sigma_q, sigma_u=sigma_u, kwargs...)
        catch err
            err isa InterruptException && rethrow()
            continue
        end
        params[i, j, :] .= result.params
        stderr[i, j, :] .= result.stderr
        chi2_red[i, j] = result.chi2_red
    end

    log_progress && @info "QU fitting finished" model
    return params, stderr, chi2_red
end
