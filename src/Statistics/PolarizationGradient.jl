"""
Spatial polarization gradient maps (Gaensler et al. 2011, Nature 478, 214).

The magnitude of the spatial gradient of the complex linear polarization
`P = Q + iU`,

    |∇P| = √( (∂Q/∂x)² + (∂Q/∂y)² + (∂U/∂x)² + (∂U/∂y)² ),

traces sharp changes in the magneto-ionic medium (shocks, cusps, turbulent
density/field fluctuations) and, unlike the polarized intensity `|P|`, is
invariant under both a global angle rotation and the addition of a uniform
polarized foreground/background. Filaments of high `|∇P|` are a classic
diagnostic of turbulence in the warm ionized medium.
"""

"""
    polarization_gradient_map(Q, U; pixel_size=1.0, normalized=false)

Compute the polarization gradient magnitude `|∇P|` of Stokes maps `Q` and `U`
(Gaensler et al. 2011).

# Arguments
- `Q::AbstractMatrix`, `U::AbstractMatrix`: Stokes maps of identical size,
  in any (consistent) unit. 3D cubes `(nx, ny, nchan)` are also accepted;
  the gradient is then computed independently for every channel.

# Keywords
- `pixel_size::Real = 1.0`: linear size of a pixel; gradients are returned
  per unit of `pixel_size` (e.g. pass the pixel scale in pc or arcmin to get
  `|∇P|` per pc or per arcmin).
- `normalized::Bool = false`: when `true`, return `|∇P| / |P|` (dimensionless,
  emphasizes angle structure over intensity structure). Pixels where
  `|P| = 0` are set to `NaN`.

Derivatives use central differences in the interior and one-sided differences
on the edges. `NaN` pixels (masked/UNSEEN) propagate only to the derivative
estimates that touch them; a pixel whose finite neighbours allow a one-sided
difference still gets a value.

Returns an array of the same size as `Q`.
"""
function polarization_gradient_map(Q::AbstractMatrix, U::AbstractMatrix;
                                   pixel_size::Real = 1.0,
                                   normalized::Bool = false)
    size(Q) == size(U) ||
        throw(ArgumentError("Q and U must have the same size (got $(size(Q)) and $(size(U)))."))
    pixel_size > 0 ||
        throw(ArgumentError("pixel_size must be positive (got $pixel_size)."))

    nx, ny = size(Q)
    grad = fill(NaN, nx, ny)
    inv_h = 1.0 / float(pixel_size)

    @inbounds for j in 1:ny, i in 1:nx
        # Masked pixels stay masked even when central differences could be
        # formed from their finite neighbours.
        (isfinite(Q[i, j]) && isfinite(U[i, j])) || continue
        dqx = _difference_1d(Q, i, j, 1)
        dqy = _difference_1d(Q, i, j, 2)
        dux = _difference_1d(U, i, j, 1)
        duy = _difference_1d(U, i, j, 2)
        if isfinite(dqx) && isfinite(dqy) && isfinite(dux) && isfinite(duy)
            g = sqrt(dqx^2 + dqy^2 + dux^2 + duy^2) * inv_h
            if normalized
                pnorm = sqrt(Q[i, j]^2 + U[i, j]^2)
                grad[i, j] = pnorm > 0 ? g / pnorm : NaN
            else
                grad[i, j] = g
            end
        end
    end

    return grad
end

function polarization_gradient_map(Q::AbstractArray{<:Real,3}, U::AbstractArray{<:Real,3};
                                   kwargs...)
    size(Q) == size(U) ||
        throw(ArgumentError("Q and U must have the same size (got $(size(Q)) and $(size(U)))."))
    grad = fill(NaN, size(Q))
    for k in 1:size(Q, 3)
        grad[:, :, k] .= polarization_gradient_map(view(Q, :, :, k), view(U, :, :, k); kwargs...)
    end
    return grad
end

# Central difference along dimension `dim` when both neighbours are finite,
# else one-sided, else NaN. Returns the derivative in pixel units.
function _difference_1d(A::AbstractMatrix, i::Int, j::Int, dim::Int)
    n = size(A, dim)
    k = dim == 1 ? i : j
    at(m) = dim == 1 ? A[m, j] : A[i, m]

    prev = k > 1 ? at(k - 1) : NaN
    next = k < n ? at(k + 1) : NaN
    here = at(k)

    if isfinite(prev) && isfinite(next)
        return (next - prev) / 2
    elseif isfinite(next) && isfinite(here)
        return next - here
    elseif isfinite(prev) && isfinite(here)
        return here - prev
    else
        return NaN
    end
end
