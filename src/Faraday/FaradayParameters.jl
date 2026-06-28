function rmsynthesis_parameters(ν_min, ν_max, N)
    @debug "Input frequencies must be provided in Hz."

    λ²_min = (C_m / ν_max)^2
    λ²_max = (C_m / ν_min)^2
    Δλ² = (λ²_max - λ²_min)
    δν = (ν_max - ν_min) / N
    δλ² = abs((2 * C_m^2) / ν_min^3) * δν

    δϕ = 2 * sqrt(3) / Δλ²
    max_scale = pi / λ²_min
    ϕ_max = sqrt(3) / δλ²

    return δϕ, max_scale, ϕ_max
end

function rmsynthesis_parameters(ν_range)
    ν_min, ν_max = extrema(float.(ν_range))
    N = length(ν_range)
    return rmsynthesis_parameters(ν_min, ν_max, N)
end
