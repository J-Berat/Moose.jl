gauss(x, μ, σ) = @. 1 / (σ * sqrt(2π))* exp(-0.5 * ((x - μ) / σ)^2)
gauss(x, Amp, μ, σ) = @. Amp / (σ * sqrt(2π))* exp(-0.5 * ((x - μ) / σ)^2)