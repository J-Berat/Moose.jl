using StatsBase, Plots, LaTeXStrings, Random, Distributions
include("/Users/jb270005/Desktop/Codes/Julia/Script_obs_synth/RM_synthesis_1D.jl")

"""
    delta_Faraday(Phi_range::StepRangeLen, phi_peak::Vector,
                  amplitude::Vector, phase::Vector, observation::Bool=false, outputs::String="all")

This function simulates and analyzes Faraday dispersion data based on delta function components.

# Arguments
- `Phi_range::StepRangeLen`: The range of Faraday depth values (in radians per square meter) over which to calculate the FDF.
- `phi_peak::Vector`: A vector containing the central Faraday depths of the delta function components.
- `amplitude::Vector`: A vector containing the amplitudes of the delta function components.
- `phase::Vector`: A vector containing the phases (in degrees) of the delta function components.
- `observation::Bool`: A boolean flag to simulate noisy observations. Defaults to `false`.
- `outputs::String`: Specifies the type of output to return ("all", "P", "angle", "FDF"). Defaults to "all".

# Output
- `FDF::Vector{Complex{Float64}}`: The calculated Faraday dispersion function.
- `reFDF::Vector{Float64}`: The real part of the Faraday dispersion function.
- `imFDF::Vector{Float64}`: The imaginary part of the Faraday dispersion function.
- `angle::Vector{Float64}`: The angle of polarization.
- `plot_layout::Layout`: A layout object containing one or more plots based on the specified `outputs`.

# Details
- The function generates a synthetic Faraday dispersion function based on the input delta function parameters.
- If `observation` is `true`, it adds Gaussian noise to the generated data.
- The angle of polarization is calculated from the complex Faraday dispersion function.
- Depending on the value of `outputs`, the function returns one or more of the following: P_lambda spectrum plot, angle of polarization plot, and Faraday dispersion function plot.

# Example
```julia
Phi_range = -5.0:0.1:5.0
phi_peak = [0.0, -2.0, 1.0]
amplitude = [0.8, 0.6, 0.7]
phase = [45.0, 60.0, 30.0]

FDF, reFDF, imFDF, angle, plot_layout = delta_Faraday(Phi_range, phi_peak, amplitude, phase)
"""
# Sampling LOFAR data in lambda^2 space
const c = 299792458  #speed of light in m/s
const num_SB = 308   # number of SubBands observed
const chan_per_SB = 8   # channels in a SB
const SB_bandwidth = 1.831054687500e5   # Subband bandwidth (in Hz)
const SB0_frequency = 1.150375360000e8  # Lowest frequency (in Hz)

function delta_Faraday(Phi_range::StepRangeLen, phi_peak::Vector,
        amplitude::Vector, phase::Vector, observation::Bool=false, outputs::String="all")
    
    #tests dimensions
    if length(phi_peak) != length(amplitude)
        throw(DimensionMismatch("phi_peak should have the same dimension as amplitude"))
    elseif length(phi_peak) != length(phase)
        throw(DimensionMismatch("phi_peak should have the same dimension as phase"))
    elseif length(amplitude) != length(phase)
        throw(DimensionMismatch("amplitude should have the same dimension as phase"))
    end

    freq_array = ((0:(num_SB * chan_per_SB - 1)) .* SB_bandwidth / chan_per_SB) .+ SB0_frequency
    l2_array = (c ./ freq_array).^2
    l2_0 = mean(l2_array)
    
    nPhi = length(phi_peak)
    
    #parameters peak [phi,amplitude,phase(deg)] 
    parms = [phi_peak, amplitude, phase]
    
    #definition P_lambda, ideal spectrum and angle of polarisation
    P_lambda = complex(zeros(length(l2_array)))
    absP_lambda_noized = zeros(length(l2_array))
    reP_lambda_noized = zeros(length(l2_array))
    imP_lambda_noized = zeros(length(l2_array))
    idealspectrum = complex(zeros(length(Phi_range)))
    angle = zeros(length(l2_array))
    
    #Calculation of P_lambda
    for i in 1:nPhi
        @inbounds @. P_lambda = P_lambda + (parms[2][i] * exp(2im * (parms[1][i] * l2_array + (pi/180)*parms[3][i])))
        if observation == true
            @inbounds @. absP_lambda_noized = rand(Normal(abs.(P_lambda), 0.3))
            @inbounds @. reP_lambda_noized = rand(Normal(real.(P_lambda), 0.3))
            @inbounds @. imP_lambda_noized = rand(Normal(imag.(P_lambda), 0.3))
        end
        @inbounds w = findmin(abs.(Phi_range .- parms[1][i]))[2]
        @inbounds idealspectrum[w] = idealspectrum[w] .+ parms[2][i] .* exp.(2im .* (pi/180).*(parms[3][i]))
    end
    
    #Calculation of the FDF and real/imaginary parts
    FDF = rm_synthesis_1d(P_lambda, l2_array, Phi_range; lambda_mean = l2_0)[1] #abs of FDF
    reFDF = rm_synthesis_1d(P_lambda, l2_array, Phi_range; lambda_mean = l2_0)[2] #real(FDF)
    imFDF = rm_synthesis_1d(P_lambda, l2_array, Phi_range; lambda_mean = l2_0)[3] #imag(FDF)
    
    #Calculation of the angle of polarisation
    angle = (1/2 .* atan.(imag.(P_lambda),real.(P_lambda))) .* 180/pi
    
    #Plotting P_lambda and F(phi)
    layout = @layout [a; b; c; d]
    
    p1 = plot(sqrt.(l2_array),abs.(P_lambda),label = L"|P|",(size = (600, 400)))
    plot!(sqrt.(l2_array),real.(P_lambda),label = L"\mathbf{Re}(P)")
    plot!(sqrt.(l2_array),imag.(P_lambda),label = L"\mathbf{Im}(P)")
        if observation == true
            scatter!(sqrt.(l2_array),absP_lambda_noized,ms=[1],mc=:blue,ma=0.1,label="")
            scatter!(sqrt.(l2_array),reP_lambda_noized,ls=:dot,ms=[1],mc=:orange,ma=0.1,label="")
            scatter!(sqrt.(l2_array),imP_lambda_noized,ls=:dot,ms=[1],mc=:green,ma=0.1,label="")
        end
    xlims!(findmin(sqrt.(l2_array))[1],findmax(sqrt.(l2_array))[1])
    xlabel!(L"\lambda \, \, [\textbf{m}]")
    ylabel!("Amp[arb.u.]")
    
    p2 = plot(sqrt.(l2_array),angle,label=L"\mathbf{\chi}",legend=false,(size = (600, 400)))
    xlims!(findmin(sqrt.(l2_array))[1],findmax(sqrt.(l2_array))[1])
    xlabel!(L"\lambda \, \, [\textbf{m}]")
    ylabel!(L"\chi [°]")
    
    p3 = plot(Phi_range,FDF,label=L"F(\Phi)",(size = (600, 400)))
    plot!(Phi_range,reFDF,ls=:dot,label=L"\mathbf{Re}(F)")
    plot!(Phi_range,imFDF,ls=:dot,label=L"\mathbf{Im}(F)")
    plot!(Phi_range,abs.(idealspectrum),label=L"F_{th}(\Phi)")
    xlims!(findmin(Phi_range)[1],findmax(Phi_range)[1])
    xlabel!(L"\Phi \,\, [\textbf{rad.m}^{-2}]")
    ylabel!("Amp[arb.u.]")
    
    if outputs == "all"
        plot_layout = plot(p1, p2, p3, layout=layout)
    elseif outputs == "P"
        plot_layout = plot(p1)
    elseif outputs == "angle"
        plot_layout = plot(p2)
    elseif outputs == "FDF"
        plot_layout = plot(p3)
    end

    return(FDF, reFDF, imFDF, angle, plot_layout)
end