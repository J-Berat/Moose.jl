using StatsBase, Plots, LaTeXStrings
include("/Users/jb270005/Desktop/Julia/Script_obs_synth/RM_synthesis_1D.jl")

"""
    slab_Faraday(Phi_range::StepRangeLen, phi_min::Vector{Float64}, phi_max::Vector{Float64},
                 amplitude::Vector{Float64}, phase::Vector{Float64})

This function calculates the Faraday dispersion function (FDF) and related quantities for a slab model based on a given set of input parameters.

# Arguments
- `Phi_range::StepRangeLen`: The range of Faraday depth values (in radians per square meter) over which to calculate the FDF.
- `phi_min::Vector`: A vector containing the minimum Faraday depths of the slabs.
- `phi_max::Vector`: A vector containing the maximum Faraday depths of the slabs.
- `amplitude::Vector`: A vector containing the amplitudes of the slabs.
- `phase::Vector`: A vector containing the phases (in degrees) of the slabs.

# Output
- `FDF::Vector{Complex{Float64}}`: The calculated Faraday dispersion function.
- `reFDF::Vector{Float64}`: The real part of the Faraday dispersion function.
- `imFDF::Vector{Float64}`: The imaginary part of the Faraday dispersion function.
- `plot_layout::Layout`: A layout object containing two plots: one for the P_lambda spectrum and one for the FDF.

# Details
- It then calculates the P_lambda spectrum and the ideal spectrum for the slab model based on the input parameters.
- The function uses the `rm_synthesis_1d` function to compute the FDF and its real and imaginary parts.
- Two plots are generated: one for the P_lambda spectrum and one for the FDF, with appropriate labels and axis limits.

# Dependencies
- The function relies on the `rm_synthesis_1d` function to perform the RM synthesis calculations.

# Example
```julia
Phi_range = -5.0:0.1:5.0
phi_min = [-2.0, -1.0, 0.0]
phi_max = [-1.0, 0.0, 1.0]
amplitude = [0.8, 0.4, 0.6]
phase = [45.0, 60.0, 30.0]

FDF, reFDF, imFDF, angle, plot_layout = slab_Faraday(Phi_range, phi_min, phi_max, amplitude, phase)
"""

function slab_Faraday(Phi_range::StepRangeLen, phi_min::Vector, phi_max::Vector,
        amplitude::Vector, phase::Vector, observation::Bool=false)
    
    #tests dimensions
    if length(phi_min) != length(phi_max)
        throw(DimensionMismatch("phi_min should have the same dimension as phi_max"))
    elseif length(phi_min) != length(amplitude)
        throw(DimensionMismatch("phi_min should have the same dimension as amplitude"))
    elseif length(phi_min) != length(phase)
        throw(DimensionMismatch("phi_min should have the same dimension as phase"))
    elseif length(phi_max) != length(amplitude)
        throw(DimensionMismatch("phi_max should have the same dimension as amplitude"))
    elseif length(amplitude) != length(phase)
        throw(DimensionMismatch("amplitude should have the same dimension as phase"))
    end
                            
    # Sampling LOFAR data in lambda^2 space
    c = 299792458  #speed of light in m/s
    num_SB = 308   # number of SubBands observed
    chan_per_SB = 8   # channels in a SB
    SB_bandwidth = 1.831054687500e5   # Subband bandwidth (in Hz)
    SB0_frequency = 1.150375360000e8  # Lowest frequency (in Hz)

    freq_array = ((0:(num_SB * chan_per_SB - 1)) .* SB_bandwidth / chan_per_SB) .+ SB0_frequency
    l2_array = (c ./ freq_array).^2
    l2_0 = mean(l2_array)

    #parameters slab [phi_min, phi_max, amplitude, phase(deg)] 
    parms = [phi_min, phi_max, amplitude, phase]
    nslabs = length(phi_min)
    
    #definition P_lambda and ideal spectrum 
    P_lambda = complex(zeros(length(l2_array)))
    idealspectrum = complex(zeros(length(Phi_range)))
    angle = zeros(length(l2_array))
    absP_lambda_noized = zeros(length(l2_array))
    reP_lambda_noized = zeros(length(l2_array))
    imP_lambda_noized = zeros(length(l2_array))
    
    #Calculation of P_lambda
    for i in 1:nslabs
        @inbounds @. P_lambda = P_lambda + parms[3][i] * exp(2im * (pi/180)*parms[4][i])*(exp(2im*parms[2][i]*l2_array) - exp(2im*parms[1][i]*l2_array))/(2im*l2_array)/(parms[2][i]-parms[1][i])
        if observation == true
            @inbounds @. absP_lambda_noized = rand(Normal(abs.(P_lambda), 0.3))
            @inbounds @. reP_lambda_noized = rand(Normal(real.(P_lambda), 0.3))
            @inbounds @. imP_lambda_noized = rand(Normal(imag.(P_lambda), 0.3))
        end
        @inbounds w = (Phi_range .>= parms[1][i]) .& (Phi_range .<= parms[2][i])
        @inbounds idealspectrum[w] = idealspectrum[w] .+ parms[3][i].*exp.(2im .* (pi/180).* (parms[3][i]))./(parms[2][i].-parms[1][i])
    end
    
    #Calculation of the FDF and real/imaginary parts
    FDF = rm_synthesis_1d(P_lambda, l2_array, Phi_range; lambda_mean = l2_0)[1] #abs of FDF
    reFDF = rm_synthesis_1d(P_lambda, l2_array, Phi_range; lambda_mean = l2_0)[2] #real(FDF)
    imFDF = rm_synthesis_1d(P_lambda, l2_array, Phi_range; lambda_mean = l2_0)[3] #imag(FDF)
    
    #Calculation of the angle of polarisation
    @. angle = (1/2 * atan(imag.(P_lambda),real.(P_lambda))).*180/pi
    
    #Plotting P_lambda and F(phi)
    layout = @layout [a; b; c]
    
    p1 = plot(sqrt.(l2_array),abs.(P_lambda),label = L"|P|",)
    plot!(sqrt.(l2_array),real.(P_lambda),label = L"\mathbf{Re}(P)",ls=:dot)
    plot!(sqrt.(l2_array),imag.(P_lambda),label = L"\mathbf{Im}(P)",ls=:dot)
    if observation == true
        scatter!(sqrt.(l2_array),absP_lambda_noized,ms=[1],mc=:blue,ma=0.1,label="")
        scatter!(sqrt.(l2_array),reP_lambda_noized,ls=:dot,ms=[1],mc=:orange,ma=0.1,label="")
        scatter!(sqrt.(l2_array),imP_lambda_noized,ls=:dot,ms=[1],mc=:green,ma=0.1,label="")
    end
    xlims!(findmin(sqrt.(l2_array))[1],findmax(sqrt.(l2_array))[1])
    xlabel!(L"\lambda \, \, [\textbf{m}]")
    ylabel!("Amp [arb.u.]")
    
    p2 = plot(sqrt.(l2_array),angle,label=L"\mathbf{\chi}")
    xlims!(findmin(sqrt.(l2_array))[1],findmax(sqrt.(l2_array))[1])
    xlabel!(L"\lambda \, \, [\textbf{m}]")
    ylabel!(L"\chi [°]")
    
    p3 = plot(Phi_range,FDF,label=L"F(\Phi)")
    plot!(Phi_range,reFDF,label=L"\mathbf{Re}(F)")
    plot!(Phi_range,imFDF,label=L"\mathbf{Im}(F)")
    plot!(Phi_range,abs.(idealspectrum),label=L"F_{th}(\Phi)")
    xlims!(findmin(Phi_range)[1],findmax(Phi_range)[1])
    xlabel!(L"\Phi \,\, [\textbf{rad.m}^{-2}]")
    ylabel!("Amp [arb.u.]")
    
    plot_layout = plot(p1, p2, p3, layout=layout)
    #returns FDF, Re(FDF), Im(FDF)
    return(FDF, reFDF, imFDF, angle, plot_layout)
end