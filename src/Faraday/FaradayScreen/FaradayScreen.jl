"""
    faraday_screen(Q, U, RM, rangenu)

Perform Faraday screen rotation on Stokes Q and U parameters.

# Arguments
- `Q`: Stokes Q parameter.
- `U`: Stokes U parameter.
- `RM`: Rotation Measure.
- `rangenu`: Array of frequencies.

# Returns
- Arrays `Qrot` and `Urot` representing the Faraday screen rotated Stokes Q and U parameters.
"""

include("/Users/jb270005/Desktop/Codes/Julia/PhysicalParameters/PolarizationAngle.jl")
include("/Users/jb270005/Desktop/Codes/Julia/PhysicalParameters/Pnu.jl")

const c = 3e8 # m.s^-1

function faraday_screen(Q,U,RM,rangenu)
    
    Nfreq = length(rangenu)
    
    Prot = Pnu(Q,U)
    psi = polarization_angle(Q,U) 

    Qrot = zeros(size(Q,1),size(Q,2),Nfreq)
    Urot = zeros(size(U,1),size(U,2),Nfreq)

    for i in 1:Nfreq
        nui = rangenu[i]
        arg = @. 2 * (psi[:,:,i] + RM * (c / nui)^2)
        cos_arg = cos(arg)
        sin_arg = sin(arg)

        Qrot[:,:,i] .= Prot[:,:,i] .* cos_arg
        Urot[:,:,i] .= Prot[:,:,i] .* sin_arg
    end
end