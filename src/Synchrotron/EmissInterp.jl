
#  Equations from Padovani et. 2021 (https://doi.org/10.1051/0004-6361/202140799) for the interpolation code 

using ProgressMeter
using QuadGK
using SpecialFunctions

nu_c(E,BField) = PRE_NU_C * BField * (E / ELECTRON_ENERGY_AT_REST_eV)^2

freq_norm(E,nu,BField) = nu / (nu_c(E,BField))

F_integrand(x) = besselk(5/3,x)

F(x) = x * quadgk(F_integrand,x,Inf, rtol=1.e-3)[1]

G(x) = x * (besselk(2/3,x))

relativistic_ve(E) = C*(sqrt(1-(1/(1+E/ELECTRON_ENERGY_AT_REST_eV)^2)))

j_e(E,a=-1.3,b=1.9) = J_0 * E^a / (E+E_0)^b #Padovani 2021

je_ve_ratio(E) = j_e(E) / relativistic_ve(E)

power_par(E,nu,BField) = PRE_P * BField * (F(freq_norm(E,nu,BField)) - G(freq_norm(E,nu,BField)))
power_perp(E,nu,BField)= PRE_P * BField * (F(freq_norm(E,nu,BField)) + G(freq_norm(E,nu,BField)))

par_integrand(E,nu,BField) = je_ve_ratio(E) * power_par(E,nu,BField)
perp_integrand(E,nu,BField) = je_ve_ratio(E) * power_perp(E,nu,BField)

function par_emissivity(nu_MHz,BField_microG)
    nu = 1e6 * nu_MHz
    BField = 1e-6 * BField_microG
    return quadgk(x -> par_integrand(x,nu,BField),ELECTRON_ENERGY_AT_REST_eV,EMAX)[1]
end

function perp_emissivity(nu_MHz,BField_microG)
    nu = 1e6 * nu_MHz
    BField = 1e-6 * BField_microG
    return quadgk(x -> perp_integrand(x,nu,BField),ELECTRON_ENERGY_AT_REST_eV,EMAX)[1]
end

# F and G are evaluated once per energy and both emissivities are integrated in a single quadgk
function emissivities_integrand(E,nu,BField)
    x = freq_norm(E,nu,BField)
    Fx, Gx = F(x), G(x)
    w = je_ve_ratio(E) * PRE_P * BField
    return [w * (Fx - Gx), w * (Fx + Gx)]
end

function emissivities(nu_MHz,BField_microG)
    nu = 1e6 * nu_MHz
    BField = 1e-6 * BField_microG
    e_para, e_perp = quadgk(x -> emissivities_integrand(x,nu,BField),ELECTRON_ENERGY_AT_REST_eV,1e10)[1]
    return e_para, e_perp
end

"""
    EmissInterp(BArray::AbstractArray, nuArray::AbstractArray)

Calculate the emissivity for a range of magnetic fields and frequencies, and write the results to a file.

# Arguments
- `BArray::AbstractArray`: Array of magnetic field strengths in microGauss.
- `nuArray::AbstractArray`: Array of frequencies in MHz.

# Returns
- `String`: The path to the file where the results are saved.

# Description
This function calculates the parallel and perpendicular emissivities for each combination of magnetic field strengths and frequencies provided in `BArray` and `nuArray`. The results are written to a file named "emissivity.dat" in a tab-separated format with columns for magnetic field strength (`B`), frequency (`nu`), parallel emissivity (`e_para`), and perpendicular emissivity (`e_perp`).

# Example
```julia
BArray = [1.0, 2.0, 3.0]
nuArray = [100, 200, 300]
EmissInterp(BArray, nuArray)
```
"""
function EmissInterp(BArray::AbstractArray, nuArray::AbstractArray; fname="emissivity.dat")
    p = Progress(length(nuArray) * length(BArray); desc="Émissivités : ", showspeed=true)
    open(fname, "w") do f
        write(f, "B\tnu\te_para\te_perp\n")
        for nui in nuArray, Bi in BArray
            e_para, e_perp = emissivities(nui, Bi)
            write(f, "$Bi\t$nui\t$e_para\t$e_perp\n")
            next!(p)
        end
    end
    println(" saved in: ", abspath(fname))
    return abspath(fname)
end