function compute_energies_from_simulation(NameSimulation; box_size_pc=50.0, Npix=256, dens_factor=1.67e-24, vel_factor=1e5, B_factor=1e-6)
    B1, B2, BLOS, V1, V2, VLOS, T, n = ReadSimulation(NameSimulation, "x", 1, 1, 1, 1e3)
    return compute_energies(n, B1, B2, BLOS, V1, V2, VLOS;
        box_size_pc=box_size_pc, Npix=Npix, dens_factor=dens_factor, vel_factor=vel_factor, B_factor=B_factor)
end

function compute_energies(n, Bx, By, Bz, Vx, Vy, Vz; box_size_pc=50.0, Npix=256, dens_factor=1.67e-24, vel_factor=1e5, B_factor=1e-6)
    dx = box_size_pc / Npix * 3.086e18
    @printf("dx = %.3e cm\n", dx)

    G = 6.674e-8
    dV = dx^3

    Vx, Vy, Vz = Vx .* vel_factor, Vy .* vel_factor, Vz .* vel_factor
    Bx, By, Bz = Bx .* B_factor, By .* B_factor, Bz .* B_factor
    n_mass = n .* dens_factor

    M = sum(n_mass) * dV
    @printf("Total mass M = %.3e g\n", M)

    nx, ny, nz = size(n_mass)
    cx, cy, cz = (nx+1)/2, (ny+1)/2, (nz+1)/2
    X = (collect(1:nx) .- cx) .* dx
    Y = (collect(1:ny) .- cy) .* dx
    Z = (collect(1:nz) .- cz) .* dx
    R² = [x^2 + y^2 + z^2 for x in X, y in Y, z in Z]
    R = sqrt(sum(n_mass .* R²) * dV / M)

    @printf("R = %.3e cm\n", R)

    E_grav = -(3/5) * G * M^2 / R

    B2 = Bx.^2 .+ By.^2 .+ Bz.^2
    B_energy = B2 ./ (8π)
    E_mag = sum(B_energy) * dV
    std_E_mag = std(B_energy) * dV * sqrt(length(B_energy))

    v₀ = (mean(Vx), mean(Vy), mean(Vz))
    dv2 = @. (Vx - v₀[1])^2 + (Vy - v₀[2])^2 + (Vz - v₀[3])^2
    turb_energy = 0.5 .* n_mass .* dv2
    E_turb = sum(turb_energy) * dV
    std_E_turb = std(turb_energy) * dV * sqrt(length(turb_energy))

    grav_energy_density = -(3/5) * G * (n_mass .* dV).^2 ./ R
    std_E_grav = std(grav_energy_density)

    @printf("|E_grav| = %.3e erg\n", abs(E_grav))
    @printf("E_mag        = %.3e erg\n", E_mag)
    @printf("E_turb       = %.3e erg\n", E_turb)

    return E_grav, E_mag, E_turb, std_E_mag, std_E_turb, std_E_grav
end