"""
    make_demo_data(dir::AbstractString; kwargs...) -> NamedTuple

Generate a tiny, fully synthetic demo dataset with **analytically known
results**, ready to be processed by `MOOSE_from_config`. This provides a
quickstart that exercises the whole pipeline without any real simulation
data, and lets you check the outputs against closed-form expectations.

# Geometry: a Faraday screen in front of a uniform emitter

Along the line of sight (LOS = z, cube axis 3), the box is split in two:
- slices `1:h` (front half, `h = npix ÷ 2`) form a **Faraday screen**:
  uniform `B_los`, zero plane-of-sky field, hence Faraday rotation but no
  synchrotron emission;
- slices `h+1:npix` (back half) form a **uniform emitter**: plane-of-sky
  field `B_perp` along axis 1, zero LOS field, hence emission but no
  additional rotation.

Density and ionization fraction are uniform, so every emitting cell sees the
same foreground rotation measure `RM_s`. The emission is therefore
Faraday-thin: the FDF is a single peak at `φ = RM_s`, and Q/U follow the
textbook rotation `χ(λ²) = ψ_src + RM_s·λ²` exactly.

The emissivity table is an exact power law in frequency, linear in B, and is
sampled by the pipeline exactly on its `(B, ν)` nodes, so interpolation
introduces no error: every expected value below is exact to floating-point
precision, except the FDF peak position which is quantized on the `dphi`
grid.

# Arguments
- `dir`: directory in which the demo dataset is created (created if needed).

# Keywords
- `npix = 16`: cube size (the cube is `npix × npix × npix`, `npix ≥ 2`).
- `box_length_pc = 10.0`: physical box length in parsec.
- `B_perp_uG = 3.0`: plane-of-sky field of the emitting half (µG). Must be a
  node of `B_nodes_uG`.
- `B_los_uG = 1.0`: LOS field of the Faraday screen (µG).
- `density_cm3 = 1.0`: gas density (cm⁻³), uniform.
- `ionization_fraction = 0.1`: constant ionization fraction (`ne_option = 2`).
- `alpha = -0.7`: spectral index of the demo emission, `S_ν ∝ ν^alpha`.
- `pol_fraction = 0.75`: intrinsic polarization fraction `(e⊥−e∥)/(e⊥+e∥)`.
- `eps0 = 1.0e-39`: total emissivity at `B = B_perp_uG` and the central
  frequency (erg s⁻¹ cm⁻³ Hz⁻¹ sr⁻¹).
- `nu_MHz = 100.0:25.0:200.0`: frequency channels (MHz), uniformly spaced.
  Also used as the ν nodes of the emissivity table, so the pipeline samples
  the table exactly.
- `B_nodes_uG = 0.0:1.0:6.0`: B nodes of the emissivity table (µG). Must
  contain at least 4 values (cubic spline), include `B_perp_uG`, and start
  at 0 so the screen half emits exactly nothing.
- `faraday = true`: enable Faraday rotation (and RM synthesis) in the config.
- `phimin = -5.0, phimax = 5.0, dphi = 0.25`: Faraday-depth grid (rad/m²).

# Returns
A NamedTuple with `base_dir`, `simulation_dir`, `config_path`,
`emissivity_path`, and `expected`, a NamedTuple of analytic expectations:
- `rm`: uniform value of `RMmap.fits` (rad/m²):
  `RM_s = 0.81 · x_e·n · B_los · L_screen_pc` (exact).
- `alpha`: uniform value of `alpha.fits` (exact).
- `Tnu`: brightness temperature per channel (K), uniform over the map:
  `Tν = c²/(2 k_B ν²) · ε_I(B_perp, ν) · L_emit_cm` (exact).
- `nu_MHz`, `lambda2_m2`: channel frequencies and λ² values.
- `qnu_over_tnu`, `unu_over_tnu`: per-channel `Q/T` and `U/T` (exact):
  `Q/T = p·cos 2(ψ_src + RM λ²)`, `U/T = p·sin 2(ψ_src + RM λ²)` with
  `ψ_src = π` (equivalent to 0 modulo π), `p = pol_fraction`, and `RM = RM_s` when `faraday = true`,
  `RM = 0` otherwise.
- `pol_fraction`: `√(Q²+U²)/T` per channel (exact, no depolarization since
  the emitter is Faraday-thin).
- `intrinsic_pol_angle = 0`: `ψ_src` modulo π (rad); equals `½·atan(U, Q)` mod π
  when `faraday = false`.
- `intne`, `intBLOS`: integrated maps (cm⁻², µG·cm), exact and uniform.
- `fdf_peak_phi`: Faraday depth of the |FDF| peak = `RM_s` (`nothing` when
  `faraday = false`); the measured grid peak lies within `dphi` of it.

To run the demo end to end:
```julia
using Moose
demo = make_demo_data("demo")
MOOSE_from_config(demo.config_path; quiet = true)
demo.expected.rm       # compare with RMmap.fits
demo.expected.alpha    # compare with alpha.fits
```
The expected values are also written to `expected_results.json` next to the
config so they can be inspected without a Julia session.
"""
function make_demo_data(dir::AbstractString;
    npix::Integer = 16,
    box_length_pc::Real = 10.0,
    B_perp_uG::Real = 3.0,
    B_los_uG::Real = 1.0,
    density_cm3::Real = 1.0,
    ionization_fraction::Real = 0.1,
    alpha::Real = -0.7,
    pol_fraction::Real = 0.75,
    eps0::Real = 1.0e-39,
    nu_MHz = 100.0:25.0:200.0,
    B_nodes_uG = 0.0:1.0:6.0,
    faraday::Bool = true,
    phimin::Real = -5.0,
    phimax::Real = 5.0,
    dphi::Real = 0.25,
)
    npix >= 2 || error("npix must be >= 2, got $npix.")
    box_length_pc > 0 || error("box_length_pc must be > 0.")
    eps0 > 0 || error("eps0 must be > 0.")
    0 <= pol_fraction <= 1 || error("pol_fraction must be in [0, 1].")
    0 < ionization_fraction <= 1 || error("ionization_fraction must be in (0, 1].")

    nu = collect(Float64, nu_MHz)
    length(nu) >= 4 || error("At least 4 frequency channels are required (cubic spline nodes), got $(length(nu)).")
    issorted(nu) && allunique(nu) || error("nu_MHz must be strictly increasing.")
    all(>(0), nu) || error("All frequencies must be positive.")
    step_MHz = nu[2] - nu[1]
    all(i -> isapprox(nu[i + 1] - nu[i], step_MHz; rtol = 1e-10), 1:(length(nu) - 1)) ||
        error("nu_MHz must be uniformly spaced (the config schema uses start/end/step).")

    Bnodes = collect(Float64, B_nodes_uG)
    length(Bnodes) >= 4 || error("At least 4 B nodes are required (cubic spline), got $(length(Bnodes)).")
    issorted(Bnodes) && allunique(Bnodes) || error("B_nodes_uG must be strictly increasing.")
    Bnodes[1] == 0.0 || error("B_nodes_uG must start at 0 so the screen half emits exactly nothing.")
    Float64(B_perp_uG) in Bnodes ||
        error("B_perp_uG = $(B_perp_uG) must be exactly one of the emissivity B nodes $(Bnodes), so the demo results stay exact.")

    base_dir = abspath(dir)
    simulation_dir = joinpath(base_dir, "demo_simu")
    mkpath(simulation_dir)

    # --- Input cubes: Faraday screen (front) + uniform emitter (back) ----
    # LOS = z maps (B1, B2, BLOS) to (Bx, By, Bz) with no axis permutation,
    # and the RM cube is a cumulative sum from slice 1 to slice npix. Putting
    # the screen in slices 1:h means every emitting cell sees the full screen
    # RM. The pipeline evaluates IntrinsicAngle(B1, B2), so By = 0 and Bx > 0
    # gives psi_src = π, equivalent to 0 modulo π, in the emitter.
    n = Int(npix)
    h = fld(n, 2)

    bx = zeros(n, n, n)
    bx[:, :, (h + 1):n] .= Float64(B_perp_uG)
    bz = zeros(n, n, n)
    bz[:, :, 1:h] .= Float64(B_los_uG)

    _write_demo_cube(joinpath(simulation_dir, "Bx.fits"), bx)
    _write_demo_cube(joinpath(simulation_dir, "By.fits"), zeros(n, n, n))
    _write_demo_cube(joinpath(simulation_dir, "Bz.fits"), bz)
    _write_demo_cube(joinpath(simulation_dir, "density.fits"), fill(Float64(density_cm3), n, n, n))
    _write_demo_cube(joinpath(simulation_dir, "temperature.fits"), fill(100.0, n, n, n))

    # --- Emissivity table: exact power law ------------------------------
    # eps_I(B, ν) = eps0 · (B / B_perp) · (ν / ν_ref)^alpha, split so that
    # (e_perp − e_para) / (e_perp + e_para) = pol_fraction. The table is
    # linear in B (zero at B = 0) and the pipeline samples it exactly on the
    # (B, ν) nodes, so interpolation is exact.
    nu_ref = nu[cld(length(nu), 2)]
    emissivity_path = joinpath(base_dir, "emissivity.csv")
    open(emissivity_path, "w") do io
        write(io, "B\tnu\te_para\te_perp\n")
        for nui in nu, B in Bnodes
            eps_I = Float64(eps0) * (B / Float64(B_perp_uG)) * (nui / nu_ref)^Float64(alpha)
            e_perp = 0.5 * (1 + Float64(pol_fraction)) * eps_I
            e_para = 0.5 * (1 - Float64(pol_fraction)) * eps_I
            write(io, "$(B)\t$(nui)\t$(e_para)\t$(e_perp)\n")
        end
    end

    # --- Config ----------------------------------------------------------
    config_path = joinpath(base_dir, "demo_config.json")
    cfg = Dict{String, Any}(
        "base_dir" => base_dir,
        "simulations" => [simulation_dir],
        "chosen_LOS" => ["z"],
        "interpolation_file_path" => "emissivity.csv",
        "faraday" => Dict(
            "enabled" => faraday ? "Y" : "N",
            "phimin" => Float64(phimin),
            "phimax" => Float64(phimax),
            "dphi" => Float64(dphi),
        ),
        "responseSynchrotron" => "N",
        "add_noise" => "N",
        "ne_option" => "2",
        "IonizationFraction" => Float64(ionization_fraction),
        "freq" => Dict("start" => nu[1], "end" => nu[end], "step" => step_MHz),
        "BoxLength_pc" => Float64(box_length_pc),
        "BoxLength_pix" => n,
        "log_progress" => false,
        "rng_seed" => 42,
    )
    write(config_path, JSON.json(cfg, 2))

    # --- Analytic expectations -------------------------------------------
    ne = Float64(ionization_fraction) * Float64(density_cm3)
    dl_pc = Float64(box_length_pc) / n
    screen_pc = h * dl_pc
    emit_cm = (n - h) * dl_pc * PARSEC_TO_CM
    rm_screen = RM_PREFACTOR * ne * Float64(B_los_uG) * screen_pc

    Tnu_expected = [
        (C^2 / (2 * K_B * (nui * 1e6)^2)) *
        (Float64(eps0) * (nui / nu_ref)^Float64(alpha)) * emit_cm
        for nui in nu
    ]

    psi_src = pi
    intrinsic_pol_angle = 0.0
    lambda2 = [(C_m / (nui * 1e6))^2 for nui in nu]
    rm_seen = faraday ? rm_screen : 0.0
    qnu_over_tnu = [Float64(pol_fraction) * cos(2 * (psi_src + rm_seen * l2)) for l2 in lambda2]
    unu_over_tnu = [Float64(pol_fraction) * sin(2 * (psi_src + rm_seen * l2)) for l2 in lambda2]

    expected = (;
        rm = rm_screen,
        alpha = Float64(alpha),
        Tnu = Tnu_expected,
        nu_MHz = nu,
        lambda2_m2 = lambda2,
        qnu_over_tnu,
        unu_over_tnu,
        pol_fraction = Float64(pol_fraction),
        intrinsic_pol_angle,
        intne = ne * Float64(box_length_pc) * PARSEC_TO_CM,
        intBLOS = Float64(B_los_uG) * screen_pc * PARSEC_TO_CM,
        fdf_peak_phi = faraday ? rm_screen : nothing,
    )

    expected_path = joinpath(base_dir, "expected_results.json")
    expected_doc = Dict{String, Any}(
        "description" => "Analytic expectations for the MOOSE demo dataset (Faraday screen + uniform power-law emitter).",
        "RMmap_rad_m2" => expected.rm,
        "alpha" => expected.alpha,
        "Tnu_K_per_channel" => expected.Tnu,
        "nu_MHz" => expected.nu_MHz,
        "lambda2_m2" => expected.lambda2_m2,
        "Qnu_over_Tnu_per_channel" => expected.qnu_over_tnu,
        "Unu_over_Tnu_per_channel" => expected.unu_over_tnu,
        "pol_fraction" => expected.pol_fraction,
        "intrinsic_pol_angle_rad" => expected.intrinsic_pol_angle,
        "intne_cm2" => expected.intne,
        "intBLOS_uG_cm" => expected.intBLOS,
        "FDF_peak_phi_rad_m2" => expected.fdf_peak_phi,
        "notes" => [
            "All values are exact to floating-point precision (emitter sampled on the emissivity table nodes).",
            "The FDF peak position is quantized on the Faraday-depth grid: the measured peak lies within dphi of FDF_peak_phi.",
            "Q/T and U/T include the screen rotation chi(lambda^2) = pi + RM * lambda^2 when Faraday rotation is enabled (pi is equivalent to 0 modulo pi).",
        ],
    )
    write(expected_path, JSON.json(expected_doc, 2))

    return (; base_dir, simulation_dir, config_path, emissivity_path, expected)
end

function _write_demo_cube(path::AbstractString, data::AbstractArray)
    FITS(path, "w") do fits
        write(fits, data)
    end
    return path
end
