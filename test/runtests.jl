using Test
using MOOSE
using FITSIO
using JSON
using CSV
include(joinpath(@__DIR__, "..", "src", "MOOSE_cli.jl"))
using MOOSE.MOOSEFromConfig: build_config

function write_test_fits(path, data)
    FITS(path, "w") do fits
        write(fits, data)
    end
end

function write_test_emissivity(path)
    open(path, "w") do io
        write(io, "B\tnu\te_para\te_perp\n")
        for nu in (99.0, 100.0, 101.0, 102.0)
            for B in (0.0, 1.0, 2.0, 3.0)
                write(io, "$(B)\t$(nu)\t$(1.0e-40 * (1 + B))\t$(2.0e-40 * (1 + B + nu / 1000))\n")
            end
        end
    end
end

@testset "RMS" begin
    x = [1.0, 2.0, 3.0]
    @test isapprox(MOOSE.RMS(x), sqrt(2 / 3))

    x2 = [1.0, 2.0]
    y2 = [3.0, 4.0]
    @test isapprox(MOOSE.RMS(x2, y2), sqrt(0.5))

    z2 = [5.0, 6.0]
    @test isapprox(MOOSE.RMS(x2, y2, z2), sqrt(0.75))
end

@testset "Pnu" begin
    q = [1.0, 2.0]
    u = [3.0, 4.0]
    expected = sqrt.(q .^ 2 .+ u .^ 2)
    @test MOOSE.Pnu(q, u) == expected
end

@testset "Conversion Jy/beam ↔ K" begin
    intensity = [0.0, 1.5]
    nu = 1.0e9
    theta = 10.0
    expected = 1.222e3 .* intensity ./ (nu^2 * theta^2)
    @test all(isapprox.(MOOSE.ConversionJyK(intensity, nothing, nu, theta), expected; rtol = 1e-12))
    @test iszero(MOOSE.ConversionJyK(0.0, nothing, nu, theta))
end

@testset "Rotation Measure" begin
    BLOS = [1.0, -2.0, 0.5]
    ne = [0.5, 1.0, 0.0]
    pixel_length = 2.0
    expected_delta = 0.81 .* ne .* BLOS .* pixel_length
    @test MOOSE.deltaRM(BLOS, ne, pixel_length) ≈ expected_delta atol = 0 rtol = 1e-12

    @test MOOSE.RM([1.0, 2.0, 3.0]) == [1.0, 3.0, 6.0]

    cube = reshape(1.0:8.0, 2, 2, 2)
    rm_cube = MOOSE.RM(cube)
    @test rm_cube[:, :, 1] == cube[:, :, 1]
    @test rm_cube[:, :, 2] == cube[:, :, 1] .+ cube[:, :, 2]
end

@testset "Moments" begin
    y = [1.0, 2.0, 3.0]
    x = [1.0, 2.0, 3.0]
    m0, m1, m2 = MOOSE.moments(y, x = x)
    @test isapprox(m0, 6.0; atol = 1e-12)
    @test isapprox(m1, 14 / 6; atol = 1e-12)
    expected_m2 = sqrt(sum(y .* (x .- m1) .^ 2) / m0)
    @test isapprox(m2, expected_m2; atol = 1e-12)

    empty_moments = MOOSE.moments([0.0, 0.0]; x = [1.0, 2.0], threshold = 0.5)
    @test all(isnan, empty_moments)
end

@testset "Power spectrum" begin
    delta_field = [1.0 0.0; 0.0 0.0]

    kx, ky, psd2d = MOOSE.power_spectrum_2d(delta_field; detrend_mean = false, normalize = true)
    @test kx == [-0.5, 0.0]
    @test ky == [-0.5, 0.0]
    @test all(isapprox.(psd2d, 0.25; atol = 1e-12))

    k, psd1d = MOOSE.radial_psd(delta_field; detrend_mean = false, normalize = true, nbins = 1)
    @test length(k) == 1
    @test isapprox(psd1d[1], 0.25; atol = 1e-12)
end

@testset "Interferometric filtering" begin
    H, Hshift = MOOSE.instrument_bandpass_L(8, 8; Δx = 1.0, Δy = 1.0, Lcut_small = 1.0, Llarge = 4.0, fNy = 0.5)
    @test size(H) == (8, 8)
    @test size(Hshift) == (8, 8)
    @test Set(vec(H)) ⊆ Set(Float32[0, 1])
    @test H[1, 1] == 0f0
    @test any(==(1f0), H)

    img = reshape(Float64.(1:64), 8, 8)
    out = MOOSE.apply_instrument_2d(img, H)
    @test size(out) == size(img)
    @test eltype(out) <: Real

    cube = repeat(img, 1, 1, 3)
    filtered = MOOSE.apply_to_array_xy(cube, H; n = 8, m = 8)
    @test size(filtered) == size(cube)
    @test filtered[:, :, 1] ≈ out
end

@testset "Regression — LOS basis is cyclic for all three LOS (BUG-1)" begin
    @test MOOSE.los_basis(:Ax, :Ay, :Az, "z") == (:Ax, :Ay, :Az)
    @test MOOSE.los_basis(:Ax, :Ay, :Az, "x") == (:Ay, :Az, :Ax)
    @test MOOSE.los_basis(:Ax, :Ay, :Az, "y") == (:Az, :Ax, :Ay)
    @test_throws ErrorException MOOSE.los_basis(1, 2, 3, "w")
end

@testset "Regression — noise σ derives from SNR (BUG-2)" begin
    rng = MOOSE.Random.MersenneTwister(1234)
    Q = fill(3.0, 32, 32, 2)
    U = fill(4.0, 32, 32, 2)
    Q0, U0 = copy(Q), copy(U)
    snr = 10.0
    MOOSE._add_noise!(Q, U, snr, rng)
    # P_rms = sqrt(3² + 4²) = 5 per channel → expected σ = 5 / 10 = 0.5
    sigmaQ = sqrt(sum(abs2, Q .- Q0) / length(Q))
    sigmaU = sqrt(sum(abs2, U .- U0) / length(U))
    @test isapprox(sigmaQ, 0.5; rtol = 0.15)
    @test isapprox(sigmaU, 0.5; rtol = 0.15)
    @test_throws ErrorException MOOSE._add_noise!(Q, U, -1.0, rng)
end

@testset "Regression — band-pass filter selects requested scales (BUG-3)" begin
    n = 64
    # Pure sine of wavelength 8 pixels along axis 1 (exact FFT bin).
    img = [sin(2π * i / 8) for i in 1:n, j in 1:n]
    Hkeep, _ = MOOSE.instrument_bandpass_L(n, n; Δx = 1.0, Δy = 1.0,
                                           Lcut_small = 2.0, Llarge = 16.0, fNy = 0.5)
    Hkill, _ = MOOSE.instrument_bandpass_L(n, n; Δx = 1.0, Δy = 1.0,
                                           Lcut_small = 2.0, Llarge = 4.0, fNy = 0.5)
    out_keep = MOOSE.apply_instrument_2d(img, Hkeep)
    out_kill = MOOSE.apply_instrument_2d(img, Hkill)
    @test out_keep ≈ img atol = 1e-8       # 2 ≤ 8 ≤ 16 → survives
    @test maximum(abs, out_kill) < 1e-8    # 8 > 4 → annihilated
end

@testset "Regression — emissivity grid handles rectangular/shuffled tables (BUG-4)" begin
    Bvals = [0.0, 1.0, 2.0, 3.0, 4.0]
    nuvals = [10.0, 20.0, 30.0, 40.0]
    rows = [(b, ν) for ν in nuvals for b in Bvals]
    reverse!(rows)  # scramble row order: result must not depend on it
    df = MOOSE.DataFrame(
        B = [r[1] for r in rows],
        nu = [r[2] for r in rows],
        e_perp = [r[1] + 1000 * r[2] for r in rows],
        e_para = zeros(length(rows)),
    )
    B, nu, eps = MOOSE.emissivity_grid(df, df.e_perp .- df.e_para)
    @test size(eps) == (length(Bvals), length(nuvals))
    @test all(eps[i, j] == B[i] + 1000 * nu[j] for i in eachindex(B), j in eachindex(nu))

    df_incomplete = df[1:end-1, :]
    @test_throws ErrorException MOOSE.emissivity_grid(df_incomplete, df_incomplete.e_perp)
end

@testset "Regression — header cache is thread-safe (BUG-5)" begin
    empty!(MOOSE._HEADER_PARAMS_CACHE)
    tasks = [Threads.@spawn MOOSE._header_params_cached("Qnu") for _ in 1:8]
    results = fetch.(tasks)
    @test all(r -> r["bunit"] == "K" && r["cunit3"] == "Hz", results)
end

@testset "Regression — FITS header units (BUG-6)" begin
    @test MOOSE.DictHeader["Pnumax"]["bunit"] == "K"
    @test MOOSE.DictHeader["intBLOS"]["bunit"] == "muG cm"
    @test MOOSE.DictHeader["intBtotal"]["bunit"] == "muG cm"
    @test MOOSE.DictHeader["intBperp"]["bunit"] == "muG cm"
    @test MOOSE.DictHeader["ne"]["ctype3"] == "DIST"
    # Single-channel spectral axis must not throw (CDELT3 falls back to 0.0).
    h = MOOSE.buildHeader3D(3, (2, 2, 1), "", "", "FREQ", "", "", "Hz", "K", [1.0e8])
    @test h["CDELT3"] == 0.0
    @test h["CRVAL3"] == 1.0e8

    h_meta = MOOSE.buildHeader3D(
        3,
        (2, 2, 1),
        "",
        "",
        "FREQ",
        "",
        "",
        "Hz",
        "K",
        [1.0e8];
        metadata = Dict("MOOSEV" => "1.0.0", "CFGHASH" => repeat("a", 64), "SKIPME" => nothing),
    )
    @test h_meta["MOOSEV"] == "1.0.0"
    @test h_meta["CFGHASH"] == repeat("a", 64)
    @test_throws KeyError h_meta["SKIPME"]
end

@testset "Unit constants" begin
    @test MOOSE.C_m == 2.99792458e8
    @test MOOSE.C == 2.99792458e10
end

@testset "Regression — FITS cube reading fails early with actionable validation" begin
    mktempdir() do dir
        two_d_path = joinpath(dir, "map.fits")
        write_test_fits(two_d_path, ones(2, 2))
        err = try
            MOOSE.read_file(two_d_path, 1.0; expected_ndims = 3)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("expected 3D", sprint(showerror, err))

        nan_path = joinpath(dir, "nan_cube.fits")
        cube = ones(2, 2, 2)
        cube[1, 1, 1] = NaN
        write_test_fits(nan_path, cube)
        err = try
            MOOSE.read_file(nan_path, 1.0; expected_ndims = 3)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("non-finite", sprint(showerror, err))
    end
end

@testset "Regression — Tnu3D cached path preserves column values" begin
    mktempdir() do dir
        emissivity_path = joinpath(dir, "emissivity.csv")
        write_test_emissivity(emissivity_path)
        df = CSV.File(emissivity_path) |> MOOSE.DataFrame
        Bperpcube = reshape(collect(range(0.1, 2.5, length = 12)), 2, 3, 2)
        nuArray = [99.0, 100.0, 101.0]
        pixel_length_cm = 2.0

        interpolator = MOOSE.TemperatureInterpolator(df)
        emissivity_cache = MOOSE.build_emissivity_frequency_cache(interpolator, nuArray)
        t3d = MOOSE.Tnu3D(Bperpcube, nuArray, df, pixel_length_cm)

        for i in axes(Bperpcube, 1), j in axes(Bperpcube, 2)
            expected = MOOSE.Tnu(@view(Bperpcube[i, j, :]), nuArray, df, pixel_length_cm;
                precomputed_interp = interpolator,
                emissivity_cache = emissivity_cache,
            )
            @test t3d[i, j, :] ≈ expected
        end
    end
end

@testset "CLI argument validation" begin
    err = try
        parse_cli_args(["--faraday", "maybe"])
        nothing
    catch e
        e
    end

    @test err isa MOOSE.MooseError
    @test err.code == :cli_invalid_argument
    @test occursin("expects Y or N", err.message)

    _, _, _, overrides = parse_cli_args(["--faraday", "true", "--noise", "0"])
    @test overrides["FaradayRotation"] == "Y"
    @test overrides["add_noise"] == "N"

    err = try
        parse_cli_args(["--write-back"])
        nothing
    catch e
        e
    end

    @test err isa MOOSE.MooseError
    @test err.code == :cli_invalid_argument
end

@testset "Config validation" begin
    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu1")
        mkdir(sim_dir)

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => joinpath(base_dir, "emissivity.csv"),
            "chosen_LOS" => ["x", "invalid"],
        )

        err = try
            build_config(cfg, "config.json")
            nothing
        catch e
            e
        end

        @test err isa MOOSE.MooseError
        @test occursin("Invalid line(s) of sight", err.message)
    end

    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu2")
        mkdir(sim_dir)

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => joinpath(base_dir, "emissivity.csv"),
            "freq" => Dict("start" => 150.0, "end" => 140.0, "step" => -1.0),
        )

        err = try
            build_config(cfg, "config.json")
            nothing
        catch e
            e
        end

        @test err isa MOOSE.MooseError
        @test err.code == :invalid_frequency
        @test occursin("Frequency step", err.message) || occursin("greater than the start frequency", err.message)
    end

    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu3")
        mkdir(sim_dir)
        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write(interpolation_path, "B\tnu\te_para\te_perp\n1.0\t120.0\t1.0\t2.0\n")

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => interpolation_path,
            "box" => Dict(
                "x" => 30.0,
                "y" => 60.0,
                "z" => 90.0,
                "npix" => 128,
            ),
            "rng_seed" => 1234,
        )

        run_cfg, _ = build_config(cfg, "config.json")
        @test run_cfg.BoxLength_pc.x == 30.0
        @test run_cfg.BoxLength_pc.y == 60.0
        @test run_cfg.BoxLength_pc.z == 90.0
        @test run_cfg.BoxLength_pix == (; x = 128, y = 128, z = 128)
        @test run_cfg.rng_seed == 1234
    end

    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu4")
        mkdir(sim_dir)
        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write(interpolation_path, "B\tnu\te_para\te_perp\n1.0\t120.0\t1.0\t2.0\n")

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => interpolation_path,
            "faraday" => Dict("enabled" => "true", "phimin" => 20.0, "phimax" => -20.0, "dphi" => 0.1),
        )

        err = try
            build_config(cfg, "config.json")
            nothing
        catch e
            e
        end

        @test err isa MOOSE.MooseError
        @test err.code == :invalid_faraday_range
    end
end

@testset "Regression — BoxLength_pix must match FITS cube dimensions" begin
    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "tiny_simu")
        mkdir(sim_dir)

        cube = fill(1.0, 2, 2, 2)
        write_test_fits(joinpath(sim_dir, "Bx.fits"), cube)
        write_test_fits(joinpath(sim_dir, "By.fits"), fill(0.5, 2, 2, 2))
        write_test_fits(joinpath(sim_dir, "Bz.fits"), fill(0.25, 2, 2, 2))
        write_test_fits(joinpath(sim_dir, "density.fits"), cube)
        write_test_fits(joinpath(sim_dir, "temperature.fits"), fill(100.0, 2, 2, 2))

        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write_test_emissivity(interpolation_path)

        config_path = joinpath(base_dir, "moose_config.json")
        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "chosen_LOS" => ["z"],
            "interpolation_file_path" => interpolation_path,
            "faraday" => Dict("enabled" => "N", "phimin" => -10.0, "phimax" => 10.0, "dphi" => 1.0),
            "responseSynchrotron" => "N",
            "add_noise" => "N",
            "ne_option" => "2",
            "IonizationFraction" => 0.1,
            "freq" => Dict("start" => 100.0, "end" => 101.0, "step" => 1.0),
            "BoxLength_pc" => 2.0,
            "BoxLength_pix" => Dict("x" => 3, "y" => 2, "z" => 2),
            "log_progress" => false,
        )
        write(config_path, JSON.json(cfg))

        err = try
            MOOSE.MOOSE_from_config(config_path; quiet = true)
            nothing
        catch e
            e
        end

        @test err isa MOOSE.MooseError
        @test err.code == :cube_shape_mismatch
        @test occursin("BoxLength_pix", err.message)
    end
end

@testset "Minimal config end-to-end FITS run" begin
    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "tiny_simu")
        mkdir(sim_dir)

        cube = fill(1.0, 2, 2, 2)
        write_test_fits(joinpath(sim_dir, "Bx.fits"), cube)
        write_test_fits(joinpath(sim_dir, "By.fits"), fill(0.5, 2, 2, 2))
        write_test_fits(joinpath(sim_dir, "Bz.fits"), fill(0.25, 2, 2, 2))
        write_test_fits(joinpath(sim_dir, "density.fits"), cube)
        write_test_fits(joinpath(sim_dir, "temperature.fits"), fill(100.0, 2, 2, 2))

        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write_test_emissivity(interpolation_path)

        config_path = joinpath(base_dir, "moose_config.json")
        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "chosen_LOS" => ["z"],
            "interpolation_file_path" => interpolation_path,
            "faraday" => Dict("enabled" => "N", "phimin" => -10.0, "phimax" => 10.0, "dphi" => 1.0),
            "responseSynchrotron" => "false",
            "add_noise" => "N",
            "ne_option" => "2",
            "IonizationFraction" => 0.1,
            "freq" => Dict("start" => 100.0, "end" => 101.0, "step" => 1.0),
            "BoxLength_pc" => 2.0,
            "BoxLength_pix" => 2,
            "log_progress" => false,
            "rng_seed" => 7,
        )
        write(config_path, JSON.json(cfg))

        MOOSE.MOOSE_from_config(config_path; quiet = true)

        result_dir = joinpath(sim_dir, "z", "Synchrotron", "noFaraday")
        @test isdir(result_dir)
        @test isfile(joinpath(result_dir, "Qnu.fits"))
        @test isfile(joinpath(result_dir, "Unu.fits"))
        @test isfile(joinpath(result_dir, "Pnu.fits"))
        @test isfile(joinpath(result_dir, "Tnu.fits"))
        @test isfile(joinpath(result_dir, "Pnumax.fits"))
        @test isfile(joinpath(sim_dir, "z", "Synchrotron", "ne.fits"))
        @test isfile(joinpath(base_dir, "MOOSE_summary.log"))

        qnu = read(FITS(joinpath(result_dir, "Qnu.fits"))[1])
        @test size(qnu) == (2, 2, 2)
        @test all(isfinite, qnu)

        q_header = FITS(joinpath(result_dir, "Qnu.fits")) do fits
            read_header(fits[1])
        end
        @test q_header["MOOSEV"] == MOOSE.moose_version()
        @test q_header["CFGHASH"] == MOOSE.moose_config_hash(cfg)
        @test q_header["LOS"] == "z"
        @test q_header["FARADAY"] == "N"
        @test q_header["FILTER"] == "N"
        @test q_header["NOISE"] == "N"
        @test q_header["NUNIT"] == "MHz input; Hz FITS"

        summary = read(joinpath(base_dir, "MOOSE_summary.log"), String)
        @test occursin("Config read: $(config_path)", summary)
        @test occursin("Config effective: $(config_path)", summary)
        @test occursin("Config saved: $(config_path)", summary)
        @test occursin("Config hash: $(MOOSE.moose_config_hash(cfg))", summary)

        config_before_cli = read(config_path, String)
        run_with_config(config_path, true, false, Dict{String, Any}("chosen_LOS" => ["z"]))
        @test read(config_path, String) == config_before_cli

        summary_after_cli = read(joinpath(base_dir, "MOOSE_summary.log"), String)
        @test occursin("Config effective: $(config_path) + CLI overrides", summary_after_cli)
        @test occursin("Config saved: <not written>", summary_after_cli)
    end
end

@testset "CLI reproducibility options" begin
    _, _, _, overrides = parse_cli_args(["--rng-seed", "42"])
    @test overrides["rng_seed"] == 42
end

@testset "HEALPix support" begin
    q_maps = [
        MOOSE.healpix_map(fill(1.0, 12); order=:ring),
        MOOSE.healpix_map(fill(0.5, 12); order=:ring),
    ]
    u_maps = [
        MOOSE.healpix_map(fill(0.0, 12); order=:ring),
        MOOSE.healpix_map(fill(0.25, 12); order=:ring),
    ]

    q_stack = MOOSE.HealpixStack(q_maps)
    @test size(q_stack) == (12, 2)
    @test q_stack.nside == 1
    @test q_stack.order == :ring

    result = MOOSE.RMSynthesisHealpix(q_stack, u_maps, [1.0e9, 1.1e9], [-10.0, 0.0, 10.0])
    @test size(result.fdf) == (12, 3)
    @test result.nside == 1
    @test result.order == :ring
    @test result.phi == [-10.0, 0.0, 10.0]

    mktempdir() do dir
        map_path = joinpath(dir, "q.fits")
        MOOSE.write_healpix_map(map_path, q_stack[:, 1]; nside=q_stack.nside, order=q_stack.order)
        reread = MOOSE.read_healpix_map(map_path)
        @test collect(reread) == q_stack[:, 1]

        paths = MOOSE.write_healpix_rm_result(joinpath(dir, "rm"), result; prefix="test")
        @test length(paths.fdf) == 3
        @test all(isfile, paths.fdf)
        @test length(MOOSE.read_healpix_map(first(paths.fdf))) == 12
    end
end
