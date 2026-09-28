@testset "Aperture filter and pipeline selection" begin
    args = (; Δx=1.0, Lcut_small=2.0, Llarge=8.0, fNy=0.5)
    hard, _ = Moose.instrument_bandpass(16, 12; args...)
    @test hard == first(Moose.instrument_bandpass_L(16, 12; args...))
    smooth, _ = Moose.instrument_bandpass(16, 12; edge=:aperture, kD=0.08, args...)
    @test all(x -> 0 <= x <= 1, smooth)
    @test any(x -> 0 < x < 1, smooth)
    @test smooth[1, 1] == 0
    sharp, _ = Moose.instrument_bandpass_aperture(16, 12; Δx=1.0,
        Lcut_small=4.0, Llarge=8.0, fNy=0.5, kD=0.02, smooth_high=false, nlut=16)
    @test sharp[5, 1] == 1 # exactly at the sharp upper cut
    @test sharp[6, 1] == 0
    empty_mask, _ = Moose.instrument_bandpass_aperture(4, 4; Δx=1.0,
        Lcut_small=2.0, Llarge=3.0, fNy=0.1, kD=0.02)
    @test all(iszero, empty_mask)
    @test all(isfinite, first(Moose.instrument_bandpass_aperture(1, 1; kD=0.08, args...)))
    @test_throws ErrorException Moose.instrument_bandpass_aperture(4, 4; kD=0.08, nlut=1, args...)
    @test_throws ErrorException Moose.aperture_edge_profile([0.1], 0.1, 0.5, 0.05; nquad=0)
    @test Moose.aperture_edge_profile(Float32[0.1], 0.1f0, 0.5f0, 0.05f0)[1] > 0
    @test Moose.aperture_kD(144e6; d=100) ≈ 30.75 * 144e6 / (2.99792458e8 * 100)
    @test_throws ErrorException Moose.aperture_kD(144e6; d=0)
    @test_throws Moose.MooseError Moose.normalize_filter_options(Dict("edge"=>"aperture"))
    @test_throws Moose.MooseError Moose.normalize_filter_options(Dict("edge"=>"typo"))
    @test_throws Moose.MooseError Moose.normalize_filter_options(Dict("nlut"=>1))
    options = Dict("edge"=>"aperture", "kD"=>0.08, "chromatic"=>true,
                   "reference_frequency_mhz"=>144.0)
    original = repeat(reshape(sin.(1:192), 16, 12), 1, 1, 2)
    q, u, t = copy(original), 2original, 3original
    Moose._apply_synchrotron_filter!(q, u, t, 8; filter_options=options, frequencies_mhz=[72., 144.])
    @test u ≈ 2q
    @test t ≈ 3q
    @test q[:, :, 2] ≈ Moose.apply_instrument_2d(original[:, :, 2], smooth)
    @test !(q[:, :, 1] ≈ q[:, :, 2])
    q, u, t = copy(original), copy(original), copy(original)
    Moose._apply_synchrotron_filter!(q, u, t, 8)
    @test q ≈ Moose.apply_to_array_xy(original, hard)
    # Channel-first support of the public filter API.
    @test Moose.apply_to_array_xy(permutedims(original, (3,1,2)), _ -> hard; n=16, m=12) ≈ permutedims(q, (3,1,2))
    mktempdir() do dir
        demo = Moose.make_demo_data(dir; npix=4)
        cfg = JSON.parsefile(demo.config_path)
        cfg["responseSynchrotron"] = "Y"
        cfg["kernel_size_synchrotron"] = 8.0
        cfg["filter"] = options
        cfg["outputs"] = ["stokes"]
        run_cfg, _ = build_config(cfg, demo.config_path)
        @test run_cfg.filter_options["edge"] == "aperture"
        write(demo.config_path, JSON.json(cfg))
        Moose.MOOSE_from_config(demo.config_path; quiet=true)
        root = joinpath(demo.simulation_dir, "z", "Synchrotron", "WithFaraday", "filtered")
        @test isdir(root)
        @test isfile(joinpath(root, "Qnu.fits"))
        FITS(joinpath(root, "Qnu.fits")) do f
            @test read_header(f[1])["FILTEDGE"] == "aperture"
            @test read_header(f[1])["FILTCHR"] == true
            @test all(isfinite, read(f[1]))
        end
        saved = JSON.parsefile(demo.config_path)
        @test saved["filter"]["chromatic"] == true
    end
end

@testset "rfft filtering matches the complex FFT and reuses plans" begin
    for (n, m) in ((16, 12), (15, 11))
        img = reshape(cos.(1:n*m) .+ 0.3 .* (1:n*m) ./ (n*m), n, m)
        H = first(Moose.instrument_bandpass_aperture(n, m; Δx=1.0,
            Lcut_small=2.0, Llarge=6.0, fNy=0.5, kD=0.05))
        reference = real.(Moose.FFTW.ifft(Moose.FFTW.fft(img) .* H))
        @test Moose.apply_instrument_2d(img, H) ≈ reference atol=1e-12
        # One workspace serves several images, in place.
        ws = Moose.FilterWorkspace(Float64, n, m)
        a, b = copy(img), 2 .* img
        Moose.apply_instrument_2d!(a, a, H, ws)
        Moose.apply_instrument_2d!(b, b, H, ws)
        @test a ≈ reference atol=1e-12
        @test b ≈ 2 .* reference atol=1e-12
    end
    img32 = rand(Float32, 8, 8)
    H = first(Moose.instrument_bandpass_L(8, 8; Δx=1.0, Lcut_small=2.0, Llarge=4.0, fNy=0.5))
    @test eltype(Moose.apply_instrument_2d(img32, H)) == Float32
    asymmetric = ones(8, 8); asymmetric[2, 1] = 0
    @test_throws ArgumentError Moose.apply_instrument_2d(rand(8, 8), asymmetric)
    @test_throws ErrorException Moose.apply_instrument_2d!(zeros(8, 8), rand(8, 8), H,
        Moose.FilterWorkspace(Float64, 4, 4))
end
