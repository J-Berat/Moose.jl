using Test
using Logging
using MOOSE
include(joinpath(@__DIR__, "..", "src", "MOOSE_cli.jl"))
using MOOSE.MOOSEFromConfig: build_config

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
end
